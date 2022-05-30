/*
  Transport PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Major transport algorithms.
*/

#include <algorithm>
#include <vector>

#include "Epetra_Vector.h"
#include "Epetra_IntVector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Import.h"
#include "Teuchos_RCP.hpp"

// Amanzi
#include "BCs.hh"
#include "errors.hh"
#include "EvaluatorMultiplicativeReciprocal.hh"
#include "Mesh.hh"
#include "PDE_Accumulation.hh"
#include "PDE_Diffusion.hh"
#include "PDE_DiffusionFactory.hh"
#include "PK_DomainFunctionFactory.hh"
#include "PK_Utils.hh"
#include "ReconstructionCellLinear.hh"
#include "UniqueLocalIndex.hh"
#include "WhetStoneDefs.hh"

// Amanzi::Transport
#include "MultiscaleTransportPorosityFactory.hh"
#include "Transport_PK.hh"
#include "TransportBoundaryFunction_Alquimia.hh"
#include "TransportBoundaryFunction_Chemistry.hh"
#include "TransportDomainFunction.hh"
#include "TransportSourceFunction_Alquimia.hh"

namespace Amanzi {
namespace Transport {

using CV_t = CompositeVector;
using CVS_t = CompositeVectorSpace;

/* ******************************************************************
* New constructor compatible with new MPC framework.
****************************************************************** */
Transport_PK::Transport_PK(Teuchos::ParameterList& pk_tree,
                           const Teuchos::RCP<Teuchos::ParameterList>& glist,
                           const Teuchos::RCP<State>& S,
                           const Teuchos::RCP<TreeVector>& soln) :
  soln_(soln)
{
  S_ = S;

  std::string pk_name = pk_tree.name();
  auto found = pk_name.rfind("->");
  if (found != std::string::npos) pk_name.erase(0, found + 2);

  if (glist->isSublist("cycle driver")) {
    if (glist->sublist("cycle driver").isParameter("component names")) {
      // grab the component names
      component_names_ = glist->sublist("cycle driver")
          .get<Teuchos::Array<std::string> >("component names").toVector();
    } else {
      Errors::Message msg("Transport PK: parameter \"component names\" is missing.");
      Exceptions::amanzi_throw(msg);
    }
  } else {
    Errors::Message msg("Transport PK: sublist \"cycle driver\" is missing.");
    Exceptions::amanzi_throw(msg);
  }

  // Create miscaleneous lists.
  Teuchos::RCP<Teuchos::ParameterList> pk_list = Teuchos::sublist(glist, "PKs", true);
  tp_list_ = Teuchos::sublist(pk_list, pk_name, true);

  preconditioner_list_ = Teuchos::sublist(glist, "preconditioners");
  linear_solver_list_ = Teuchos::sublist(glist, "solvers");
  nonlinear_solver_list_ = Teuchos::sublist(glist, "nonlinear solvers");

  subcycling_ = tp_list_->get<bool>("transport subcycling", true);
   
  // domain name
  domain_ = tp_list_->template get<std::string>("domain name", "domain");

  // initialize io
  Teuchos::RCP<Teuchos::ParameterList> units_list = Teuchos::sublist(glist, "units");
  units_.Init(*units_list);

  vo_ = Teuchos::null;
}


/* ******************************************************************
* Old constructor for unit tests.
****************************************************************** */
Transport_PK::Transport_PK(const Teuchos::RCP<Teuchos::ParameterList>& glist,
                           Teuchos::RCP<State> S, 
                           const std::string& pk_list_name,
                           std::vector<std::string>& component_names) :
    component_names_(component_names)
{
  S_ = S;

  // Create miscaleneous lists.
  Teuchos::RCP<Teuchos::ParameterList> pk_list = Teuchos::sublist(glist, "PKs", true);
  tp_list_ = Teuchos::sublist(pk_list, pk_list_name, true);

  preconditioner_list_ = Teuchos::sublist(glist, "preconditioners");
  linear_solver_list_ = Teuchos::sublist(glist, "solvers");
  nonlinear_solver_list_ = Teuchos::sublist(glist, "nonlinear solvers");

  // domain name
  domain_ = tp_list_->template get<std::string>("domain name", "domain");

  // initialize io
  Teuchos::RCP<Teuchos::ParameterList> units_list = Teuchos::sublist(glist, "units");
  units_.Init(*units_list);

  vo_ = Teuchos::null;
}


/* ******************************************************************
* Setup for Alquimia.
****************************************************************** */
#ifdef ALQUIMIA_ENABLED
void Transport_PK::SetupAlquimia()
{
  if (chem_pk_ == Teuchos::null) return;

  alquimia_pk_ = Teuchos::rcp_dynamic_cast<AmanziChemistry::Alquimia_PK>(chem_pk_);
  chem_engine_ = chem_pk_->chem_engine();

  if (chem_engine_ != Teuchos::null) {
    // Retrieve the component names (primary and secondary) from the chemistry 
    // engine.
    std::vector<std::string> component_names;
    chem_engine_->GetPrimarySpeciesNames(component_names);
    component_names_ = component_names;
    for (int i = 0; i < chem_engine_->NumAqueousComplexes(); ++i) {
      char secondary_name[128];
      snprintf(secondary_name, 127, "secondary_%d", i);
      component_names_.push_back(secondary_name);
    }
  }
}
#endif


/* ******************************************************************
* Define structure of this PK.
****************************************************************** */
void Transport_PK::Setup()
{
  passwd_ = "state";  // owner's password

  mesh_ = S_->GetMesh(domain_);
  dim = mesh_->space_dimension();

  // cross-coupling of PKs
  auto physical_models = Teuchos::sublist(tp_list_, "physical models and assumptions");
  bool abs_perm = physical_models->get<bool>("permeability field is required", false);
  std::string multiscale_model = physical_models->get<std::string>("multiscale model", "single continuum");
  transport_on_manifold_ = physical_models->get<bool>("transport in fractures", false);
  use_transport_porosity_ = physical_models->get<bool>("effective transport porosity", false);
  use_effective_diffusion_ = physical_models->get<bool>("effective transport diffusion", false);
  use_dispersion_ = physical_models->get<bool>("use dispersion solver", true);

  // generate keys here to be available for setup of the base class
  tcc_key_ = Keys::getKey(domain_, "total_component_concentration"); 

  permeability_key_ = Keys::getKey(domain_, "permeability"); 
  porosity_key_ = Keys::getKey(domain_, "porosity"); 
  transport_porosity_key_ = Keys::getKey(domain_, "transport_porosity"); 

  std::string tmp = physical_models->get<std::string>("volumetric flow rate key", "volumetric_flow_rate");
  vol_flowrate_key_ = Keys::getKey(domain_, tmp); 

  tmp = physical_models->get<std::string>("saturation key", "saturation_liquid");
  saturation_liquid_key_ = Keys::getKey(domain_, tmp); 

  water_content_key_ = Keys::getKey(domain_, "water_content"); 
  prev_water_content_key_ = Keys::getKey(domain_, "prev_water_content"); 

  porosity_msp_key_ = Keys::getKey(domain_, "porosity_msp"); 
  water_content_msp_key_ = Keys::getKey(domain_, "water_content_msp"); 
  prev_water_content_msp_key_ = Keys::getKey(domain_, "prev_water_content_msp"); 

  aperture_key_ = Keys::getKey(domain_, "aperture"); 

  // require state fields when Flow PK is off
  S_->Require<double>("const_fluid_density", Tags::DEFAULT, "state");

  if (!S_->HasRecord(permeability_key_) && abs_perm) {
    S_->Require<CV_t, CVS_t>(permeability_key_, Tags::DEFAULT, permeability_key_).
      SetMesh(mesh_)->SetGhosted(true)->SetComponent("cell", AmanziMesh::CELL, dim);
  }
  if (!S_->HasRecord(vol_flowrate_key_)) {
    if (transport_on_manifold_) {
      auto cvs = Operators::CreateNonManifoldCVS(mesh_);
      *S_->Require<CV_t, CVS_t>(vol_flowrate_key_, Tags::DEFAULT, passwd_).SetMesh(mesh_)->SetGhosted(true) = *cvs;
    } else {
      S_->Require<CV_t, CVS_t>(vol_flowrate_key_, Tags::DEFAULT, passwd_)
        .SetMesh(mesh_)->SetGhosted(true)->SetComponent("face", AmanziMesh::FACE, 1);
    }
  }

  if (!S_->HasRecord(saturation_liquid_key_)) {
    S_->Require<CV_t, CVS_t>(saturation_liquid_key_, Tags::DEFAULT, saturation_liquid_key_)
      .SetMesh(mesh_)->SetGhosted(true)->SetComponent("cell", AmanziMesh::CELL, 1);
    S_->RequireEvaluator(saturation_liquid_key_, Tags::DEFAULT);
  }

  if (transport_on_manifold_) {
    S_->Require<CV_t, CVS_t>(aperture_key_, Tags::DEFAULT, aperture_key_)
      .SetMesh(mesh_)->SetGhosted(true)->SetComponent("cell", AmanziMesh::CELL, 1);
    S_->RequireEvaluator(aperture_key_, Tags::DEFAULT);
  }

  // require state fields when Transport PK is on
  if (component_names_.size() == 0) {
    Errors::Message msg;
    msg << "Transport PK: list of solutes is empty.\n";
    Exceptions::amanzi_throw(msg);  
  }

  int ncomponents = component_names_.size();
  if (!S_->HasRecord(tcc_key_)) {
    S_->Require<CV_t, CVS_t>(tcc_key_, Tags::DEFAULT, passwd_, component_names_)
      .SetMesh(mesh_)->SetGhosted(true)->SetComponent("cell", AmanziMesh::CELL, ncomponents);
  }

  // porosity evaluators
  if (!S_->HasRecord(porosity_key_)) {
    S_->Require<CV_t, CVS_t>(porosity_key_, Tags::DEFAULT, porosity_key_)
      .SetMesh(mesh_)->SetGhosted(true)->SetComponent("cell", AmanziMesh::CELL, 1);
    S_->RequireEvaluator(porosity_key_, Tags::DEFAULT);
  }

  if (use_transport_porosity_) {
    if (!S_->HasRecord(transport_porosity_key_)) {
      S_->Require<CV_t, CVS_t>(transport_porosity_key_, Tags::DEFAULT, transport_porosity_key_)
        .SetMesh(mesh_)->SetGhosted(true)->SetComponent("cell", AmanziMesh::CELL, 1);
      S_->RequireEvaluator(transport_porosity_key_, Tags::DEFAULT);
    }
  }

  // -- water content (should we move it Flow PK ???)
  if (!S_->HasRecord(water_content_key_)) {
    S_->Require<CV_t, CVS_t>(water_content_key_, Tags::DEFAULT, water_content_key_)
      .SetMesh(mesh_)->SetGhosted(true)->SetComponent("cell", AmanziMesh::CELL, 1);

    std::vector<std::string> listm({ Keys::getVarName(porosity_key_), Keys::getVarName(saturation_liquid_key_) });
    if (transport_on_manifold_) listm.push_back(Keys::getVarName(aperture_key_));

    Teuchos::ParameterList elist(water_content_key_);
    elist.set<std::string>("my key", water_content_key_)
         .set<Teuchos::Array<std::string> >("multiplicative dependencies", listm)
         .set<std::string>("tag", "");
    auto eval = Teuchos::rcp(new EvaluatorMultiplicativeReciprocal(elist));
    S_->SetEvaluator(water_content_key_, Tags::DEFAULT, eval);
  }
  if (!S_->HasRecord(prev_water_content_key_)) {
    S_->Require<CV_t, CVS_t>(prev_water_content_key_, Tags::DEFAULT, passwd_)
      .SetMesh(mesh_)->SetGhosted(true)->SetComponent("cell", AmanziMesh::CELL, 1);
    S_->GetRecordW(prev_water_content_key_, passwd_).set_io_vis(false);
  }

  // require multiscale fields
  multiscale_porosity_ = false;
  if (multiscale_model == "dual continuum discontinuous matrix") {
    multiscale_porosity_ = true;
    msp_ = CreateMultiscaleTransportPorosityPartition(mesh_, tp_list_);

    if (!S_->HasRecord("total_component_concentration_msp")) {
      S_->Require<CV_t, CVS_t>("total_component_concentration_msp", Tags::DEFAULT, passwd_, component_names_)
        .SetMesh(mesh_)->SetGhosted(false)
        ->SetComponent("cell", AmanziMesh::CELL, ncomponents);

      // Secondary matrix nodes are collected here. We assume same order of species.
      int nnodes, nnodes_tmp = NumberMatrixNodes(msp_);
      mesh_->get_comm()->MaxAll(&nnodes_tmp, &nnodes, 1);
      if (nnodes > 1) {
        S_->Require<CV_t, CVS_t>("total_component_concentration_msp_aux", Tags::DEFAULT, passwd_)
          .SetMesh(mesh_)->SetGhosted(false)
          ->SetComponent("cell", AmanziMesh::CELL, ncomponents * (nnodes - 1));
        S_->GetRecordW("total_component_concentration_msp_aux", passwd_).set_io_vis(false);
      }
    }

    // -- water content for multiscale models (a dual porosity model)
    if (!S_->HasRecord(water_content_msp_key_)) {
      S_->Require<CV_t, CVS_t>(water_content_msp_key_, Tags::DEFAULT, passwd_)
        .SetMesh(mesh_)->SetGhosted(true)->SetComponent("cell", AmanziMesh::CELL, 1);
    }
    if (!S_->HasRecord(prev_water_content_msp_key_)) {
      S_->Require<CV_t, CVS_t>(prev_water_content_msp_key_, Tags::DEFAULT, passwd_)
        .SetMesh(mesh_)->SetGhosted(true)->SetComponent("cell", AmanziMesh::CELL, 1);
      S_->GetRecordW(prev_water_content_msp_key_, passwd_).set_io_vis(false);
    }

    // -- porosity of matrix
    if (!S_->HasRecord(porosity_msp_key_)) {
      S_->Require<CV_t, CVS_t>(porosity_msp_key_, Tags::DEFAULT, passwd_)
        .SetMesh(mesh_)->SetGhosted(true)->SetComponent("cell", AmanziMesh::CELL, 1);
    }
  }

#ifdef ALQUIMIA_ENABLED
  SetupAlquimia();
#endif

  // temporary fields
  S_->Require<CV_t, CVS_t>(tcc_key_, Tags::COPY, passwd_, component_names_);
}


/* ******************************************************************
* Routine processes parameter list. It needs to be called only once
* on each processor.                                                     
****************************************************************** */
void Transport_PK::Initialize()
{
  // Set initial values for transport variables.
  dt_ = dt_debug_ = t_physics_ = 0.0;
  double time = S_->get_time();
  if (time >= 0.0) t_physics_ = time;

  dispersion_preconditioner = "identity";

  internal_tests_ = 0;
  internal_tests_tol_ = TRANSPORT_CONCENTRATION_OVERSHOOT;

  // Create verbosity object.
  Teuchos::ParameterList vlist;
  vlist.sublist("verbose object") = tp_list_->sublist("verbose object");
  vo_ = Teuchos::rcp(new VerboseObject("Transport-" + domain_, vlist)); 

  MyPID = mesh_->get_comm()->MyPID();

  // initialize missed fields
  InitializeFields_();

  // Check input parameters. Due to limited amount of checks, we can do it earlier.
  Policy(S_.ptr());

  ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  ncells_wghost = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);

  nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  nfaces_wghost = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);

  nnodes_wghost = mesh_->num_entities(AmanziMesh::NODE, AmanziMesh::Parallel_type::ALL);

  // extract control parameters
  InitializeAll_();
 
  // pointers to state variables (move in subroutines for consistency)
  S_->Get<CV_t>(vol_flowrate_key_).ScatterMasterToGhosted("face");

  phi = S_->Get<CV_t>(porosity_key_).ViewComponent("cell");
  if (use_transport_porosity_) {
    transport_phi = S_->Get<CV_t>(transport_porosity_key_).ViewComponent("cell");
  } else {
    transport_phi = phi;
  }

  // initialize copy of component
  tcc = S_->GetPtrW<CV_t>(tcc_key_, Tags::DEFAULT, passwd_);
  tcc_tmp = S_->GetPtrW<CV_t>(tcc_key_, Tags::COPY, passwd_);
  *tcc_tmp = *tcc;
  S_->GetRecordW(tcc_key_, Tags::COPY, passwd_).set_initialized();

  // upwind structures
  IdentifyUpwindCells();

  // advection block initialization
  current_component_ = -1;

  const Epetra_Map& cmap_owned = mesh_->cell_map(false);
  wc_subcycle_start = Teuchos::rcp(new Epetra_Vector(cmap_owned));
  wc_subcycle_end = Teuchos::rcp(new Epetra_Vector(cmap_owned));

  // reconstruction initialization
  lifting_ = Teuchos::rcp(new Operators::ReconstructionCellLinear(mesh_));
  limiter_ = Teuchos::rcp(new Operators::LimiterCell(mesh_));
  fct_ = Teuchos::rcp(new Operators::FCT(mesh_, mesh_, limiter_, Teuchos::null));

  // create boundary conditions
  if (tp_list_->isSublist("boundary conditions")) {
    // -- try simple Dirichlet conditions for species
    PK_DomainFunctionFactory<TransportDomainFunction> factory(mesh_, S_);
    Teuchos::ParameterList& clist = tp_list_->sublist("boundary conditions").sublist("concentration");

    for (auto it = clist.begin(); it != clist.end(); ++it) {
      std::string name = it->first;
      if (name == "coupling") {
        Teuchos::ParameterList& bc_list = clist.sublist(name);        
        Teuchos::ParameterList::ConstIterator it1 = bc_list.begin();
        std::string specname = it1->first;
        Teuchos::ParameterList& spec = bc_list.sublist(specname);
        Teuchos::RCP<TransportDomainFunction> 
          bc = factory.Create(spec, "boundary concentration", AmanziMesh::FACE, Kxy);
        
        for (int i = 0; i < component_names_.size(); i++) {
          bc->tcc_names().push_back(component_names_[i]);
          bc->tcc_index().push_back(i);
        }
        bc->set_state(S_);
        bc->set_location("interface");
        bcs_.push_back(bc);
      } else {        
        if (clist.isSublist(name)) {
          Teuchos::ParameterList& bc_list = clist.sublist(name);
          for (auto it1 = bc_list.begin(); it1 != bc_list.end(); ++it1) {
            std::string specname = it1->first;
            Teuchos::ParameterList& spec = bc_list.sublist(specname);
            Teuchos::RCP<TransportDomainFunction> 
              bc = factory.Create(spec, "boundary concentration", AmanziMesh::FACE, Kxy);

            bc->tcc_names().push_back(name);
            bc->tcc_index().push_back(FindComponentNumber(name));
            bc->set_state(S_);
            bcs_.push_back(bc);
          }
        }
      }
    }
    // -- try geochemical Dirichlet conditions for species
    PK_DomainFunctionFactory<TransportBoundaryFunction_Chemistry> factory2(mesh_, S_);
    Teuchos::ParameterList& cvlist = tp_list_->sublist("boundary conditions").sublist("constraints");

    for (auto it = cvlist.begin(); it != cvlist.end(); ++it) {
      std::string name = it->first;
      if (cvlist.isSublist(name)) {
        Teuchos::ParameterList& spec = cvlist.sublist(name);
        spec.set<Teuchos::RCP<AmanziChemistry::Chemistry_PK> >("chemical pk", chem_pk_);

        Teuchos::RCP<TransportBoundaryFunction_Chemistry> 
          bc = factory2.Create(spec, "boundary constraints", AmanziMesh::FACE, Kxy);

        for (int i = 0; i < component_names_.size(); i++) {
          bc->tcc_names().push_back(component_names_[i]);
          bc->tcc_index().push_back(i);
        }
        bc->set_state(S_);
        bcs_.push_back(bc);
      }
    }
#ifdef ALQUIMIA_ENABLED
    // -- try geochemical conditions
    Teuchos::ParameterList& glist = tp_list_->sublist("boundary conditions").sublist("geochemical");

    for (auto it = glist.begin(); it != glist.end(); ++it) {
      std::string specname = it->first;
      Teuchos::ParameterList& spec = glist.sublist(specname);

      Teuchos::RCP<TransportBoundaryFunction_Alquimia> 
          bc = Teuchos::rcp(new TransportBoundaryFunction_Alquimia(spec, mesh_, alquimia_pk_, chem_engine_));

      std::vector<int>& tcc_index = bc->tcc_index();
      std::vector<std::string>& tcc_names = bc->tcc_names();

      for (int i = 0; i < tcc_names.size(); i++) {
        tcc_index.push_back(FindComponentNumber(tcc_names[i]));
      }

      bcs_.push_back(bc);
    }
#endif
  } else {
    if (vo_->getVerbLevel() > Teuchos::VERB_NONE) {
      Teuchos::OSTab tab = vo_->getOSTab();
      *vo_->os() << vo_->color("yellow") << "No BCs were specified." << vo_->reset() << std::endl;
    }
  }

  // -- initialization
  time = t_physics_;
  for (int i = 0; i < bcs_.size(); i++) {
    bcs_[i]->Compute(time, time);
    bcs_[i]->ComputeSubmodel(mesh_, tcc);
  }

  VV_CheckInfluxBC();

  // source term initialization: so far only "concentration" is available.
  if (tp_list_->isSublist("source terms")) {
    PK_DomainFunctionFactory<TransportDomainFunction> factory(mesh_, S_);
    PKUtils_CalculatePermeabilityFactorInWell(S_.ptr(), Kxy);

    Teuchos::ParameterList& clist = tp_list_->sublist("source terms").sublist("concentration");
    for (auto it = clist.begin(); it != clist.end(); ++it) {
      std::string name = it->first;
      if (name == "coupling") {
        Teuchos::ParameterList& src_list = clist.sublist(name);
        Teuchos::ParameterList::ConstIterator it1 = src_list.begin();
        std::string specname = it1->first;
        Teuchos::ParameterList& spec = src_list.sublist(specname);
        Teuchos::RCP<TransportDomainFunction> src = factory.Create(spec, "sink", AmanziMesh::CELL, Kxy);
        
        for (int i = 0; i < component_names_.size(); i++) {
          src->tcc_names().push_back(component_names_[i]);
          src->tcc_index().push_back(i);
        }
          
        src->set_state(S_);
        srcs_.push_back(src);

      } else {
        if (clist.isSublist(name)) {
          Teuchos::ParameterList& src_list = clist.sublist(name);
          for (auto it1 = src_list.begin(); it1 != src_list.end(); ++it1) {
            std::string specname = it1->first;
            Teuchos::ParameterList& spec = src_list.sublist(specname);
            Teuchos::RCP<TransportDomainFunction> src = factory.Create(spec, AmanziMesh::CELL, Kxy);

            src->tcc_names().push_back(name);
            src->tcc_index().push_back(FindComponentNumber(name));

            srcs_.push_back(src);
          }
        }
      }
    }
#ifdef ALQUIMIA_ENABLED
    // -- try geochemical conditions
    Teuchos::ParameterList& glist = tp_list_->sublist("source terms").sublist("geochemical");

    for (auto it = glist.begin(); it != glist.end(); ++it) {
      std::string specname = it->first;
      Teuchos::ParameterList& spec = glist.sublist(specname);

      Teuchos::RCP<TransportSourceFunction_Alquimia> 
          src = Teuchos::rcp(new TransportSourceFunction_Alquimia(spec, mesh_, alquimia_pk_, chem_engine_));

      std::vector<int>& tcc_index = src->tcc_index();
      std::vector<std::string>& tcc_names = src->tcc_names();

      for (int i = 0; i < tcc_names.size(); i++) {
        tcc_index.push_back(FindComponentNumber(tcc_names[i]));
      }

      srcs_.push_back(src);
    }
#endif
  }

  // support of the Henry law
  PrepareAirWaterPartitioning_();

  if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "Number of components: " << tcc->size() << std::endl
               << "cfl=" << cfl_ << " spatial/temporal discretization: " 
               << spatial_disc_order << " " << temporal_disc_order << std::endl
               << "using transport porosity: " << use_transport_porosity_ << std::endl;
    *vo_->os() << vo_->color("green") << "Initialization of PK is complete." 
               << vo_->reset() << std::endl << std::endl;
  }
}


/* ******************************************************************
* Initialized fields not touched by State and other PKs.
****************************************************************** */
void Transport_PK::InitializeFields_()
{
  Teuchos::OSTab tab = vo_->getOSTab();

  // set popular default values when flow PK is off
  InitializeField_(saturation_liquid_key_, Tags::DEFAULT, passwd_, 1.0);
  InitializeField_(aperture_key_, Tags::DEFAULT, passwd_, 1.0);

  // we assume that liquid saturation is 1 if wwater content was uninitialized
  InitializeFieldFromField_(prev_water_content_key_, water_content_key_, false);

  InitializeFieldFromField_(water_content_msp_key_, porosity_msp_key_, false);
  InitializeFieldFromField_(prev_water_content_msp_key_, water_content_msp_key_, false);

  InitializeFieldFromField_("total_component_concentration_msp", tcc_key_, false);
  InitializeField_("total_component_concentration_msp_aux", Tags::DEFAULT, passwd_, 0.0);
}


/* ****************************************************************
* Auxiliary initialization technique.
**************************************************************** */
void Transport_PK::InitializeFieldFromField_(
    const std::string& field0, const std::string& field1, bool call_evaluator)
{
  if (S_->HasRecord(field0)) {
    if (S_->GetRecord(field0).owner() == passwd_) {
      if (!S_->GetRecord(field0).initialized()) {
        if (call_evaluator)
            S_->GetEvaluator(field1).Update(*S_, passwd_);

        const auto& f1 = S_->Get<CV_t>(field1);
        auto& f0 = S_->GetW<CV_t>(field0, passwd_);
        f0 = f1;

        S_->GetRecordW(field0, passwd_).set_initialized();

        if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM)
            *vo_->os() << "initialized " << field0 << " to " << field1 << std::endl;
      }
    }
  }
}


/* *******************************************************************
* Estimation of the time step based on T.Barth (Lecture Notes   
* presented at VKI Lecture Series 1994-05, Theorem 4.2.2.       
* Routine must be called every time we update a flow field.
*
* Warning: Barth calculates influx, we calculate outflux. The methods
* are equivalent for divergence-free flows and guarantee the extrema
* diminishing principle. Outflux takes into account sinks and 
* sources but preserves only positivity of an advected mass.
* ***************************************************************** */
double Transport_PK::StableTimeStep()
{
  IdentifyUpwindCells();

  const auto& wc = *S_->Get<CompositeVector>(water_content_key_, Tags::DEFAULT).ViewComponent("cell");
  const auto& wc_prev = *S_->Get<CompositeVector>(prev_water_content_key_, Tags::DEFAULT).ViewComponent("cell");

  // Accumulate upwinding fluxes.
  std::vector<double> total_outflux(ncells_wghost, 0.0);

  for (int f = 0; f < nfaces_wghost; f++) {
    for (int k = 0; k < upwind_cells_[f].size(); k++) {
      int c = upwind_cells_[f][k];
      if (c >= 0) {
        total_outflux[c] += fabs(upwind_flux_[f][k]);
      }
    }
  }

  // Account for extraction of solute in production wells.
  // We assume one well per cell (FIXME).
  double t_old = S_->intermediate_time();
  for (int m = 0; m < srcs_.size(); m++) {
    if (srcs_[m]->keyword() == "producer") {
      srcs_[m]->Compute(t_old, t_old); 

      for (auto it = srcs_[m]->begin(); it != srcs_[m]->end(); ++it) {
        int c = it->first;
        std::vector<double>& values = it->second; 

        for (int i = 0; i < values.size(); ++i) {
          double value = fabs(values[i]) * mesh_->cell_volume(c);
          total_outflux[c] = std::max(total_outflux[c], value);
        }
      }
    }
    else if (srcs_[m]->keyword() == "sink") {  // FIXME
      srcs_[m]->Compute(t_old, t_old); 

      const auto& to_matrix = srcs_[m]->linear_term();
      for (auto it = to_matrix.begin(); it != to_matrix.end(); ++it) {
        int c = it->first;
        total_outflux[c] += std::fabs(it->second);
      }
    }
  }

  // loop over cells and calculate minimal time step
  double vol, outflux, dt_cell;
  dt_ = dt_cell = TRANSPORT_LARGE_TIME_STEP;
  int cmin_dt = -1;
  for (int c = 0; c < ncells_owned; c++) {
    outflux = total_outflux[c];
    if (outflux) {
      vol = mesh_->cell_volume(c);
      dt_cell = vol * std::min(wc_prev[0][c], wc[0][c]) / outflux;
    }
    if (dt_cell < dt_) {
      dt_ = dt_cell;
      cmin_dt = c;
    }
  }

  // correct time step for high-order schemes
  if (spatial_disc_order == 2) dt_ /= 2;

  // no CFL update forsinks, since their are flow dependent.

  // communicate global time step
  double dt_tmp = dt_;
  wc_prev.Comm().MinAll(&dt_tmp, &dt_, 1);

  // incorporate developers and CFL constraints
  dt_ = std::min(dt_, dt_debug_);
  dt_ *= cfl_;

  // print optional diagnostics using maximum cell id as the filter
  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
    int cmin_dt_unique = (fabs(dt_tmp * cfl_ - dt_) < 1e-6 * dt_) ? cmin_dt : -2;
 
    int cmin_dt_tmp = cmin_dt_unique;
    wc_prev.Comm().MaxAll(&cmin_dt_tmp, &cmin_dt_unique, 1);
    if (cmin_dt == cmin_dt_unique && cmin_dt >= 0) {
      const AmanziGeometry::Point& p = mesh_->cell_centroid(cmin_dt);

      Teuchos::OSTab tab = vo_->getOSTab();
      *vo_->os() << "cell " << cmin_dt << " has smallest dt, (" << p[0] << ", " << p[1];
      if (p.dim() == 3) *vo_->os() << ", " << p[2];
      *vo_->os() << ")" << std::endl;
    }
  }
  return dt_;
}


/* ******************************************************************* 
* Estimate returns last time step unless it is zero.     
******************************************************************* */
double Transport_PK::get_dt()
{
  if (subcycling_) {
    return 1e+99;
  } else {
    StableTimeStep();
    return dt_;
  }
}


/* ******************************************************************* 
* Add implict time discretization of the multiscale porosity model on 
* time sub-interval [t_int1, t_int2]:
*   d(VWC_f)/dt -= G_s, d(VWC_m) = G_s 
*   G_s = G_w C^* + omega_s (C_f - C_m).
******************************************************************* */
void Transport_PK::AddMultiscalePorosity_(
    double t_old, double t_new, double t_int1, double t_int2)
{
  Epetra_MultiVector& tcc_next = *tcc_tmp->ViewComponent("cell");
  auto& tcc_matrix = *S_->GetW<CV_t>("total_component_concentration_msp", passwd_).ViewComponent("cell");

  const auto& wcf_prev = *S_->Get<CV_t>(prev_water_content_key_).ViewComponent("cell");
  const auto& wcf = *S_->Get<CV_t>(water_content_key_).ViewComponent("cell");

  const auto& wcm_prev = *S_->Get<CV_t>(prev_water_content_msp_key_).ViewComponent("cell");
  const auto& wcm = *S_->Get<CV_t>(water_content_msp_key_).ViewComponent("cell");

  // multi-node matrix requires more input data
  const auto& phi_matrix = *S_->Get<CV_t>(porosity_msp_key_).ViewComponent("cell");

  int nnodes(1);
  Teuchos::RCP<Epetra_MultiVector> tcc_matrix_aux;
  if (S_->HasRecord("total_component_concentration_msp_aux")) {
    tcc_matrix_aux = S_->GetW<CV_t>("total_component_concentration_msp_aux", passwd_).ViewComponent("cell");
    nnodes = tcc_matrix_aux->NumVectors() + 1; 
  }
  WhetStone::DenseVector tcc_m(nnodes);

  double flux_liquid, wcm0, wcm1, wcf0, wcf1;
  double dtg, dts, t1, t2, tmp0, tmp1, tfp0, tfp1, a, b;
  std::vector<AmanziMesh::Entity_ID> block;

  dtg = t_new - t_old;
  dts = t_int2 - t_int1;
  t1 = t_int1 - t_old;
  t2 = t_int2 - t_old;

  for (int c = 0; c < ncells_owned; ++c) {
    wcm0 = wcm_prev[0][c];
    wcm1 = wcm[0][c];
    flux_liquid = (wcm1 - wcm0) / dtg;
  
    wcf0 = wcf_prev[0][c];
    wcf1 = wcf[0][c];
  
    a = t1 / dtg;
    b = t2 / dtg;

    tfp0 = a * wcf1 + (1.0 - a) * wcf0;
    tfp1 = b * wcf1 + (1.0 - b) * wcf0;

    tmp0 = a * wcm1 + (1.0 - a) * wcm0;
    tmp1 = b * wcm1 + (1.0 - b) * wcm0;

    double phim = phi_matrix[0][c];

    for (int i = 0; i < num_aqueous; ++i) {
      tcc_m(0) = tcc_matrix[i][c];
      if (tcc_matrix_aux != Teuchos::null) {
        for (int n = 0; n < nnodes - 1; ++n) 
          tcc_m(n + 1) = (*tcc_matrix_aux)[n][c];
      }

      msp_->second[(*msp_->first)[c]]->ComputeSoluteFlux(
          flux_liquid, tcc_next[i][c], tcc_m, 
          i, dts, tfp0, tfp1, tmp0, tmp1, phim);

      tcc_matrix[i][c] = tcc_m(0);
      if (tcc_matrix_aux != Teuchos::null) {
        for (int n = 0; n < nnodes - 1; ++n) 
          (*tcc_matrix_aux)[n][c] = tcc_m(n + 1);
      }
    }
  }
}


/* ******************************************************************* 
* Copy the advected tcc field to the state.
******************************************************************* */
void Transport_PK::CommitStep(double t_old, double t_new, const Tag& tag)
{
  auto tcc_aux = S_->GetPtrW<CV_t>(tcc_key_, Tags::DEFAULT, passwd_);
  *tcc_aux = *tcc_tmp;
}


/* ******************************************************************
* Adds source terms to conservative quantity tcc [mol]. Producers
* use the initial concentration vector tcc_prev. 
* The routine treats two cases of tcc with one and all components.
****************************************************************** */
void Transport_PK::ComputeSources_(
    double tp, double dtp, Epetra_MultiVector& tcc,
    const Epetra_MultiVector& tcc_prev, int n0, int n1)
{
  int num_vectors = tcc.NumVectors();
  int nsrcs = srcs_.size();

  for (int m = 0; m < nsrcs; m++) {
    double t0 = tp - dtp;
    srcs_[m]->Compute(t0, tp); 

    std::vector<int> index = srcs_[m]->tcc_index();
    for (auto it = srcs_[m]->begin(); it != srcs_[m]->end(); ++it) {
      int c = it->first;
      std::vector<double>& values = it->second; 

      for (int k = 0; k < index.size(); ++k) {
        int i = index[k];
        if (i < n0 || i > n1) continue;

        int imap = i;
        if (num_vectors == 1) imap = 0;

        double value = mesh_->cell_volume(c) * values[k];
        if (srcs_[m]->keyword() == "producer") {
          // correction for an extraction well
          value *= tcc_prev[imap][c];
        } else if (srcs_[m]->name() == "domain coupling") {
          value = values[k];         
        } else {
          // correction for non-SI concentration units
          if (srcs_[m]->name() == "volume" || srcs_[m]->name() == "weight")
              value /= units_.concentration_factor();
        }

        tcc[imap][c] += dtp * value;
        mass_solutes_source_[i] += value;
      }
    }
  }
}


/* *******************************************************************
* Populates operators' boundary data for given component.
* Returns true if at least one face was populated.
******************************************************************* */
bool Transport_PK::ComputeBCs_(
    std::vector<int>& bc_model, std::vector<double>& bc_value, int component)
{
  bool flag = false;

  for (int i = 0; i < bc_model.size(); i++) {
    bc_model[i] = Operators::OPERATOR_BC_NONE;
    bc_value[i] = 0.0;
  }

  AmanziMesh::Entity_ID_List cells;
  for (int f = 0; f < nfaces_wghost; f++) {
    mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
    if (cells.size() == 1) bc_model[f] = Operators::OPERATOR_BC_NEUMANN;
  }

  for (int m = 0; m < bcs_.size(); m++) {
    if (bcs_[m]->get_location() == "boundary") {
      std::vector<int>& tcc_index = bcs_[m]->tcc_index();
      int ncomp = tcc_index.size();

      for (auto it = bcs_[m]->begin(); it != bcs_[m]->end(); ++it) {
        int f = it->first;

        std::vector<double>& values = it->second; 

        for (int i = 0; i < ncomp; i++) {
          int k = tcc_index[i];
          if (k == component) {
            bc_model[f] = Operators::OPERATOR_BC_DIRICHLET;
            bc_value[f] = values[i];
            flag = true;
          }
        }
      }
    }
  }

  return flag;
}


/* *******************************************************************
* Identify flux direction based on orientation of the face normal 
* and sign of the Darcy velocity.                               
******************************************************************* */
void Transport_PK::IdentifyUpwindCells()
{
  const auto& flowrate = S_->Get<CV_t>(vol_flowrate_key_);
  flowrate.ScatterMasterToGhosted("face");

  const auto& flowrate_f = *flowrate.ViewComponent("face", true);
  const auto& map = flowrate.Map().Map("face", true);

  upwind_cells_.clear();
  downwind_cells_.clear();

  upwind_cells_.resize(nfaces_wghost);
  downwind_cells_.resize(nfaces_wghost);

  upwind_flux_.clear();
  downwind_flux_.clear();

  upwind_flux_.resize(nfaces_wghost);
  downwind_flux_.resize(nfaces_wghost);

  AmanziMesh::Entity_ID_List cells;

  // the case of fluxes that use unique face normal even if there
  // exists more than one flux on a face
  if (mesh_->space_dimension() == mesh_->manifold_dimension()) {
    for (int f = 0; f < nfaces_wghost; f++) {
      int ndofs = map->ElementSize(f);
      upwind_cells_[f].assign(ndofs, -1);
      downwind_cells_[f].assign(ndofs, -1);
      upwind_flux_[f].assign(ndofs, 0.0);
      downwind_flux_[f].assign(ndofs, 0.0);
    }
    
    for (int c = 0; c < ncells_wghost; c++) {
      const auto& faces = mesh_->cell_get_faces(c);
      const auto& dirs = mesh_->cell_get_face_dirs(c);
      
      for (int i = 0; i < faces.size(); i++) {
        int f = faces[i];
        mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);

        int g = map->FirstPointInElement(f);
        int ndofs = map->ElementSize(f);

        // We assume that two DOFs are placed only on internal faces
        if (ndofs == 2) {
          int pos = Operators::UniqueIndexFaceToCells(*mesh_, f, c);

          // define only upwind cell
          double tmp = flowrate_f[0][g + pos] * dirs[i];

          if (tmp >= 0.0) {
            upwind_cells_[f][pos] = c;
            upwind_flux_[f][pos] = flowrate_f[0][g + pos];
          } else {
            downwind_cells_[f][pos] = c;
            downwind_flux_[f][pos] = flowrate_f[0][g + pos];
          }
        }
        else {
          double tmp = flowrate_f[0][g] * dirs[i];
          if (tmp >= 0.0) {
            upwind_cells_[f][0] = c;
            upwind_flux_[f][0] = flowrate_f[0][g];
          } else if (tmp < 0.0) {
            downwind_cells_[f][0] = c;
            downwind_flux_[f][0] = flowrate_f[0][g];
          } else if (dirs[i] > 0) {
            upwind_cells_[f][0] = c;
            upwind_flux_[f][0] = flowrate_f[0][g];            
          } else {
            downwind_cells_[f][0] = c;
            downwind_flux_[f][0] = flowrate_f[0][g];
          }
        }
      }
    }

  // the case of fluxes that use the external face normal for each
  // flux (could be more than one) on a face
  } else {
    for (int c = 0; c < ncells_wghost; c++) {
      const auto& faces = mesh_->cell_get_faces(c);
      const auto& dirs = mesh_->cell_get_face_dirs(c);

      for (int i = 0; i < faces.size(); i++) {
        int f = faces[i];
        int g = map->FirstPointInElement(f);
        int ndofs = map->ElementSize(f);
        if (ndofs > 1) g += Operators::UniqueIndexFaceToCells(*mesh_, f, c);

        double u = flowrate_f[0][g] * dirs[i];  // external flux for cell c
        if (u >= 0.0) {
          upwind_cells_[f].push_back(c);
          upwind_flux_[f].push_back(u);
        } else {
          downwind_cells_[f].push_back(c);
          downwind_flux_[f].push_back(u);
        }
      }
    }
  }
}


/* *******************************************************************
* Interpolate linearly in time between two values v0 and v1. The time 
* is measuared relative to value v0; so that v1 is at time dt. The
* interpolated data are at time dt_int.            
******************************************************************* */
void Transport_PK::InterpolateCellVector(
    const Epetra_MultiVector& v0, const Epetra_MultiVector& v1, 
    double dt_int, double dt, Epetra_MultiVector& v_int) 
{
  double a = dt_int / dt;
  double b = 1.0 - a;
  v_int.Update(b, v0, a, v1, 0.0);
}

}  // namespace Transport
}  // namespace Amanzi

