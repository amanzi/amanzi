/*
  Flow PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Neil Carlson (nnc@lanl.gov), 
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <vector>

// TPLs
#include "Epetra_IntVector.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

// Amanzi
#include "BoundaryFlux.hh"
#include "EvaluatorMultiplicativeReciprocal.hh"
#include "EvaluatorPrimary.hh"
#include "CommonDefs.hh"
#include "dbc.hh"
#include "exceptions.hh"
#include "InverseFactory.hh"
#include "Mesh.hh"
#include "Mesh_Algorithms.hh"
#include "OperatorDefs.hh"
#include "PDE_DiffusionFactory.hh"
#include "PK_Utils.hh"
#include "Point.hh"
#include "UpwindFactory.hh"
#include "XMLParameterListWriter.hh"

// Amanzi::Flow
#include "DarcyVelocityEvaluator.hh"
#include "PorosityModelEvaluator.hh"
#include "Richards_PK.hh"
#include "WaterStorage.hh"
#include "WRMEvaluator.hh"
#include "WRM.hh"

namespace Amanzi {
namespace Flow {

using CV_t = CompositeVector;
using CVS_t = CompositeVectorSpace;

/* ******************************************************************
* Simplest possible constructor: extracts lists and requires fields.
****************************************************************** */
Richards_PK::Richards_PK(Teuchos::ParameterList& pk_tree,
                         const Teuchos::RCP<Teuchos::ParameterList>& glist,
                         const Teuchos::RCP<State>& S,
                         const Teuchos::RCP<TreeVector>& soln) :
  PK(pk_tree, glist, S, soln),
  Flow_PK(pk_tree, glist, S, soln),
  glist_(glist),
  soln_(soln)
{
  S_ = S;

  std::string pk_name = pk_tree.name();
  auto found = pk_name.rfind("->");
  if (found != std::string::npos) pk_name.erase(0, found + 2);

  // We need the flow list
  Teuchos::RCP<Teuchos::ParameterList> pk_list = Teuchos::sublist(glist, "PKs", true);
  fp_list_ = Teuchos::sublist(pk_list, pk_name, true);
  
  // We also need miscaleneous sublists
  preconditioner_list_ = Teuchos::sublist(glist, "preconditioners", true);
  linear_operator_list_ = Teuchos::sublist(glist, "solvers", true);
  ti_list_ = Teuchos::sublist(fp_list_, "time integrator");

  // domain name
  domain_ = fp_list_->template get<std::string>("domain name", "domain");

  vo_ = Teuchos::null;
}


/* ******************************************************************
* Old constructor for unit tests.
****************************************************************** */
Richards_PK::Richards_PK(const Teuchos::RCP<Teuchos::ParameterList>& glist,
                         const std::string& pk_list_name,
                         Teuchos::RCP<State> S,
                         const Teuchos::RCP<TreeVector>& soln) :
    Flow_PK(),
    glist_(glist),
    soln_(soln)
{
  S_ = S;

  // We need the flow list
  Teuchos::RCP<Teuchos::ParameterList> pk_list = Teuchos::sublist(glist, "PKs", true);
  fp_list_ = Teuchos::sublist(pk_list, pk_list_name, true);
 
  // We also need miscaleneous sublists
  preconditioner_list_ = Teuchos::sublist(glist, "preconditioners", true);
  linear_operator_list_ = Teuchos::sublist(glist, "solvers", true);
  ti_list_ = Teuchos::sublist(fp_list_, "time integrator");

  // domain name
  domain_ = fp_list_->template get<std::string>("domain name", "domain");

  ms_itrs_ = 0;
  ms_calls_ = 0;

  vo_ = Teuchos::null;
}


/* ******************************************************************
* Define structure of this PK. We request physical fields and their
* evaluators. Selection of a few models is available and driven by
* model factories, evaluator factories, and parameters of the list
* "physical models and assumptions".
****************************************************************** */
void Richards_PK::Setup()
{
  dt_ = 0.0;
  mesh_ = S_->GetMesh(domain_);
  dim = mesh_->space_dimension();

  // generate keys here to be available for setup of the base class
  pressure_key_ = Keys::getKey(domain_, "pressure"); 
  hydraulic_head_key_ = Keys::getKey(domain_, "hydraulic_head"); 
  darcy_velocity_key_ = Keys::getKey(domain_, "darcy_velocity"); 

  water_storage_key_ = Keys::getKey(domain_, "water_storage"); 
  prev_water_storage_key_ = Keys::getKey(domain_, "prev_water_storage"); 

  pressure_msp_key_ = Keys::getKey(domain_, "pressure_msp"); 
  porosity_msp_key_ = Keys::getKey(domain_, "porosity_msp"); 
  water_storage_msp_key_ = Keys::getKey(domain_, "water_storage_msp"); 
  prev_water_storage_msp_key_ = Keys::getKey(domain_, "prev_water_storage_msp"); 

  viscosity_liquid_key_ = Keys::getKey(domain_, "viscosity_liquid"); 
  mol_density_liquid_key_ = Keys::getKey(domain_, "molar_density_liquid"); 
  mass_density_liquid_key_ = Keys::getKey(domain_, "mass_density_liquid"); 
 
  relperm_key_ = Keys::getKey(domain_, "relative_permeability"); 
  alpha_key_ = Keys::getKey(domain_, "alpha_coef"); 

  temperature_key_ = Keys::getKey(domain_, "temperature"); 

  // set up the base class 
  key_ = pressure_key_;
  Flow_PK::Setup();

  // Our decision can be affected by the list of models
  auto physical_models = Teuchos::sublist(fp_list_, "physical models and assumptions");
  vapor_diffusion_ = physical_models->get<bool>("vapor diffusion", false);
  std::string multiscale_model = physical_models->get<std::string>("multiscale model", "single continuum");

  // Require primary field for this PK, which is pressure
  std::vector<std::string> names({ "cell" });
  std::vector<AmanziMesh::Entity_kind> locations({ AmanziMesh::CELL });
  std::vector<int> ndofs(1, 1);

  Teuchos::RCP<Teuchos::ParameterList> list1 = Teuchos::sublist(fp_list_, "operators", true);
  Teuchos::RCP<Teuchos::ParameterList> list2 = Teuchos::sublist(list1, "diffusion operator", true);
  Teuchos::RCP<Teuchos::ParameterList> list3 = Teuchos::sublist(list2, "matrix", true);
  std::string name = list3->get<std::string>("discretization primary");

  if (name != "fv: default" && name != "nlfv: default") {
    names.push_back("face");
    locations.push_back(AmanziMesh::FACE);
    ndofs.push_back(1);
  } else {
    names.push_back("boundary_face");
    locations.push_back(AmanziMesh::BOUNDARY_FACE);
    ndofs.push_back(1);
  }

  if (!S_->HasRecord(pressure_key_)) {
    S_->Require<CV_t, CVS_t>(pressure_key_, Tags::DEFAULT, passwd_)
      .SetMesh(mesh_)->SetGhosted(true)->SetComponents(names, locations, ndofs);
    AddDefaultPrimaryEvaluator_(pressure_key_);
  }

  // Require conserved quantity.
  // -- water storage
  if (!S_->HasRecord(water_storage_key_)) {
    S_->Require<CV_t, CVS_t>(water_storage_key_, Tags::DEFAULT, water_storage_key_)
      .SetMesh(mesh_)->SetGhosted(true)->SetComponent("cell", AmanziMesh::CELL, 1);

    Teuchos::ParameterList elist(water_storage_key_);
    elist.set<std::string>("water storage key", water_storage_key_)
         .set<std::string>("tag", "")
         .set<std::string>("pressure key", pressure_key_)
         .set<std::string>("saturation key", saturation_liquid_key_)
         .set<std::string>("porosity key", porosity_key_)
         .set<bool>("water vapor", vapor_diffusion_);
    if (flow_on_manifold_) 
      elist.set<std::string>("aperture key", aperture_key_);

    S_->RequireDerivative<CV_t, CVS_t>(water_storage_key_, Tags::DEFAULT,
                                       pressure_key_, Tags::DEFAULT, water_storage_key_).SetGhosted();

    auto eval = Teuchos::rcp(new WaterStorage(elist));
    S_->SetEvaluator(water_storage_key_, Tags::DEFAULT, eval);
  }

  // -- water storage from the previous time step
  if (!S_->HasRecord(prev_water_storage_key_)) {
    S_->Require<CV_t, CVS_t>(prev_water_storage_key_, Tags::DEFAULT, passwd_)
      .SetMesh(mesh_)->SetGhosted(true)->SetComponent("cell", AmanziMesh::CELL, 1);
    S_->GetRecordW(prev_water_storage_key_, passwd_).set_io_vis(false);
  }

  // -- multiscale extension: secondary (immobile water storage)
  if (multiscale_model == "dual continuum discontinuous matrix") {
    if (!S_->HasRecord(pressure_msp_key_)) {
      S_->Require<CV_t, CVS_t>(pressure_msp_key_, Tags::DEFAULT, passwd_)
        .SetMesh(mesh_)->SetGhosted(false)->SetComponent("cell", AmanziMesh::CELL, 1);

      Teuchos::ParameterList elist(pressure_msp_key_);
      elist.set<std::string>("evaluator name", pressure_msp_key_);
      pressure_msp_eval_ = Teuchos::rcp(new EvaluatorPrimary<CV_t, CVS_t>(elist));
      S_->SetEvaluator(pressure_msp_key_, Tags::DEFAULT, pressure_msp_eval_);
    }

    Teuchos::RCP<Teuchos::ParameterList> msp_list = Teuchos::sublist(fp_list_, "multiscale models", true);
    msp_ = CreateMultiscaleFlowPorosityPartition(mesh_, msp_list);

    if (!S_->HasRecord(water_storage_msp_key_)) {
      S_->Require<CV_t, CVS_t>(water_storage_msp_key_, Tags::DEFAULT, passwd_)
        .SetMesh(mesh_)->SetGhosted(true)->SetComponent("cell", AmanziMesh::CELL, 1);
    }
    if (!S_->HasRecord(prev_water_storage_msp_key_)) {
      S_->Require<CV_t, CVS_t>(prev_water_storage_msp_key_, Tags::DEFAULT, passwd_)
        .SetMesh(mesh_)->SetGhosted(true)->SetComponent("cell", AmanziMesh::CELL, 1);
      S_->GetRecordW(prev_water_storage_msp_key_, passwd_).set_io_vis(false);
    }

    S_->Require<CV_t, CVS_t>(porosity_msp_key_, Tags::DEFAULT, porosity_msp_key_)
      .SetMesh(mesh_)->SetGhosted(false)->SetComponent("cell", AmanziMesh::CELL, 1);

    Teuchos::ParameterList elist(porosity_msp_key_);
    elist.set<std::string>("porosity key", porosity_msp_key_)
         .set<std::string>("pressure key", pressure_msp_key_)
         .set<std::string>("tag", "");

    Teuchos::RCP<PorosityModelPartition> pom = CreatePorosityModelPartition(mesh_, msp_list);
    auto eval = Teuchos::rcp(new PorosityModelEvaluator(elist, pom));
    S_->SetEvaluator(porosity_msp_key_, Tags::DEFAULT, eval);

    // Secondary matrix nodes are collected here.
    int nnodes, nnodes_tmp = NumberMatrixNodes(msp_);
    mesh_->get_comm()->MaxAll(&nnodes_tmp, &nnodes, 1);
    if (nnodes > 1) {
      S_->Require<CV_t, CVS_t>("water_storage_ms_aux", Tags::DEFAULT, passwd_)
        .SetMesh(mesh_)->SetGhosted(false)->SetComponent("cell", AmanziMesh::CELL, nnodes - 1);
      S_->GetRecordW("water_storage_ms_aux", passwd_).set_io_vis(false);

      S_->Require<CV_t, CVS_t>("pressure_msp_aux", Tags::DEFAULT, passwd_)
        .SetMesh(mesh_)->SetGhosted(false)->SetComponent("cell", AmanziMesh::CELL, nnodes - 1);
      S_->GetRecordW("pressure_msp_aux", passwd_).set_io_vis(false);

      S_->Require<CV_t, CVS_t>("porosity_msp_aux", Tags::DEFAULT, "porosity_msp_aux")
        .SetMesh(mesh_)->SetGhosted(false)->SetComponent("cell", AmanziMesh::CELL, nnodes - 1);
      S_->SetEvaluator("porosity_msp_aux", Tags::DEFAULT, eval);
      S_->GetRecordW("porosity_msp_aux", passwd_).set_io_vis(false);
    }
  }

  // Require additional fields and evaluators for this PK.
  // -- porosity
  if (!S_->HasRecord(porosity_key_)) {
    S_->Require<CV_t, CVS_t>(porosity_key_, Tags::DEFAULT, porosity_key_)
      .SetMesh(mesh_)->SetGhosted(true)->SetComponent("cell", AmanziMesh::CELL, 1);

    std::string pom_name = physical_models->get<std::string>("porosity model", "constant porosity");

    if (pom_name == "compressible: pressure function") {
      Teuchos::RCP<Teuchos::ParameterList>
          pom_list = Teuchos::sublist(fp_list_, "porosity models", true);
      Teuchos::RCP<PorosityModelPartition> pom = CreatePorosityModelPartition(mesh_, pom_list);

      Teuchos::ParameterList elist(porosity_key_);
      elist.set<std::string>("porosity key", porosity_key_)
           .set<std::string>("pressure key", pressure_key_)
           .set<std::string>("tag", "");
      // elist.sublist("verbose object").set<std::string>("verbosity level", "extreme");

      S_->RequireDerivative<CV_t, CVS_t>(porosity_key_, Tags::DEFAULT,
                                         pressure_key_, Tags::DEFAULT, porosity_key_).SetGhosted();

      auto eval = Teuchos::rcp(new PorosityModelEvaluator(elist, pom));
      S_->SetEvaluator(porosity_key_, Tags::DEFAULT, eval);
    } else {
      S_->RequireEvaluator(porosity_key_, Tags::DEFAULT);
    }
  }

  // -- viscosity: if not requested by any PK, we request its constant value.
  if (!S_->HasRecord(viscosity_liquid_key_)) {
    if (!S_->HasRecord("const_fluid_viscosity")) {
      S_->Require<double>("const_fluid_viscosity", Tags::DEFAULT, "state");
    }
    S_->Require<CV_t, CVS_t>(viscosity_liquid_key_, Tags::DEFAULT, viscosity_liquid_key_)
      .SetMesh(mesh_)->SetGhosted(true)->AddComponent("cell", AmanziMesh::CELL, 1)
      ->AddComponent("boundary_face", AmanziMesh::BOUNDARY_FACE, 1);

    double mu = glist_->sublist("state").sublist("initial conditions")
                       .sublist("const_fluid_viscosity").get<double>("value");

    std::vector<std::string> components({ "cell", "boundary_face" });
    Teuchos::ParameterList& eval = S_->FEList().sublist(viscosity_liquid_key_);
    eval.sublist("function").sublist("DOMAIN")
        .set<std::string>("region", "All")
        .set<Teuchos::Array<std::string> >("components", components)
        .sublist("function").sublist("function-constant")
        .set<double>("value", mu);
    eval.set<std::string>("evaluator type", "independent variable");

    S_->RequireEvaluator(viscosity_liquid_key_, Tags::DEFAULT);
  }

  // -- model for liquid density is constant density unless specified otherwise
  //    in high-level PKs.
  if (!S_->HasRecord(mol_density_liquid_key_)) {
    S_->Require<CV_t, CVS_t>(mol_density_liquid_key_, Tags::DEFAULT, mol_density_liquid_key_)
      .SetMesh(mesh_)->SetGhosted(true)
      ->AddComponent("cell", AmanziMesh::CELL, 1)
      ->AddComponent("boundary_face", AmanziMesh::BOUNDARY_FACE, 1);

    S_->RequireEvaluator(mol_density_liquid_key_, Tags::DEFAULT);
  }

  if (!S_->HasRecord(mass_density_liquid_key_)) {
    S_->Require<CV_t, CVS_t>(mass_density_liquid_key_, Tags::DEFAULT, mass_density_liquid_key_)
      .SetMesh(mesh_)->SetGhosted(true)
      ->AddComponent("cell", AmanziMesh::CELL, 1)
      ->AddComponent("boundary_face", AmanziMesh::BOUNDARY_FACE, 1);
    S_->RequireEvaluator(mass_density_liquid_key_, Tags::DEFAULT);
  }
  
  // -- saturation
  auto wrm_list = Teuchos::sublist(fp_list_, "water retention models", true);
  wrm_ = CreateWRMPartition(mesh_, wrm_list);

  if (!S_->HasRecord(saturation_liquid_key_)) {
    S_->Require<CV_t, CVS_t>(saturation_liquid_key_, Tags::DEFAULT, saturation_liquid_key_)
      .SetMesh(mesh_)->SetGhosted(true)->SetComponent("cell", AmanziMesh::CELL, 1);

    Teuchos::ParameterList elist(saturation_liquid_key_);
    elist.set<std::string>("saturation key", saturation_liquid_key_)
         .set<std::string>("pressure key", pressure_key_)
         .set<std::string>("tag", "");

    // elist.sublist("verbose object").set<std::string>("verbosity level", "extreme");
    auto eval = Teuchos::rcp(new WRMEvaluator(elist, wrm_));
    S_->SetEvaluator(saturation_liquid_key_, Tags::DEFAULT, eval);
  }

  if (!S_->HasRecord(prev_saturation_liquid_key_)) {
    S_->Require<CV_t, CVS_t>(prev_saturation_liquid_key_, Tags::DEFAULT, passwd_)
      .SetMesh(mesh_)->SetGhosted(true)->SetComponent("cell", AmanziMesh::CELL, 1);
    S_->GetRecordW(prev_saturation_liquid_key_, passwd_).set_io_vis(false);
  }

  // -- relative permeability
  if (!S_->HasRecord(relperm_key_)) {
    S_->Require<CV_t, CVS_t>(relperm_key_, Tags::DEFAULT, relperm_key_)
      .SetMesh(mesh_)->SetGhosted(true)->AddComponent("cell", AmanziMesh::CELL, 1)
      ->AddComponent("boundary_face", AmanziMesh::BOUNDARY_FACE, 1);

    Teuchos::ParameterList elist(relperm_key_);
    elist.set<std::string>("relative permeability key", relperm_key_)
         .set<std::string>("tag", "");

    S_->RequireDerivative<CV_t, CVS_t>(relperm_key_, Tags::DEFAULT,
                                       pressure_key_, Tags::DEFAULT, relperm_key_).SetGhosted();

    auto eval = Teuchos::rcp(new Flow::RelPermEvaluator(elist, S_.ptr(), wrm_));
    S_->SetEvaluator(relperm_key_, Tags::DEFAULT, eval);
  }

  // -- inverse of kinematic viscosity
  if (!S_->HasRecord(alpha_key_)) {
    Key kkey;
    if (flow_on_manifold_) {
      S_->Require<CV_t, CVS_t>(alpha_key_, Tags::DEFAULT, alpha_key_)
        .SetMesh(mesh_)->SetGhosted(true)->AddComponent("cell", AmanziMesh::CELL, 1);
      kkey = permeability_key_;
    } else {
      S_->Require<CV_t, CVS_t>(alpha_key_, Tags::DEFAULT, alpha_key_)
        .SetMesh(mesh_)->SetGhosted(true)
        ->AddComponent("cell", AmanziMesh::CELL, 1)
        ->AddComponent("boundary_face", AmanziMesh::BOUNDARY_FACE, 1);
      kkey = relperm_key_;
    }

    std::vector<std::string> listm({ Keys::getVarName(kkey), Keys::getVarName(mol_density_liquid_key_) });
    std::vector<std::string> listr({ Keys::getVarName(viscosity_liquid_key_) });
    if (flow_on_manifold_) listm.push_back(Keys::getVarName(aperture_key_));

    Teuchos::ParameterList elist(alpha_key_);
    elist.set<std::string>("my key", alpha_key_)
         .set<Teuchos::Array<std::string> >("multiplicative dependencies", listm)
         .set<Teuchos::Array<std::string> >("reciprocal dependencies", listr)
         .set<std::string>("tag", "");

    S_->RequireDerivative<CV_t, CVS_t>(alpha_key_, Tags::DEFAULT,
                                       pressure_key_, Tags::DEFAULT, alpha_key_).SetGhosted();

    auto eval = Teuchos::rcp(new EvaluatorMultiplicativeReciprocal(elist));
    S_->SetEvaluator(alpha_key_, Tags::DEFAULT, eval);
  }

  // Local fields and evaluators.
  // -- hydraulic head
  if (!S_->HasRecord(hydraulic_head_key_)) {
    S_->Require<CV_t, CVS_t>(hydraulic_head_key_, Tags::DEFAULT, passwd_)
      .SetMesh(mesh_)->SetGhosted(true)->SetComponent("cell", AmanziMesh::CELL, 1);
  }

  // -- Darcy velocity vector
  if (!S_->HasRecord(darcy_velocity_key_)) {
    S_->Require<CV_t, CVS_t>(darcy_velocity_key_, Tags::DEFAULT, darcy_velocity_key_)
      .SetMesh(mesh_)->SetGhosted(true)->SetComponent("cell", AmanziMesh::CELL, dim);

    Teuchos::ParameterList elist(darcy_velocity_key_);
    elist.set<std::string>("domain name", domain_);
    elist.set<std::string>("darcy velocity key", darcy_velocity_key_)
         .set<std::string>("volumetric flow rate key", vol_flowrate_key_)
         .set<std::string>("tag", "");

    auto eval = Teuchos::rcp(new DarcyVelocityEvaluator(elist));
    S_->SetEvaluator(darcy_velocity_key_, Tags::DEFAULT, eval);
  }

  // -- temperature for EOSs
  if (!S_->HasRecord(temperature_key_)) {
    auto cvs = S_->GetRecordSetW(pressure_key_).GetFactory<CV_t, CVS_t>();
    *S_->Require<CV_t, CVS_t>(temperature_key_, Tags::DEFAULT, passwd_)
      .SetMesh(mesh_)->SetGhosted(true) = cvs;
  }

  if (!S_->HasEvaluator(temperature_key_, Tags::DEFAULT)) {
    AddDefaultIndependentEvaluator_(temperature_key_, Tags::DEFAULT, 298.15);
  }

  // Require additional components for the existing fields
  Teuchos::ParameterList abs_perm = fp_list_->sublist("absolute permeability");
  coordinate_system_ = abs_perm.get<std::string>("coordinate system", "cartesian");
  int noff = abs_perm.get<int>("off-diagonal components", 0);
 
  if (noff > 0) {
    CompositeVectorSpace& cvs = S_->Require<CV_t, CVS_t>(permeability_key_, Tags::DEFAULT, permeability_key_);
    cvs.SetOwned(false);
    cvs.AddComponent("offd", AmanziMesh::CELL, noff)->SetOwned(true);
  }

  // Since high-level PK may own some fields, we have to populate 
  // frequently used evaluators outside of field registration
  pressure_eval_ = Teuchos::rcp_dynamic_cast<EvaluatorPrimary<CV_t, CVS_t> >(
      S_->GetEvaluatorPtr(pressure_key_, Tags::DEFAULT));
  vol_flowrate_eval_ = Teuchos::rcp_dynamic_cast<EvaluatorPrimary<CV_t, CVS_t> >(
      S_->GetEvaluatorPtr(vol_flowrate_key_, Tags::DEFAULT));
}


/* ******************************************************************
* This is a long but simple routine. It goes through flow parameter
* list and initializes various objects including those created during 
* the setup step.
****************************************************************** */
void Richards_PK::Initialize()
{
  // Initialize miscalleneous defaults.
  // -- times
  double t_ini = S_->get_time(); 
  dt_desirable_ = dt_;
  dt_next_ = dt_;

  // -- others
  error_control_ = FLOW_TI_ERROR_CONTROL_PRESSURE;

  mass_bc = 0.0;
  seepage_mass_ = 0.0;
  mass_initial = 0.0;
  initialize_with_darcy_ = true;
  num_itrs_ = 0;

  // Create verbosity object to print out initialiation statistics.
  Teuchos::ParameterList vlist;
  vlist.sublist("verbose object") = fp_list_->sublist("verbose object");

  std::string ioname = "RichardsPK";
  if (domain_ != "domain") ioname += "-" + domain_;
  vo_ = Teuchos::rcp(new VerboseObject(ioname, vlist)); 

  if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os()<< "\nPK initialization started...\n";
  }

  // Initilize various base class data.
  Flow_PK::Initialize();

  // Initialize local fields and evaluators.
  InitializeFields_();
  UpdateLocalFields_(S_.ptr());

  // Create BCs and source terms.
  InitializeBCsSources_(*fp_list_);

  // relative permeability
  // -- create basic fields, factories and control variables
  Teuchos::RCP<Teuchos::ParameterList> upw_list = Teuchos::sublist(fp_list_, "relative permeability", true);

  Operators::UpwindFactory upwind_factory;
  upwind_ = upwind_factory.Create(mesh_, *upw_list);

  std::string upw_upd = upw_list->get<std::string>("upwind frequency", "every timestep");
  if (upw_upd == "every nonlinear iteration") upwind_frequency_ = FLOW_UPWIND_UPDATE_ITERATION;
  else upwind_frequency_ = FLOW_UPWIND_UPDATE_TIMESTEP;  

  // face component of upwind field matches that of the flow field 
  auto cvs = S_->Get<CV_t>(vol_flowrate_key_).Map();
  cvs.SetOwned(false);
  cvs.AddComponent("cell", AmanziMesh::CELL, 1);
  alpha_upwind_ = Teuchos::rcp(new CompositeVector(cvs));
  alpha_upwind_dP_ = Teuchos::rcp(new CompositeVector(cvs));

  // Process models and assumptions.
  flux_units_ = molar_rho_ / rho_;

  // -- coupling with other physical PKs
  Teuchos::RCP<Teuchos::ParameterList> physical_models = 
      Teuchos::sublist(fp_list_, "physical models and assumptions");
  multiscale_porosity_ = (physical_models->get<std::string>(
      "multiscale model", "single continuum") != "single continuum");

  // Select a proper matrix class. 
  const Teuchos::ParameterList& tmp_list = fp_list_->sublist("operators")
                                                    .sublist("diffusion operator");
  Teuchos::ParameterList oplist_matrix = tmp_list.sublist("matrix");
  Teuchos::ParameterList oplist_pc = tmp_list.sublist("preconditioner");

  std::string name = fp_list_->sublist("relative permeability").get<std::string>("upwind method");
  std::string nonlinear_coef("standard: cell");
  if (flow_on_manifold_) {
    nonlinear_coef = "standard: cell";
  } else if (name == "upwind: darcy velocity") {
    nonlinear_coef = "upwind: face";
  } else if (name == "upwind: gravity") {
    nonlinear_coef = "upwind: face";
  } else if (name == "upwind: amanzi" || name == "upwind: amanzi new") {
    nonlinear_coef = "divk: cell-face";
    // nonlinear_coef = "divk: face";
  } else if (name == "other: arithmetic average") {
    nonlinear_coef = "upwind: face";
  }
  if (!oplist_matrix.isParameter("nonlinear coefficient")) {
    oplist_matrix.set<std::string>("nonlinear coefficient", nonlinear_coef);
    oplist_pc.set<std::string>("nonlinear coefficient", nonlinear_coef);
  }

  // auto rho_cv = S_->GetPtr<CV_t>(mass_density_liquid, Tags::DEFAULT);

  Operators::PDE_DiffusionFactory opfactory(oplist_matrix, mesh_);
  opfactory.SetConstantGravitationalTerm(gravity_, rho_);

  if (!flow_on_manifold_) {
    SetAbsolutePermeabilityTensor();
    Teuchos::RCP<std::vector<WhetStone::Tensor> > Kptr = Teuchos::rcpFromRef(K);
    opfactory.SetVariableTensorCoefficient(Kptr);
    opfactory.SetVariableScalarCoefficient(alpha_upwind_, alpha_upwind_dP_);
  } else {
    auto kptr = S_->GetPtrW<CV_t>(alpha_key_, Tags::DEFAULT, alpha_key_);
    opfactory.SetVariableScalarCoefficient(kptr);
  }

  op_matrix_diff_ = opfactory.Create();
  op_matrix_ = op_matrix_diff_->global_operator();

  opfactory.SetPList(oplist_pc);
  op_preconditioner_diff_ = opfactory.Create();
  op_preconditioner_ = op_preconditioner_diff_->global_operator();

  op_acc_ = Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::CELL, op_preconditioner_));

  if (vapor_diffusion_) {
    Teuchos::ParameterList oplist_vapor = tmp_list.sublist("vapor matrix");
    op_vapor_diff_ = opfactory.Create(oplist_vapor, mesh_, op_bc_);
    op_vapor_ = op_vapor_diff_->global_operator();
    op_preconditioner_->OpPushBack(op_vapor_diff_->local_op());
  }

  // Create pointers to the primary flow field pressure.
  solution = S_->GetPtrW<CV_t>(pressure_key_, Tags::DEFAULT, passwd_);
  soln_->SetData(solution); 
  
  // Create auxiliary vectors for time history and error estimates.
  const Epetra_BlockMap& cmap_owned = mesh_->cell_map(false);
  pdot_cells_prev = Teuchos::rcp(new Epetra_Vector(cmap_owned));
  pdot_cells = Teuchos::rcp(new Epetra_Vector(cmap_owned));

  // Initialize flux copy for the upwind operator.
  vol_flowrate_copy = Teuchos::rcp(new CompositeVector(S_->Get<CV_t>(vol_flowrate_key_)));

  // Conditional initialization of lambdas from pressures.
  auto& pressure = S_->GetW<CV_t>(pressure_key_, Tags::DEFAULT, passwd_);

  if (ti_list_->isSublist("pressure-lambda constraints") && pressure.HasComponent("face")) {
    DeriveFaceValuesFromCellValues(*pressure.ViewComponent("cell"),
                                   *pressure.ViewComponent("face"));
    S_->GetRecordW(pressure_key_, passwd_).set_initialized(true);
  }

  // error control options
  AMANZI_ASSERT(ti_list_->isParameter("error control options"));

  error_control_ = 0;
  std::vector<std::string> options;
  options = ti_list_->get<Teuchos::Array<std::string> >("error control options").toVector();

  for (int i=0; i < options.size(); i++) {
    if (options[i] == "pressure") {
      error_control_ += FLOW_TI_ERROR_CONTROL_PRESSURE;
    } else if (options[i] == "saturation") {
      error_control_ += FLOW_TI_ERROR_CONTROL_SATURATION;
    } else if (options[i] == "residual") {
      error_control_ += FLOW_TI_ERROR_CONTROL_RESIDUAL;
    }
  }

  if (error_control_ == 0) {
    error_control_ = FLOW_TI_ERROR_CONTROL_PRESSURE +  // usually 1 [Pa]
                     FLOW_TI_ERROR_CONTROL_SATURATION;  // usually 1e-4;
  }

  // initialize time integrator
  std::string ti_method_name = ti_list_->get<std::string>("time integration method", "none");
  if (ti_method_name == "BDF1") {
    Teuchos::ParameterList& bdf1_list = ti_list_->sublist("BDF1");

    if (! bdf1_list.isSublist("verbose object"))
        bdf1_list.sublist("verbose object") = fp_list_->sublist("verbose object");

    bdf1_dae_ = Teuchos::rcp(new BDF1_TI<TreeVector, TreeVectorSpace>(*this, bdf1_list, soln_));
  } else {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "WARNING: BDF1 time integration list is missing..." << std::endl;
  }

  // Initialize boundary conditions and source terms.
  UpdateSourceBoundaryData(t_ini, t_ini, pressure);

  // Initialize matrix and preconditioner operators.
  // -- molar density requires to rescale gravity later.
  //    make an evaluator for alpha_upwind? FIXME
  if (!flow_on_manifold_) {
    auto& alpha = S_->GetW<CompositeVector>(alpha_key_, Tags::DEFAULT, alpha_key_);
    *alpha_upwind_->ViewComponent("cell") = *alpha.ViewComponent("cell");
  }

  op_matrix_->Init();
  op_matrix_diff_->SetBCs(op_bc_, op_bc_);
  op_matrix_diff_->UpdateMatrices(Teuchos::null, solution.ptr());
  op_matrix_diff_->ApplyBCs(true, true, true);

  op_preconditioner_->Init();
  op_preconditioner_diff_->SetBCs(op_bc_, op_bc_);
  op_preconditioner_diff_->UpdateMatrices(vol_flowrate_copy.ptr(), solution.ptr());
  op_preconditioner_diff_->UpdateMatricesNewtonCorrection(vol_flowrate_copy.ptr(), solution.ptr(), molar_rho_);
  op_preconditioner_diff_->ApplyBCs(true, true, true);

  if (vapor_diffusion_) {
    // op_vapor_diff_->SetBCs(op_bc_);
    op_vapor_diff_->SetScalarCoefficient(Teuchos::null, Teuchos::null);
  }

  // -- generic linear solver for most cases

  // -- preconditioner or encapsulated preconditioner
  std::string pc_name = ti_list_->get<std::string>("preconditioner");
  op_preconditioner_->set_inverse_parameters(pc_name, *preconditioner_list_);
  
  // Optional step: calculate hydrostatic solution consistent with BCs
  // and clip it as requested. We have to do it only once at the beginning
  // of time period.
  if (ti_list_->isSublist("initialization") && initialize_with_darcy_ 
      && S_->get_position() == Amanzi::TIME_PERIOD_START) {
    initialize_with_darcy_ = false;
    Teuchos::ParameterList& ini_list = ti_list_->sublist("initialization");
 
    std::string ini_method_name = ini_list.get<std::string>("method", "none");
    if (ini_method_name == "saturated solver") {
      name = ini_list.get<std::string>("linear solver");
      SolveFullySaturatedProblem(t_ini, *solution, name);

      bool clip(false);
      double clip_saturation = ini_list.get<double>("clipping saturation value", -1.0);
      if (clip_saturation > 0.0) {
        double pmin = atm_pressure_;
        Epetra_MultiVector& p = *solution->ViewComponent("cell");
        ClipHydrostaticPressure(pmin, clip_saturation, p);
        clip = true;
      }

      double clip_pressure = ini_list.get<double>("clipping pressure value", -1e+10);
      if (clip_pressure > -5 * atm_pressure_) {
        Epetra_MultiVector& p = *solution->ViewComponent("cell");
        ClipHydrostaticPressure(clip_pressure, p);
        clip = true;
      }

      if (clip && vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
        Teuchos::OSTab tab = vo_->getOSTab();
        *vo_->os() << "\nClipped pressure field.\n" << std::endl;
      }

      if (clip && solution->HasComponent("face")) {
        Epetra_MultiVector& p = *solution->ViewComponent("cell");
        Epetra_MultiVector& lambda = *solution->ViewComponent("face", true);
        DeriveFaceValuesFromCellValues(p, lambda);
      }
    }
    else if (ini_method_name == "picard") {
      AdvanceToSteadyState_Picard(ti_list_->sublist("initialization"));
    }
    pressure_eval_->SetChanged();

    // initialization is usually done at time 0, so we need to update other
    // fields such as prev_saturation_liquid
    S_->GetEvaluator(saturation_liquid_key_).Update(*S_, "flow");
    auto& s_l = S_->GetW<CV_t>(saturation_liquid_key_, Tags::DEFAULT, saturation_liquid_key_);
    auto& s_l_prev = S_->GetW<CV_t>(prev_saturation_liquid_key_, Tags::DEFAULT, passwd_);
    s_l_prev = s_l;

    S_->GetEvaluator(water_storage_key_).Update(*S_, "flow");
    auto& wc = S_->GetW<CV_t>(water_storage_key_, Tags::DEFAULT, water_storage_key_);
    auto& wc_prev = S_->GetW<CV_t>(prev_water_storage_key_, Tags::DEFAULT, passwd_);
    wc_prev = wc;

    // We start with pressure equilibrium
    if (multiscale_porosity_) {
      *S_->GetW<CV_t>(pressure_msp_key_, passwd_).ViewComponent("cell") =
          *S_->Get<CV_t>(pressure_key_).ViewComponent("cell");
      pressure_msp_eval_->SetChanged();
    }
  }

  // Trigger update of secondary fields depending on the primary pressure.
  pressure_eval_->SetChanged();

  // Derive volumetric flow rate (state may not have it at time 0)
  double tmp;
  vol_flowrate_copy->Norm2(&tmp);
  if (tmp == 0.0) {
    op_matrix_diff_->UpdateFlux(solution.ptr(), vol_flowrate_copy.ptr());
    vol_flowrate_copy->ScaleMasterAndGhosted(1.0 / molar_rho_);
  }

  // Subspace entering: re-initialize lambdas.
  if (ti_list_->isSublist("pressure-lambda constraints") && solution->HasComponent("face")) {
    solver_name_constraint_ = ti_list_->sublist("pressure-lambda constraints").get<std::string>("linear solver");

    if (S_->get_position() == Amanzi::TIME_PERIOD_START) {
      EnforceConstraints(t_ini, solution);
      pressure_eval_->SetChanged();

      // update mass flux
      op_matrix_->Init();
      op_matrix_diff_->UpdateMatrices(Teuchos::null, solution.ptr());
      op_matrix_diff_->UpdateFlux(solution.ptr(), vol_flowrate_copy.ptr());
      vol_flowrate_copy->ScaleMasterAndGhosted(1.0 / molar_rho_);
    }
  }

  // Development: miscalleneous
  algebraic_water_storage_balance_ = fp_list_->get<bool>("algebraic water storage balance", false);
  if (algebraic_water_storage_balance_) {
    CompositeVectorSpace cvs1; 
    cvs1.SetMesh(mesh_)->SetGhosted(false)
        ->AddComponent("cell", AmanziMesh::CELL, 1)
        ->AddComponent("dpre", AmanziMesh::CELL, 1);
    cnls_limiter_ = Teuchos::rcp(new CompositeVector(cvs1));
  }

  // NOTE: this is alternatively called "linear solver" or "preconditioner
  // enhancement".  One got stuffed into op_pc_solver_, the other gets
  // constructed in, e.g. AdvanceToSteadyState_Picard.  Must these be separate?
  // They can be... and so I kept them separate for now.  But this means that
  // all of flow PK cannot have linear solver in the Operator. --etc
  // solver_name_ = ti_list_->get<std::string>("linear solver");
  // op_preconditioner_->set_inverse_parameters(pc_name, *preconditioner_list_,
  //         solver_name_, *linear_operator_list_, true);
  std::string tmp_solver = ti_list_->get<std::string>("preconditioner enhancement", "none");
  if (tmp_solver != "none") {
    AMANZI_ASSERT(linear_operator_list_->isSublist(tmp_solver));
    Teuchos::ParameterList tmp_plist = linear_operator_list_->sublist(tmp_solver);
    op_pc_solver_ = AmanziSolvers::createIterativeMethod(tmp_plist, op_preconditioner_);
  } else {
    op_pc_solver_ = op_preconditioner_;
  }
  op_pc_solver_->InitializeInverse();

  // Verbose output of initialization statistics.
  InitializeStatistics_();
}


/* ****************************************************************
* This completes initialization of common fields that were not 
* initialized by the state.
**************************************************************** */
void Richards_PK::InitializeFields_()
{
  Teuchos::OSTab tab = vo_->getOSTab();

  // set popular default values for missed fields.
  if (S_->GetRecord(saturation_liquid_key_).owner() == passwd_) {
    if (S_->HasRecord(saturation_liquid_key_)) {
      if (!S_->GetRecord(saturation_liquid_key_).initialized()) {
        S_->GetW<CV_t>(saturation_liquid_key_, Tags::DEFAULT, passwd_).PutScalar(1.0);
        S_->GetRecordW(saturation_liquid_key_, Tags::DEFAULT, passwd_).set_initialized();

        if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM)
            *vo_->os() << "initialized saturation_liquid to default value 1.0" << std::endl;  
      }
    }
  }

  InitializeFieldFromField_(prev_saturation_liquid_key_, saturation_liquid_key_, true);
  InitializeFieldFromField_(prev_water_storage_key_, water_storage_key_, true);
  InitializeFieldFromField_(prev_aperture_key_, aperture_key_, true);

  // set matrix fields assuming pressure equilibrium
  // -- pressure
  if (S_->HasRecord(pressure_msp_key_)) {
    // if (!S_->GetField(pressure_msp_key_, passwd_)->initialized()) {
      const auto& p1 = *S_->Get<CV_t>(pressure_key_).ViewComponent("cell");
      auto& p0 = *S_->GetW<CV_t>(pressure_msp_key_, passwd_).ViewComponent("cell");
      p0 = p1;

      S_->GetRecordW(pressure_msp_key_, passwd_).set_initialized();
      pressure_msp_eval_->SetChanged();

      if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM)
          *vo_->os() << "initialized pressure_msp to pressure" << std::endl;  
    // }
  }

  // -- water storage
  if (S_->HasRecord(water_storage_msp_key_)) {
    if (!S_->GetRecord(water_storage_msp_key_, Tags::DEFAULT).initialized()) {
      CalculateWaterStorageMultiscale_();
      S_->GetRecordW(water_storage_msp_key_, Tags::DEFAULT, passwd_).set_initialized();

      if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM)
          *vo_->os() << "initialized water_storage_msp to WaterStorage(pressure_msp)" << std::endl;  
    }
  }

  InitializeFieldFromField_(prev_water_storage_msp_key_, water_storage_msp_key_, false);
}


/* ******************************************************************
* Print the header for new time period.
****************************************************************** */
void Richards_PK::InitializeStatistics_()
{
  if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
    std::string ti_method_name = ti_list_->get<std::string>("time integration method");
    std::string pc_name = ti_list_->get<std::string>("preconditioner");

    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os()<< "\nEC:" << error_control_ 
              << " Upwind:" << op_matrix_diff_->little_k()
              << " PC:\"" << pc_name.c_str() << "\"" 
              << " TI:\"" << ti_method_name.c_str() << "\"" << std::endl
              << "matrix: " << op_matrix_->PrintDiagnostics() << std::endl
              << "precon: " << op_preconditioner_->PrintDiagnostics() << std::endl;

    int missed_tmp = missed_bc_faces_;
    int dirichlet_tmp = dirichlet_bc_faces_;
#ifdef HAVE_MPI
    mesh_->get_comm()->SumAll(&missed_tmp, &missed_bc_faces_, 1);
    mesh_->get_comm()->SumAll(&dirichlet_tmp, &dirichlet_bc_faces_, 1);
#endif

    *vo_->os() << "pressure BC assigned to " << dirichlet_bc_faces_ << " faces" << std::endl;
    *vo_->os() << "default (no-flow) BC assigned to " << missed_bc_faces_ << " faces" << std::endl << std::endl;

    VV_PrintHeadExtrema(*solution);
    VV_PrintSourceExtrema();

    if (vo_->getVerbLevel() >= Teuchos::VERB_EXTREME) {
      *vo_->os() << "\nPrinting WRM files..." << std::endl;
      PlotWRMcurves_();
    }

    *vo_->os() << vo_->color("green") << "Initialization of PK is complete, T=" 
               << units_.OutputTime(S_->get_time()) << vo_->reset() << std::endl << std::endl;
  }

  if (dirichlet_bc_faces_ == 0 &&
      domain_ == "domain" && vo_->getVerbLevel() >= Teuchos::VERB_LOW) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "WARNING: no essential boundary conditions, solver may fail" << std::endl;
  }
}


/* ******************************************************************* 
* Performs one time step from time t_old to time t_new either for
* steady-state or transient simulation. If reinit=true, enforce 
* p-lambda constraints.
******************************************************************* */
bool Richards_PK::AdvanceStep(double t_old, double t_new, bool reinit)
{
  AMANZI_ASSERT(bdf1_dae_ != Teuchos::null);

  dt_ = t_new - t_old;

  // initialize statistics
  ms_itrs_ = 0;
  ms_calls_ = 0;

  // save a copy of primary and conservative fields
  std::vector<std::string> fields;
  if (flow_on_manifold_) {
    fields.push_back(prev_aperture_key_);
  }

  StateArchive archive(S_, vo_);
  archive.Add(fields, {}, {}, Tags::DEFAULT, name());
  archive.Swap("");

  std::map<std::string, CompositeVector> copies;
  copies.emplace(pressure_key_, S_->Get<CV_t>(pressure_key_, Tags::DEFAULT));

  // -- saturations, swap prev <- current
  S_->GetEvaluator(saturation_liquid_key_).Update(*S_, "flow");
  const auto& sat = S_->Get<CV_t>(saturation_liquid_key_);
  auto& sat_prev = S_->GetW<CV_t>(prev_saturation_liquid_key_, Tags::DEFAULT, passwd_);

  copies.emplace(prev_saturation_liquid_key_, sat_prev);
  sat_prev = sat;

  // -- water_storage, swap prev <- current
  S_->GetEvaluator(water_storage_key_).Update(*S_, "flow");
  auto& wc = S_->GetW<CV_t>(water_storage_key_, Tags::DEFAULT, water_storage_key_);
  auto& wc_prev = S_->GetW<CV_t>(prev_water_storage_key_, Tags::DEFAULT, passwd_);

  copies.emplace(prev_water_storage_key_, wc_prev);
  wc_prev = wc;

  // -- field for multiscale models, save and swap
  Teuchos::RCP<CompositeVector> pressure_msp_copy, wc_matrix_prev_copy;
  if (multiscale_porosity_) {
    copies.emplace(pressure_msp_key_, S_->Get<CV_t>(pressure_msp_key_));

    auto& wc_matrix = S_->GetW<CV_t>(water_storage_msp_key_, Tags::DEFAULT, passwd_);
    auto& wc_matrix_prev = S_->GetW<CV_t>(prev_water_storage_msp_key_, Tags::DEFAULT, passwd_);

    copies.emplace(prev_water_storage_msp_key_, wc_matrix_prev);
    wc_matrix_prev = wc_matrix;
  }

  // enter subspace
  if (reinit && solution->HasComponent("face")) {
    EnforceConstraints(t_new, solution);

    if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
      VV_PrintHeadExtrema(*solution);
    }
  }

  // initialization
  if (num_itrs_ == 0) {
    auto udot = Teuchos::rcp(new TreeVector(*soln_));
    udot->PutScalar(0.0);
    bdf1_dae_->SetInitialState(t_old, soln_, udot);

    UpdatePreconditioner(t_old, soln_, dt_);
    num_itrs_++;
  }

  // trying to make a step
  bool failed(false);
  failed = bdf1_dae_->TimeStep(dt_, dt_next_, soln_);
  if (failed) {
    dt_ = dt_next_;

    // revover the original primary solution, pressure
    for (auto it = copies.begin(); it != copies.end(); ++it) {
      S_->GetW<CV_t>(it->first, passwd_) = it->second;

      Teuchos::OSTab tab = vo_->getOSTab();
      *vo_->os() << "Reverted field \"" << it->first << "\"" << std::endl;
    }
    archive.Restore("");
    pressure_eval_->SetChanged();

    return failed;
  }

  // commit solution (should we do it here ?)
  bdf1_dae_->CommitSolution(dt_, soln_);
  pressure_eval_->SetChanged();

  dt_tuple times(t_old, dt_);
  dT_history_.push_back(times);
  num_itrs_++;

  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
    VV_ReportWaterBalance(S_.ptr());
    VV_ReportMultiscale();
  }
  if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
    VV_ReportSeepageOutflow(S_.ptr(), dt_);
  }

  dt_ = dt_next_;
  
  return failed;
}


/* ******************************************************************
* Save internal data needed by time integration. Calculate temporarily
* the Darcy flux.
****************************************************************** */
void Richards_PK::CommitStep(double t_old, double t_new, const Tag& tag)
{
  // calculate Darcy mass flux
  auto vol_flowrate = S_->GetPtrW<CV_t>(vol_flowrate_key_, Tags::DEFAULT, passwd_);

  if (coupled_to_matrix_ || flow_on_manifold_) {
    op_matrix_diff_->UpdateFluxManifold(solution.ptr(), vol_flowrate.ptr());
    VV_FractureConservationLaw();
  } else {
    op_matrix_diff_->UpdateFlux(solution.ptr(), vol_flowrate.ptr());
  }

  Epetra_MultiVector& flowrate = *vol_flowrate->ViewComponent("face", true);
  flowrate.Scale(1.0 / molar_rho_);
  *vol_flowrate_copy->ViewComponent("face", true) = flowrate;

  // update time derivative
  *pdot_cells_prev = *pdot_cells;

  dt_ = dt_next_;
}


/* ******************************************************************
* Returns either known pressure face value or calculates it using
* the two-point flux approximation (FV) scheme.
****************************************************************** */
double Richards_PK::BoundaryFaceValue(int f, const CompositeVector& u)
{
  double face_value;

  if (u.HasComponent("face")) {
    const Epetra_MultiVector& u_face = *u.ViewComponent("face");
    face_value = u_face[0][f];
  } else {
    int c = AmanziMesh::getFaceOnBoundaryInternalCell(*mesh_, f);
    face_value = DeriveBoundaryFaceValue(f, u, wrm_->second[(*wrm_->first)[c]]);
  }
  return face_value;
}


/* ******************************************************************
* Calculates pressure value on the boundary using the two-point flux 
* approximation (FV) scheme.
****************************************************************** */
double Richards_PK::DeriveBoundaryFaceValue(
    int f, const CompositeVector& u, Teuchos::RCP<const WRM> wrm_model) 
{
  const std::vector<int>& bc_model = op_bc_->bc_model();
  const std::vector<double>& bc_value = op_bc_->bc_value();

  if (bc_model[f] == Operators::OPERATOR_BC_DIRICHLET) {
    return bc_value[f];
  } else {
    const auto& mu_cell = *S_->Get<CV_t>(viscosity_liquid_key_).ViewComponent("cell");
    const auto& u_cell = *u.ViewComponent("cell");
    AmanziMesh::Entity_ID_List cells;
    mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
    int c = cells[0];

    double pc_shift(atm_pressure_);   
    double trans_f = op_matrix_diff_->ComputeTransmissibility(f);
    double g_f = op_matrix_diff_->ComputeGravityFlux(f);
    double lmd = u_cell[0][c];
    int dir;
    mesh_->face_normal(f, false, c, &dir);
    double bnd_flux = dir * bc_value[f] / (molar_rho_ / mu_cell[0][c]);

    double max_val(atm_pressure_), min_val;
    if (bnd_flux <= 0.0) {
      min_val = u_cell[0][c];
    } else {
      min_val= u_cell[0][c] + (g_f - bnd_flux) / (dir * trans_f);
    }
    double eps = std::max(1.0e-4 * std::abs(bnd_flux), 1.0e-8);

    const KRelFn func = &WRM::k_relative;       
    Amanzi::BoundaryFaceSolver<WRM> bnd_solver(trans_f, g_f, u_cell[0][c], lmd, bnd_flux, dir, pc_shift, 
                                               min_val, max_val, eps, wrm_model, func);
    lmd = bnd_solver.FaceValue();

    return lmd;      
  }
}


/* ******************************************************************
* This is strange.
****************************************************************** */
void Richards_PK::VV_ReportMultiscale()
{
  if (multiscale_porosity_ && ms_calls_ && 
      vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "multiscale: NS:" << double(ms_itrs_) / ms_calls_ << std::endl;
  }
}


/* ******************************************************************
* This is strange.
****************************************************************** */
void Richards_PK::CalculateDiagnostics(const Tag& tag) {
  UpdateLocalFields_(S_.ptr());
}


/* ******************************************************************
* Return a pointer to a local operator
****************************************************************** */
Teuchos::RCP<Operators::Operator> Richards_PK::my_operator(
    const Operators::OperatorType& type)
{
  if (type == Operators::OPERATOR_MATRIX) return op_matrix_;
  else if (type == Operators::OPERATOR_PRECONDITIONER_RAW) return op_preconditioner_;
  return Teuchos::null;
}

}  // namespace Flow
}  // namespace Amanzi

