/*
  Transport PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <set>
#include <string>
#include <vector>

#include "Teuchos_RCP.hpp"

#include "errors.hh"
#include "MultiFunction.hh"
#include "GMVMesh.hh"

#include "Mesh.hh"
#include "Transport_PK.hh"

namespace Amanzi {
namespace Transport {

/* ******************************************************************
* Inialization of various transport structures.
****************************************************************** */
void Transport_PK::InitializeAll_()
{
  Teuchos::OSTab tab = vo_->getOSTab();

  // global transport parameters
  cfl_ = tp_list_->get<double>("cfl", 1.0);

  std::string name = tp_list_->get<std::string>("method", "muscl");
  method_ = (name == "fct") ? Method_t::FCT : Method_t::MUSCL;

  spatial_disc_order = tp_list_->get<int>("spatial discretization order", 1);
  temporal_disc_order = tp_list_->get<int>("temporal discretization order", 1);

  num_aqueous = tp_list_->get<int>("number of aqueous components", component_names_.size());
  num_gaseous = tp_list_->get<int>("number of gaseous components", 0);

  // transport dispersion (default is none)
  dispersion_solver = tp_list_->get<std::string>("solver", "missing");
  dispersion_preconditioner = tp_list_->get<std::string>("preconditioner", "identity");

  if (tp_list_->isSublist("material properties")) {
    Teuchos::ParameterList& dlist = tp_list_->sublist("material properties");

    int nblocks = 0; 
    for (auto i = dlist.begin(); i != dlist.end(); i++) {
      if (dlist.isSublist(dlist.name(i))) nblocks++;
    }

    mat_properties_.resize(nblocks);

    int iblock = 0;
    for (auto i = dlist.begin(); i != dlist.end(); i++) {
      if (dlist.isSublist(dlist.name(i))) {
        mat_properties_[iblock] = Teuchos::rcp(new MaterialProperties());

        Teuchos::ParameterList& model_list = dlist.sublist(dlist.name(i));

        mat_properties_[iblock]->tau[0] = model_list.get<double>("aqueous tortuosity", 0.0);
        mat_properties_[iblock]->tau[1] = model_list.get<double>("gaseous tortuosity", 0.0);
        mat_properties_[iblock]->regions = model_list.get<Teuchos::Array<std::string> >("regions").toVector();
        iblock++;
      }
    }
  }

  // transport diffusion (default is none)
  diffusion_phase_.resize(TRANSPORT_NUMBER_PHASES, Teuchos::null);

  if (tp_list_->isSublist("molecular diffusion")) {
    Teuchos::ParameterList& dlist = tp_list_->sublist("molecular diffusion");
    if (dlist.isParameter("aqueous names")) { 
      diffusion_phase_[0] = Teuchos::rcp(new DiffusionPhase());
      diffusion_phase_[0]->names() = dlist.get<Teuchos::Array<std::string> >("aqueous names").toVector();
      diffusion_phase_[0]->values() = dlist.get<Teuchos::Array<double> >("aqueous values").toVector();
    }

    if (dlist.isParameter("gaseous names")) { 
      diffusion_phase_[1] = Teuchos::rcp(new DiffusionPhase());
      diffusion_phase_[1]->names() = dlist.get<Teuchos::Array<std::string> >("gaseous names").toVector();
      diffusion_phase_[1]->values() = dlist.get<Teuchos::Array<double> >("gaseous values").toVector();
    }
  }

  // do we really have mechanical dispersion?
  flag_dispersion_ = false;
  if (tp_list_->isSublist("material properties")) {
    Teuchos::RCP<Teuchos::ParameterList>
        mdm_list = Teuchos::sublist(tp_list_, "material properties");
    mdm_ = CreateMDMPartition(mesh_, mdm_list, flag_dispersion_);
    if (flag_dispersion_) CalculateAxiSymmetryDirection();
  }

  // do we really have molecular diffusion?
  flag_diffusion_ = false;
  for (int i = 0; i < 2; i++) {
    if (diffusion_phase_[i] != Teuchos::null) {
      if (diffusion_phase_[i]->values().size() != 0) flag_diffusion_ = true;
    }
  }
  if (flag_diffusion_) {
    // no molecular diffusion if all tortuosities are zero.
    double tau(0.0);
    for (int i = 0; i < mat_properties_.size(); i++) {
      tau += mat_properties_[i]->tau[0] + mat_properties_[i]->tau[1];
    }
    if (tau == 0.0) flag_diffusion_ = false;
  }

  use_dispersion_ &= flag_dispersion_ || flag_diffusion_;

  // statistics of solutes
  if (tp_list_->isParameter("runtime diagnostics: solute names")) {
    runtime_solutes_ = tp_list_->get<Teuchos::Array<std::string> >("runtime diagnostics: solute names").toVector();
  } else {
    runtime_solutes_.push_back(component_names_[0]);
    if (num_gaseous > 0) runtime_solutes_.push_back(component_names_[num_aqueous]);
  }
  mass_solutes_exact_.assign(num_aqueous + num_gaseous, 0.0);
  mass_solutes_source_.assign(num_aqueous + num_gaseous, 0.0);

  if (tp_list_->isParameter("runtime diagnostics: regions")) {
    runtime_regions_ = tp_list_->get<Teuchos::Array<std::string> >("runtime diagnostics: regions").toVector();
  }

  internal_tests_ = tp_list_->get<bool>("enable internal tests", false);
  internal_tests_tol_ = tp_list_->get<double>("internal tests tolerance", TRANSPORT_CONCENTRATION_OVERSHOOT);
  dt_debug_ = tp_list_->get<double>("maximum time step", TRANSPORT_LARGE_TIME_STEP);

  if (spatial_disc_order < 1 || spatial_disc_order > 2 ||
     temporal_disc_order < 1 || temporal_disc_order > 4) {
    Errors::Message msg;
    msg << "TransportPK: unsupported combination of spatial or temporal discretization orders.\n";
    Exceptions::amanzi_throw(msg);  
  }
}


/* ****************************************************************
* Find place of the given component in a multivector.
**************************************************************** */
int Transport_PK::FindComponentNumber(const std::string component_name)
{
  int ncomponents = component_names_.size();
  for (int i = 0; i < ncomponents; i++) {
    if (component_names_[i] == component_name) return i;
  } 
  return -1;
}

}  // namespace Transport
}  // namespace Amanzi

