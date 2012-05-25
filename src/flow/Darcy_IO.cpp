/*
This is the flow component of the Amanzi code. 

Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided in the top-level COPYRIGHT file.

Authors: Neil Carlson (version 1) 
         Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
*/

#include <string>
#include <set>

#include "Flow_BC_Factory.hpp"
#include "Flow_Source_Factory.hpp"
#include "Darcy_PK.hpp"

namespace Amanzi {
namespace AmanziFlow {

/* ******************************************************************
* Routine processes parameter list. It needs to be called only once
* on each processor.                                                     
****************************************************************** */
void Darcy_PK::ProcessParameterList()
{
  Errors::Message msg;

  // create verbosity list if it does not exist
  if (! dp_list_.isSublist("VerboseObject")) {
    Teuchos::ParameterList verbosity_list;
    verbosity_list.set<std::string>("Verbosity Level", "none");
    dp_list_.set("VerboseObject", verbosity_list);
  }

  // extract verbosity level
  Teuchos::ParameterList verbosity_list = dp_list_.get<Teuchos::ParameterList>("VerboseObject");
  std::string verbosity_name = verbosity_list.get<std::string>("Verbosity Level");
  if (verbosity_name == "none") {
    verbosity = FLOW_VERBOSITY_NONE;
  } else if (verbosity_name == "low") {
    verbosity = FLOW_VERBOSITY_LOW;
  } else if (verbosity_name == "medium") {
    verbosity = FLOW_VERBOSITY_MEDIUM;
  } else if (verbosity_name == "high") {
    verbosity = FLOW_VERBOSITY_HIGH;
  } else if (verbosity_name == "extreme") {
    verbosity = FLOW_VERBOSITY_EXTREME;
  }

  atm_pressure = dp_list_.get<double>("atmospheric pressure");

  // Create the BC objects.
  Teuchos::RCP<Teuchos::ParameterList> bc_list = Teuchos::rcpFromRef(dp_list_.sublist("boundary conditions", true));
  FlowBCFactory bc_factory(mesh_, bc_list);

  bc_pressure = bc_factory.createPressure();
  bc_head = bc_factory.createStaticHead(atm_pressure, rho_, gravity_);
  bc_flux = bc_factory.createMassFlux();
  bc_seepage = bc_factory.createSeepageFace();

  validate_boundary_conditions(bc_pressure, bc_head, bc_flux);

  double time = T_internal;
  bc_pressure->Compute(time);
  bc_head->Compute(time);
  bc_flux->Compute(time);
  bc_seepage->Compute(time);

  // Create the source object if any
  if (dp_list_.isSublist("source terms")) {
    Teuchos::RCP<Teuchos::ParameterList> src_list = Teuchos::rcpFromRef(dp_list_.sublist("source terms", true));
    FlowSourceFactory src_factory(mesh_, src_list);
    src_sink = src_factory.createSource();

    src_sink->Compute(time);
  } else {
    src_sink = NULL;
  }

  // discretization method
  string mfd3d_method_name = dp_list_.get<string>("discretization method", "optimized mfd");
  ProcessStringMFD3D(mfd3d_method_name, &mfd3d_method); 

  // Set up internal clock.
  double T_physical = FS->get_time();
  T_internal = (standalone_mode) ? T_internal : T_physical;

  // Time integrator for period I, temporary called steady state time integrator
  Teuchos::ParameterList& sss_list = dp_list_.sublist("steady state time integrator");
  Teuchos::ParameterList& solver_list = sss_list.sublist("linear solver");

  max_itrs_sss = solver_list.get<int>("maximal number of iterations", 100);
  convergence_tol_sss = solver_list.get<double>("error tolerance", 1e-12);

  if (sss_list.isParameter("preconditioner")) {
    preconditioner_name_sss_ = sss_list.get<string>("preconditioner");
  } else {
    msg << "Darcy Problem: steady state time integrator does not define a preconditioner.";
    Exceptions::amanzi_throw(msg);
  }

  if (! preconditioner_list_.isSublist(preconditioner_name_sss_)) {
    msg << "Darcy Problem: steady state preconditioner does not exist.";
    Exceptions::amanzi_throw(msg);
  }

}

}  // namespace AmanziFlow
}  // namespace Amanzi


