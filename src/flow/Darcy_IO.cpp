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
  ProcessStringVerbosity(verbosity_name, &verbosity);
 
  atm_pressure = dp_list_.get<double>("atmospheric pressure", 101325.0);

  // Create the BC objects.
  Teuchos::RCP<Teuchos::ParameterList> bc_list = Teuchos::rcpFromRef(dp_list_.sublist("boundary conditions", true));
  FlowBCFactory bc_factory(mesh_, bc_list);

  bc_pressure = bc_factory.createPressure();
  bc_head = bc_factory.createStaticHead(atm_pressure, rho_, gravity_);
  bc_flux = bc_factory.createMassFlux();
  bc_seepage = bc_factory.createSeepageFace();

  ValidateBoundaryConditions(bc_pressure, bc_head, bc_flux);

  double time = T_physics;
  bc_pressure->Compute(time);
  bc_head->Compute(time);
  bc_flux->Compute(time);
  bc_seepage->Compute(time);

  // Create the source object if any
  if (dp_list_.isSublist("source terms")) {
    string distribution_method_name = dp_list_.get<string>("source and sink distribution method", "none");
    ProcessStringSourceDistribution(distribution_method_name, &src_sink_distribution); 

    Teuchos::RCP<Teuchos::ParameterList> src_list = Teuchos::rcpFromRef(dp_list_.sublist("source terms", true));
    FlowSourceFactory src_factory(mesh_, src_list);
    src_sink = src_factory.createSource();

    src_sink->Compute(time);  // remove it from this routine (lipnikov@lanl.gov)
  } else {
    src_sink = NULL;
  }

  // discretization method
  string mfd3d_method_name = dp_list_.get<string>("discretization method", "optimized mfd");
  ProcessStringMFD3D(mfd3d_method_name, &mfd3d_method); 

  // Time integrator for period I, temporary called steady state time integrator
  Teuchos::ParameterList& sss_list = dp_list_.sublist("steady state time integrator");

  string ti_method_name = sss_list.get<string>("time integration method", "none");
  ProcessSublistTimeIntegration(sss_list, ti_method_name, ti_specs_sss);

  if (sss_list.isParameter("preconditioner")) {
    preconditioner_name_sss_ = sss_list.get<string>("preconditioner");
  } else {
    msg << "Darcy PK: steady state time integrator does not define a preconditioner.";
    Exceptions::amanzi_throw(msg);
  }

  if (! preconditioner_list_.isSublist(preconditioner_name_sss_)) {
    msg << "Darcy PK: steady state preconditioner does not exist.";
    Exceptions::amanzi_throw(msg);
  }

  std::string linear_solver_name;
  if (sss_list.isParameter("linear solver")) {
    linear_solver_name = sss_list.get<string>("linear solver");
  } else {
    msg << "Darcy PK: steady state time integrator does not define <linear solver>.";
    Exceptions::amanzi_throw(msg);
  }
  ProcessStringLinearSolver(linear_solver_name, &max_itrs_sss, &convergence_tol_sss);
}


/* ****************************************************************
* Process string for the linear solver.
**************************************************************** */
void Darcy_PK::ProcessStringLinearSolver(
    const std::string name, int* max_itrs, double* convergence_tol)
{
  Errors::Message msg;

  if (! solver_list_.isSublist(name)) {
    msg << "Richards PK: steady state linear solver does not exist.";
    Exceptions::amanzi_throw(msg);
  }

  Teuchos::ParameterList& tmp_list = solver_list_.sublist(name);
  *max_itrs = tmp_list.get<int>("maximum number of iterations", 100);
  *convergence_tol = tmp_list.get<double>("error tolerance", 1e-12);
}


/* ****************************************************************
* Process string for the linear solver.
**************************************************************** */
void Darcy_PK::ProcessStringSourceDistribution(const std::string name, int* method)
{
  Errors::Message msg;
  if (name == "none") {
    *method = AmanziFlow::FLOW_SOURCE_DISTRIBUTION_NONE;
  } else if (name == "volume") {
    *method = AmanziFlow::FLOW_SOURCE_DISTRIBUTION_VOLUME;
  } else if (name == "premeability") {
    *method = AmanziFlow::FLOW_SOURCE_DISTRIBUTION_PERMEABILITY;
  } else {
    msg << "Darcy PK: unknown source normalization method has been specified.";
    Exceptions::amanzi_throw(msg);
  }
}

}  // namespace AmanziFlow
}  // namespace Amanzi


