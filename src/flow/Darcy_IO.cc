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

#include "Flow_BC_Factory.hh"
#include "Flow_SourceFactory.hh"
#include "Darcy_PK.hh"

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

  bc_pressure = bc_factory.CreatePressure(bc_submodel);
  bc_head = bc_factory.CreateStaticHead(atm_pressure, rho_, gravity_, bc_submodel);
  bc_flux = bc_factory.CreateMassFlux(bc_submodel);
  bc_seepage = bc_factory.CreateSeepageFace(bc_submodel);

  ValidateBoundaryConditions(bc_pressure, bc_head, bc_flux);
  ProcessStaticBCsubmodels(bc_submodel, rainfall_factor);

  // Create the source object if any
  if (dp_list_.isSublist("source terms")) {
    string distribution_method_name = dp_list_.get<string>("source and sink distribution method", "none");
    ProcessStringSourceDistribution(distribution_method_name, &src_sink_distribution);

    Teuchos::RCP<Teuchos::ParameterList> src_list = Teuchos::rcpFromRef(dp_list_.sublist("source terms", true));
    FlowSourceFactory src_factory(mesh_, src_list);
    src_sink = src_factory.createSource();
    src_sink_distribution = src_sink->CollectActionsList();
  }

  // discretization method
  string mfd3d_method_name = dp_list_.get<string>("discretization method", "optimized mfd");
  ProcessStringMFD3D(mfd3d_method_name, &mfd3d_method); 

  // Time integrator for period I, called steady-state time integrator
  if (dp_list_.isSublist("steady state time integrator")) {
    Teuchos::ParameterList& sss_list = dp_list_.sublist("steady state time integrator");

    string ti_method_name = sss_list.get<string>("time integration method", "none");
    ProcessStringTimeIntegration(ti_method_name, &ti_specs_sss.ti_method);
    ProcessSublistTimeIntegration(sss_list, ti_method_name, ti_specs_sss);
    ti_specs_sss.ti_method_name = "steady state time integrator";

    if (sss_list.isParameter("preconditioner")) {
      ti_specs_sss.preconditioner_name = sss_list.get<string>("preconditioner");
      ProcessStringPreconditioner(ti_specs_sss.preconditioner_name, &ti_specs_sss.preconditioner_method);
    } else {
      msg << "Flow PK: steady state time integrator does not define a preconditioner.";
      Exceptions::amanzi_throw(msg);
    }

    std::string linear_solver_name = FindStringLinearSolver(sss_list, solver_list_);
    LinearSolver_Specs& ls_specs = ti_specs_sss.ls_specs;
    ProcessStringLinearSolver(linear_solver_name, &ls_specs.max_itrs, &ls_specs.convergence_tol);
  } else {  // fills-in the important defaults (See InputParserIS_defaults.hh)
    std::string linear_solver_name("AztecOO");  // Must equal to ST_SOLVER
    LinearSolver_Specs& ls_specs = ti_specs_sss.ls_specs;
    ProcessStringLinearSolver(linear_solver_name, &ls_specs.max_itrs, &ls_specs.convergence_tol);   
  }

  if (dp_list_.isSublist("transient time integrator")) {
    Teuchos::ParameterList& trs_list = dp_list_.sublist("transient time integrator");

    string ti_method_name = trs_list.get<string>("time integration method", "none");
    ProcessStringTimeIntegration(ti_method_name, &ti_specs_trs.ti_method);
    ProcessSublistTimeIntegration(trs_list, ti_method_name, ti_specs_trs);
    ti_specs_trs.ti_method_name = "transient time integrator";

    if (trs_list.isParameter("preconditioner")) {
      ti_specs_trs.preconditioner_name = trs_list.get<string>("preconditioner");
      ProcessStringPreconditioner(ti_specs_trs.preconditioner_name, &ti_specs_trs.preconditioner_method);
    } else {
      msg << "Flow PK: steady state time integrator does not define a preconditioner.";
      Exceptions::amanzi_throw(msg);
    }

    std::string linear_solver_name = FindStringLinearSolver(trs_list, solver_list_);
    LinearSolver_Specs& ls_specs = ti_specs_trs.ls_specs;
    ProcessStringLinearSolver(linear_solver_name, &ls_specs.max_itrs, &ls_specs.convergence_tol);
  }
}


/* ****************************************************************
* Process string for the linear solver.
**************************************************************** */
void Darcy_PK::ProcessStringLinearSolver(
    const std::string name, int* max_itrs, double* convergence_tol)
{
  Errors::Message msg;

  if (! solver_list_.isSublist(name)) {
    msg << "Flow PK: steady state linear solver does not exist.";
    Exceptions::amanzi_throw(msg);
  }

  Teuchos::ParameterList& tmp_list = solver_list_.sublist(name);
  *max_itrs = tmp_list.get<int>("maximum number of iterations", 100);
  *convergence_tol = tmp_list.get<double>("error tolerance", 1e-12);
}


/* ******************************************************************
* Prints information about status of this PK.                                                     
****************************************************************** */
void Darcy_PK::PrintStatistics() const
{
  if (MyPID == 0) {
    cout << "Flow PK:" << endl;
    cout << "    Verbosity level = " << verbosity << endl;
  }
}

}  // namespace AmanziFlow
}  // namespace Amanzi


