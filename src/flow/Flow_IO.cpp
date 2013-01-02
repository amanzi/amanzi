/*
This is the flow component of the Amanzi code. 

Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided in the top-level COPYRIGHT file.

Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <vector>

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "errors.hh"
#include "tabular-function.hh"

#include "Flow_PK.hpp"

#include "Ifpack.h"

namespace Amanzi {
namespace AmanziFlow {

/* ****************************************************************
* Process time integration sublist.
**************************************************************** */
void Flow_PK::ProcessSublistTimeIntegration(
    Teuchos::ParameterList& list, const std::string name, TI_Specs& ti_specs)
{
  Errors::Message msg;

  if (list.isSublist(name)) {
    ti_specs.initialize_with_darcy = list.get<bool>("initialize with darcy", false);
    ti_specs.clip_saturation = list.get<double>("clipping saturation value", -1.0);
    ti_specs.clip_pressure = list.get<double>("clipping pressure value", -1e+10);

    if (list.isParameter("enforce pressure-lambda constraints"))
        ti_specs.pressure_lambda_constraints = list.get<bool>("enforce pressure-lambda constraints");

    ti_specs.dT_method = 0;
    std::string dT_name = list.get<string>("time stepping strategy", "simple");
    if (dT_name == "adaptive") ti_specs.dT_method = FLOW_DT_ADAPTIVE;

    Teuchos::ParameterList& tmp_list = list.sublist(name);
    ti_specs.residual_tol = tmp_list.get<double>("convergence tolerance", FLOW_TI_NONLINEAR_RESIDUAL_TOLERANCE);
    ti_specs.max_itrs = tmp_list.get<int>("maximum number of iterations", FLOW_TI_MAX_ITERATIONS);

    ti_specs.T0 = tmp_list.get<double>("start time", 0.0);
    ti_specs.T1 = tmp_list.get<double>("end time", 0.0);
    ti_specs.dTfactor = tmp_list.get<double>("time step increase factor", 1.0);
    ti_specs.dT0 = tmp_list.get<double>("initial time step", AmanziFlow::FLOW_INITIAL_DT);
    ti_specs.dTmax = tmp_list.get<double>("maximum time step", AmanziFlow::FLOW_MAXIMUM_DT);

    ti_specs.atol = tmp_list.get<double>("absolute error tolerance", 1e-3);
    ti_specs.rtol = tmp_list.get<double>("relative error tolerance", 1e-3);

  } else if (name != "none") {
    msg << "\nFlow PK: specified time integration sublist does not exist.";
    Exceptions::amanzi_throw(msg);
  }
}


/* ****************************************************************
* Process string for the discretization method.
**************************************************************** */
void Flow_PK::ProcessStringVerbosity(const std::string name, int* verbosity)
{
  verbosity_AztecOO = AZ_none;
  if (name == "none") {
    *verbosity = FLOW_VERBOSITY_NONE;
  } else if (name == "low") {
    *verbosity = FLOW_VERBOSITY_LOW;
  } else if (name == "medium") {
    *verbosity = FLOW_VERBOSITY_MEDIUM;
  } else if (name == "high") {
    *verbosity = FLOW_VERBOSITY_HIGH;
  } else if (name == "extreme") {
    *verbosity = FLOW_VERBOSITY_EXTREME;
    verbosity_AztecOO = AZ_all;
  }
}


/* ****************************************************************
* Process string for the linear solver.
**************************************************************** */
void Flow_PK::ProcessStringSourceDistribution(const std::string name, int* method)
{
  Errors::Message msg;
  if (name == "none") {
    *method = AmanziFlow::FLOW_SOURCE_DISTRIBUTION_NONE;
  } else if (name == "volume") {
    *method = AmanziFlow::FLOW_SOURCE_DISTRIBUTION_VOLUME;
  } else if (name == "permeability") {
    *method = AmanziFlow::FLOW_SOURCE_DISTRIBUTION_PERMEABILITY;
  } else {
    msg << "Flow PK: unknown source normalization method has been specified.";
    Exceptions::amanzi_throw(msg);
  }
}


/* ****************************************************************
* Process string for the discretization method.
**************************************************************** */
void Flow_PK::ProcessStringMFD3D(const std::string name, int* method)
{
  if (name == "monotone mfd") {
    *method = FLOW_MFD3D_HEXAHEDRA_MONOTONE;
  } else if (name == "support operator") {
    *method = FLOW_MFD3D_SUPPORT_OPERATOR;
  } else if (name == "two point flux approximation") {
    *method = FLOW_MFD3D_TWO_POINT_FLUX;
  } else if (name == "optimized mfd") {
    *method = FLOW_MFD3D_OPTIMIZED;
  } else if (name == "mfd") {
    *method = FLOW_MFD3D_POLYHEDRA;
  } else {
    *method = FLOW_MFD3D_POLYHEDRA;
  }
}


/* ****************************************************************
* Process string for the preconitioner.
**************************************************************** */
void Flow_PK::ProcessStringPreconditioner(const std::string name, int* preconditioner)
{
  Errors::Message msg;

  if (name == "Trilinos ML") {
    *preconditioner = FLOW_PRECONDITIONER_TRILINOS_ML;
  } else if (name == "Hypre AMG") {
    *preconditioner = FLOW_PRECONDITIONER_HYPRE_AMG;
#ifndef HAVE_HYPRE
    Errors::Message msg;
    msg << "\nFlow PK: Hypre TPL has not been activated.";
    Exceptions::amanzi_throw(msg);   
#endif
  } else if (name == "Block ILU") {
    *preconditioner = FLOW_PRECONDITIONER_TRILINOS_BLOCK_ILU;
  } else {
    msg << "\nFlow PK: specified preconditioner does not exist.";
    Exceptions::amanzi_throw(msg);
  }
}


/* ****************************************************************
* Find string for the preconditoner.
**************************************************************** */
std::string Flow_PK::FindStringLinearSolver(const Teuchos::ParameterList& list, 
                                            const Teuchos::ParameterList& solver_list)
{   
  Errors::Message msg;
  std::string name; 

  if (list.isParameter("linear solver")) {
    name = list.get<string>("linear solver");
  } else {
    msg << "Flow PK: time integrator does not define <linear solver>.";
    Exceptions::amanzi_throw(msg);
  }

  if (! solver_list.isSublist(name)) {
    msg << "Flow PK: linear solver \"" << name.c_str() << "\" does not exist.";
    Exceptions::amanzi_throw(msg);
  }
  return name;
}


/* ****************************************************************
* Find string for the preconditoner.
**************************************************************** */
void Flow_PK::OutputTimeHistory(std::vector<dt_tuple>& dT_history)
{
  if (MyPID == 0 && verbosity >= FLOW_VERBOSITY_MEDIUM) {
    printf("Flow PK: saving krel-pc curves in file flow_dt_history.txt...\n");

    char file_name[30];
    sprintf(file_name, "flow_dt_history_%d.txt", ti_phase_counter++);

    ofstream ofile;
    ofile.open(file_name);

    for (double n = 0; n < dT_history.size(); n++) {
      ofile << setprecision(10) << dT_history[n].first / FLOW_YEAR << " " << dT_history[n].second << endl;
    }
    ofile.close();
  }
}


}  // namespace AmanziFlow
}  // namespace Amanzi

