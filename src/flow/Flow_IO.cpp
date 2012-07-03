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
    Teuchos::ParameterList& tmp_list = list.sublist(name);
    ti_specs.atol = tmp_list.get<double>("absolute error tolerance", FLOW_TI_ABSOLUTE_TOLERANCE);
    ti_specs.rtol = tmp_list.get<double>("relative error tolerance", FLOW_TI_RELATIVE_TOLERANCE);
    ti_specs.residual_tol = tmp_list.get<double>("convergence tolerance", FLOW_TI_NONLINEAR_RESIDUAL_TOLERANCE);
    ti_specs.max_itrs = tmp_list.get<int>("maximum number of iterations", FLOW_TI_MAX_ITERATIONS);

    ti_specs.T0 = tmp_list.get<double>("start time", 0.0);
    ti_specs.T1 = tmp_list.get<double>("end time", 0.0);
    ti_specs.dTfactor = tmp_list.get<double>("time step increase factor", 1.0);
    ti_specs.dT0 = tmp_list.get<double>("initial time step", AmanziFlow::FLOW_INITIAL_DT);
    ti_specs.dTmax = tmp_list.get<double>("maximum time step", AmanziFlow::FLOW_MAXIMUM_DT);

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
  if (name == "Trilinos ML") {
    *preconditioner = FLOW_PRECONDITIONER_TRILINOS_ML;
  } else if (name == "Hypre AMG") {
    *preconditioner = FLOW_PRECONDITIONER_HYPRE_AMG;
#ifndef HAVE_HYPRE_API
    Errors::Message msg;
    msg << "\nFlow PK: Hypre TPL has not been activated.";
    Exceptions::amanzi_throw(msg);   
#endif
  }
}

}  // namespace AmanziFlow
}  // namespace Amanzi

