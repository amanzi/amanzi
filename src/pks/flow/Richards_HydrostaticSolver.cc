/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Flow PK

  Hydrostatic solver for compressible water at high P/T values.
*/

#include "InverseFactory.hh"
#include "MeshAlgorithms.hh"
#include "OperatorDefs.hh"
#include "PDE_Diffusion.hh"
#include "PDE_DiffusionFactory.hh"
#include "PK_DomainFunctionSimple.hh"
#include "RemapUtils.hh"

#include "FlowBoundaryFunction.hh"
#include "Flow_SolverFnBase.hh"
#include "Richards_PK.hh"

namespace Amanzi {
namespace Flow {

/* ******************************************************************
* Solve nonlinear problem div(rho q) = 0 using Dirichlet boundary
* condition at one point.
****************************************************************** */
void
Richards_PK::SolveHydrostaticProblem(Teuchos::ParameterList& plist, const Teuchos::RCP<TreeVector>& u)
{
  // create a new set of control variables: error computing method, 
  // boundary conditions, and source terms
  int error_copy = error_control_;
  auto srcs_copy = srcs;
  auto bcs_copy = bcs_;

  error_control_ = FLOW_TI_ERROR_CONTROL_PRESSURE;
  srcs.clear();
  bcs_.clear();

  int f;
  std::vector<double> value;

  for (int i = 0; i < bcs_copy.size(); ++i) {
    if (bcs_copy[i]->get_bc_name() == "pressure") {
      bcs_.push_back(bcs_copy[i]);
      break; 
    }
  }

  auto solver_fn = Teuchos::rcp(new Flow_SolverFnBase<TreeVector>(plist, this));
  AmanziSolvers::SolverFactory<TreeVector, TreeVectorSpace> factory;
  auto solver = factory.Create(plist);

  solver->Init(solver_fn, u->Map());
  solver->Solve(u);

  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    int num_itrs = solver->num_itrs();
    int code = solver->returned_code();
    double pnorm;
    u->Norm2(&pnorm);

    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "nonlinear saturated solver (NKA): ||p,lambda||=" << pnorm
               << " itr=" << num_itrs << " code=" << code << std::endl;
  }

  // restore original variables
  error_control_ = error_copy;
  srcs = srcs_copy;
  bcs_ = bcs_copy;
}

} // namespace Flow
} // namespace Amanzi
