/*
  Flow PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <string>

// Amanzi
#include "OperatorDefs.hh"

#include "FlowDefs.hh"
#include "MultiscaleFlowPorosity_GDPM.hh"
#include "WRMFactory.hh"

namespace Amanzi {
namespace Flow {

/* ******************************************************************
* This model is minor extension of the WRM.
****************************************************************** */
MultiscaleFlowPorosity_GDPM::MultiscaleFlowPorosity_GDPM(Teuchos::ParameterList& plist)
{
  WRMFactory factory;
  wrm_ = factory.Create(plist);

  auto& sublist = plist.sublist("generalized dual porosity parameters");
  matrix_nodes_ = sublist.get<int>("number of matrix nodes");

  // depth is defined for each matrix block as A_m / V_m, so in general,
  // it depends on geometry
  depth_ = sublist.get<double>("matrix depth");
  tau_ = sublist.get<double>("matrix tortuosity");
  Ka_ = sublist.get<double>("absolute permeability");
  tol_ = plist.get<double>("tolerance", FLOW_DPM_NEWTON_TOLERANCE);
}


/* ******************************************************************
* Compute water storage.
****************************************************************** */
double
MultiscaleFlowPorosity_GDPM::ComputeField(double phi, double n_l, double pcm)
{
  return wrm_->saturation(pcm) * phi * n_l;
}


/* ******************************************************************
* Main capability: cell-based Newton solver. It returns water storage 
* and pressure in the matrix. max_itrs is input/output parameter.
****************************************************************** */
double
MultiscaleFlowPorosity_GDPM::WaterContentMatrix(double pcf0,
                                                WhetStone::DenseVector& pcm,
                                                double wcm0,
                                                double dt,
                                                double phi,
                                                double n_l,
                                                int& max_itrs)
{
  Operators::Mini_Diffusion1D op_diff;

  // make uniform mesh inside the matrix
  auto mesh = std::make_shared<WhetStone::DenseVector>(WhetStone::DenseVector(matrix_nodes_ + 1));
  double h = depth_ / matrix_nodes_;
  for (int i = 0; i < matrix_nodes_ + 1; ++i) (*mesh)(i) = h * i;

  // initialize diffusion operator
  auto kr = std::make_shared<WhetStone::DenseVector>(WhetStone::DenseVector(matrix_nodes_));
  auto dkdp = std::make_shared<WhetStone::DenseVector>(WhetStone::DenseVector(matrix_nodes_));

  op_diff.Init(mesh);
  op_diff.Setup(Ka_);
  op_diff.Setup(kr, dkdp);

  // create nonlinear problem
  Teuchos::RCP<SolverFnBase<WhetStone::DenseVector>> fn = Teuchos::rcp(this);
  auto sol = Teuchos::rcp(new WhetStone::DenseVector(matrix_nodes_));

  // create the solver
  Teuchos::ParameterList plist;
  plist.set<double>("nonlinear tolerance", 1.0e-7);
  plist.sublist("verbose object").set<std::string>("verbosity level", "high");

  Amanzi::AmanziSolvers::SolverNewton<WhetStone::DenseVector, int> newton(plist);
  newton.Init(fn, 1);

  // solve the problem
  newton.Solve(sol);

  return wrm_->saturation(pcm(0)) * phi * n_l;
}


/* ******************************************************************
* Nonlinear residual for Newton-type solvers.
****************************************************************** */
void
MultiscaleFlowPorosity_GDPM::Residual(const Teuchos::RCP<WhetStone::DenseVector>& u,
                                      const Teuchos::RCP<WhetStone::DenseVector>& f)
{
  for (int i = 0; i < u->NumRows(); ++i) { op_->k()(i) = 1.0 + (*u)(i) * (*u)(i); }

  op_->UpdateMatrices();
  op_->ApplyBCs(bcl_, Operators::OPERATOR_BC_DIRICHLET, 0.0, Operators::OPERATOR_BC_NEUMANN);

  op_->Apply(*u, *f);
  f->Update(-1.0, op_->rhs(), 1.0);
}


/* ******************************************************************
* Populate preconditioner's data.
****************************************************************** */
void
MultiscaleFlowPorosity_GDPM::UpdatePreconditioner(
  const Teuchos::RCP<const WhetStone::DenseVector>& u)
{
  for (int i = 0; i < u->NumRows(); ++i) {
    op_->k()(i) = 1.0 + (*u)(i) * (*u)(i);
    op_->dkdp()(i) = 2.0 * (*u)(i);
  }
  op_->UpdateJacobian(
    *u, bcl_, Operators::OPERATOR_BC_DIRICHLET, 0.0, Operators::OPERATOR_BC_NEUMANN);
}


/* ******************************************************************
* Error estimate for convergence criteria.
****************************************************************** */
double
MultiscaleFlowPorosity_GDPM::ErrorNorm(const Teuchos::RCP<const WhetStone::DenseVector>& u,
                                       const Teuchos::RCP<const WhetStone::DenseVector>& du)
{
  double tmp;
  du->NormInf(&tmp);
  return tmp;
}

} // namespace Flow
} // namespace Amanzi
