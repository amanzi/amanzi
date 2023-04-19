/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Flow PK

*/

#include <string>

// Amanzi
#include "OperatorDefs.hh"
#include "SolverNewton.hh"

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
  Ka_ = sublist.get<double>("absolute permeability");
  tol_ = plist.get<double>("tolerance", FLOW_DPM_NEWTON_TOLERANCE);
  clip_ = plist.get<double>("clipping factor", FLOW_DPM_NEWTON_CLIPPING_FACTOR);
  atm_pressure_ = plist.get<double>("atmospheric pressure", FLOW_PRESSURE_ATMOSPHERIC);
}


/* ******************************************************************
* Compute water storage.
****************************************************************** */
double
MultiscaleFlowPorosity_GDPM::ComputeField(double phi, double n_l, double prm)
{
  double pc = atm_pressure_ - prm;
  return wrm_->saturation(pc) * phi * n_l;
}


/* ******************************************************************
* Main capability: cell-based Newton solver. It returns water storage
* and pressure in the matrix. max_itrs is input/output parameter.
****************************************************************** */
WhetStone::DenseVector
MultiscaleFlowPorosity_GDPM::WaterContentMatrix(double prf0,
                                                WhetStone::DenseVector& prm,
                                                WhetStone::DenseVector& wcm0,
                                                double dt,
                                                double phi,
                                                double n_l,
                                                double mu_l,
                                                int& max_itrs)
{
  dt_ = dt;
  phi_ = phi;
  nl_ = n_l;
  mu_ = mu_l;

  bcl_ = prf0;
  wcm0_ = wcm0;
  wcm1_ = wcm0;

  // make uniform mesh inside the matrix
  auto mesh = std::make_shared<WhetStone::DenseVector>(WhetStone::DenseVector(matrix_nodes_ + 1));
  double h = depth_ / matrix_nodes_;
  for (int i = 0; i < matrix_nodes_ + 1; ++i) (*mesh)(i) = h * i;

  // initialize diffusion operator
  auto kr = std::make_shared<WhetStone::DenseVector>(WhetStone::DenseVector(matrix_nodes_));
  auto dkdp = std::make_shared<WhetStone::DenseVector>(WhetStone::DenseVector(matrix_nodes_));

  op_diff_ = std::make_shared<Operators::Mini_Diffusion1D>();
  op_diff_->Init(mesh);
  op_diff_->Setup(Ka_);
  op_diff_->Setup(kr, dkdp);

  // create nonlinear problem
  Teuchos::RCP<SolverFnBase<WhetStone::DenseVector>> fn = Teuchos::rcpFromRef(*this);
  auto sol = Teuchos::rcpFromRef(prm);

  // create the solver
  Teuchos::ParameterList plist;
  plist.set<double>("nonlinear tolerance", 1.0e-7)
       .set<int>("limit iterations", max_itrs)
       .sublist("verbose object").set<std::string>("verbosity level", "low");

  Amanzi::AmanziSolvers::SolverNewton<WhetStone::DenseVector, int> newton(plist);
  newton.Init(fn, 1);

  // solve the problem
  newton.Solve(sol);
  max_itrs = newton.num_itrs();

  return wcm1_;
}


/* ******************************************************************
* Nonlinear residual for Newton-type solvers.
****************************************************************** */
void
MultiscaleFlowPorosity_GDPM::Residual(const Teuchos::RCP<WhetStone::DenseVector>& u,
                                      const Teuchos::RCP<WhetStone::DenseVector>& f)
{
  double factor = nl_ / mu_;
  for (int i = 0; i < u->NumRows(); ++i) { 
    double pc = atm_pressure_ - (*u)(i);
    op_diff_->k()(i) = factor * wrm_->k_relative(pc);
  }

  op_diff_->rhs().PutScalar(0.0);
  op_diff_->UpdateMatrices();
  op_diff_->ApplyBCs(bcl_, Operators::OPERATOR_BC_DIRICHLET, 0.0, Operators::OPERATOR_BC_NEUMANN);

  op_diff_->Apply(*u, *f);
  f->Update(-1.0, op_diff_->rhs(), 1.0);

  for (int i = 0; i < u->NumRows(); ++i) { 
    double h = op_diff_->mesh_cell_volume(i);
    double pc = atm_pressure_ - (*u)(i);
    wcm1_(i) = nl_ * phi_ * wrm_->saturation(pc);
    (*f)(i) += (wcm1_(i) - wcm0_(i)) * h / dt_;
  }
}


/* ******************************************************************
* Populate preconditioner's data.
****************************************************************** */
void
MultiscaleFlowPorosity_GDPM::UpdatePreconditioner(
  const Teuchos::RCP<const WhetStone::DenseVector>& u)
{
  double factor = nl_ / mu_;
  for (int i = 0; i < u->NumRows(); ++i) {
    double pc = atm_pressure_ - (*u)(i);
    op_diff_->k()(i) = factor * wrm_->k_relative(pc);
    op_diff_->dkdp()(i) = -factor * wrm_->dKdPc(pc);
  }

  op_diff_->UpdateJacobian(
    *u, bcl_, Operators::OPERATOR_BC_DIRICHLET, 0.0, Operators::OPERATOR_BC_NEUMANN);
  // op_diff_->UpdateMatrices();
  // op_diff_->ApplyBCs(bcl_, Operators::OPERATOR_BC_DIRICHLET, 0.0, Operators::OPERATOR_BC_NEUMANN);

  WhetStone::DenseVector dwcm(wcm0_);
  for (int i = 0; i < u->NumRows(); ++i) { 
    double h = op_diff_->mesh_cell_volume(i);
    double pc = atm_pressure_ - (*u)(i);
    dwcm(i) = -wrm_->dSdPc(pc) * nl_ * phi_ * h / dt_;
  }
  op_diff_->AddAccumulationTerm(dwcm, false);
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


/********************************************************************
* Modifies nonlinear update du based on the maximum allowed change
* of saturation.
****************************************************************** */
AmanziSolvers::FnBaseDefs::ModifyCorrectionResult 
MultiscaleFlowPorosity_GDPM::ModifyCorrection(const Teuchos::RCP<const WhetStone::DenseVector>& r,
                                              const Teuchos::RCP<const WhetStone::DenseVector>& u,
                                              const Teuchos::RCP<WhetStone::DenseVector>& du)
{
  int ncells = u->NumRows();
  double s = 0.0;
  for (int i = 0; i < ncells; ++i) {
    s = std::max(s, std::fabs((*du)(i)) / (std::fabs((*u)(i)) + atm_pressure_));
  }

  if (s > clip_) {
    (*du) *= std::min(clip_, 1.0 / s);
    return AmanziSolvers::FnBaseDefs::CORRECTION_MODIFIED;
  }

  return AmanziSolvers::FnBaseDefs::CORRECTION_NOT_MODIFIED;
}

} // namespace Flow
} // namespace Amanzi
