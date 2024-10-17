/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  MechanicsElasticity PK

*/

#include "WhetStoneDefs.hh"
#include "PorosityEvaluator.hh"

#include "MechanicsElasticity_PK.hh"

namespace Amanzi {
namespace Mechanics {

using CV_t = CompositeVector;

/* ******************************************************************
* Calculate f(u, du/dt) = a d(s(u))/dt + A*u - rhs.
****************************************************************** */
void
MechanicsElasticity_PK::FunctionalResidual(double t_old,
                                           double t_new,
                                           Teuchos::RCP<const TreeVector> u_old,
                                           Teuchos::RCP<TreeVector> u_new,
                                           Teuchos::RCP<TreeVector> f)
{
  UpdateSourceBoundaryData(t_old, t_new);

  op_matrix_elas_->global_operator()->Init();

  // Add external forces
  auto rhs = op_matrix_->rhs();
  if (use_gravity_) AddGravityTerm(*rhs);
  if (poroelasticity_) AddPressureGradient(*rhs);
  if (thermoelasticity_) AddTemperatureGradient(*rhs);

  op_matrix_elas_->UpdateMatrices();
  op_matrix_elas_->ApplyBCs(true, true, true);

  // compute negative residual, A u - f
  op_matrix_->ComputeNegativeResidual(*u_new->Data(), *f->Data());
}


/* ******************************************************************
* Apply preconditioner inv(B) * X.
****************************************************************** */
int
MechanicsElasticity_PK::ApplyPreconditioner(Teuchos::RCP<const TreeVector> X,
                                            Teuchos::RCP<TreeVector> Y)
{
  Y->PutScalar(0.0);
  return op_matrix_->ApplyInverse(*X->Data(), *Y->Data());
}


/* ******************************************************************
* Update new preconditioner on the interval (tp-dtp, tp].
****************************************************************** */
void
MechanicsElasticity_PK::UpdatePreconditioner(double tp,
                                             Teuchos::RCP<const TreeVector> u,
                                             double dtp)
{
  // force component-wise coarsening
  /*
  auto block_indices = Teuchos::rcp(new std::vector<int>(nnodes_owned_ * dim_));
  for (int n = 0; n < nnodes_owned_ * dim_; ++n) (*block_indices)[n] = n % dim_;
  auto block_ids = std::make_pair(dim_, block_indices);
  op_matrix_->set_coloring(dim_, block_indices);
  */

  op_matrix_->ComputeInverse();
}


/* ******************************************************************
* Check difference du between the predicted and converged solutions.
* This is a wrapper for various error control methods.
****************************************************************** */
double
MechanicsElasticity_PK::ErrorNorm(Teuchos::RCP<const TreeVector> u,
                                  Teuchos::RCP<const TreeVector> du)
{
  if (!u->Data()->HasComponent("node")) return 0.0;

  const auto& uv = *u->Data()->ViewComponent("node");
  const auto& duv = *du->Data()->ViewComponent("node");

  double error(0.0);

  for (int v = 0; v < nnodes_owned_; v++) {
    for (int k = 0; k < dim_; ++k) {
      double tmp = fabs(duv[k][v]) / (fabs(uv[k][v]) + 1.0);
      error = std::max(error, tmp);
    }
  }

#ifdef HAVE_MPI
  double buf = error;
  du->Comm()->MaxAll(&buf, &error, 1);
#endif

  return error;
}

} // namespace Mechanics
} // namespace Amanzi
