/*
  Flow PK
 
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/


#include "Darcy_PK.hh"

namespace Amanzi {
namespace Flow {

/* ******************************************************************
* Calculate f(u, du/dt) = A*u - rhs.
****************************************************************** */
void Darcy_PK::FunctionalResidual(
    double t_old, double t_new, 
    Teuchos::RCP<TreeVector> u_old, Teuchos::RCP<TreeVector> u_new, 
    Teuchos::RCP<TreeVector> f)
{ 
  dt_ = t_new - t_old;

  // refresh data
  UpdateSourceBoundaryData(t_old, t_new, *u_new->Data());

  // calculate and assemble elemental stiffness matrices
  double factor = 1.0 / g_;
  const auto& ss = S_->Get<CompositeVector>(specific_storage_key_);
  CompositeVector ss_g(ss); 
  ss_g.Update(0.0, ss, factor);

  factor = 1.0 / (g_ * dt_);
  CompositeVector sy_g(*specific_yield_copy_); 
  sy_g.Scale(factor);

  if (flow_on_manifold_) {
    const auto& aperture = S_->Get<CompositeVector>(aperture_key_, Tags::DEFAULT);
    ss_g.Multiply(1.0, ss_g, aperture, 0.0);
    sy_g.Multiply(1.0, sy_g, aperture, 0.0);
  }

  op_->RestoreCheckPoint();
  op_acc_->AddAccumulationDelta(*u_old->Data(), ss_g, ss_g, dt_, "cell");
  op_acc_->AddAccumulationDeltaNoVolume(*u_old->Data(), sy_g, "cell");

  // Peaceman model
  if (S_->HasRecord("well_index")) {
    const auto& wi = S_->Get<CompositeVector>("well_index");
    op_acc_->AddAccumulationTerm(wi, "cell");
  }

  // optional update of diffusion operator
  if (flow_on_manifold_) {
    S_->GetEvaluator(permeability_eff_key_).Update(*S_, "flow");
    op_diff_->UpdateMatrices(Teuchos::null, Teuchos::null);
  }
  op_diff_->ApplyBCs(true, true, true);

  CompositeVector& rhs = *op_->rhs();
  AddSourceTerms(rhs);

  // if the matrix was assemble, it will be used in Apply. Due to new
  // accumulation term, we either have to destroy or reassemble it 
  op_->ComputeNegativeResidual(*u_new->Data(), *f->Data());
}


/* ******************************************************************
* Apply preconditioner inv(B) * X.                                                 
****************************************************************** */
int Darcy_PK::ApplyPreconditioner(Teuchos::RCP<const TreeVector> X, 
                                  Teuchos::RCP<TreeVector> Y)
{
  num_itrs_++;
  Y->PutScalar(0.0);

  return op_->ApplyInverse(*X->Data(), *Y->Data());
}


/* ******************************************************************
* Update preconditioner for the interval (tp-dtp, tp].
****************************************************************** */
void Darcy_PK::UpdatePreconditioner(double tp, Teuchos::RCP<const TreeVector> u, double dtp)
{
  op_->ComputeInverse();
}


/* ******************************************************************
* Check difference du between the predicted and converged solutions.
****************************************************************** */
double Darcy_PK::ErrorNorm(Teuchos::RCP<const TreeVector> u, 
                           Teuchos::RCP<const TreeVector> du)
{
  const Epetra_MultiVector& uc = *u->Data()->ViewComponent("cell");
  const Epetra_MultiVector& duc = *du->Data()->ViewComponent("cell");

  double error(0.0);

  for (int c = 0; c < ncells_owned; c++) {
    double tmp = fabs(duc[0][c]) / (fabs(uc[0][c] - atm_pressure_) + atm_pressure_);
    error = std::max(error, tmp);
  }

  /*
  const Epetra_MultiVector& uf = *u->Data()->ViewComponent("face");
  const Epetra_MultiVector& duf = *du->Data()->ViewComponent("face");

  for (int f = 0; f < nfaces_owned; f++) {
    double tmp = fabs(duf[0][f]) / (fabs(uf[0][f] - atm_pressure_) + atm_pressure_);
    error = std::max(error, tmp);
  }
  */

#ifdef HAVE_MPI
    double tmp = error;
    u->Comm()->MaxAll(&tmp, &error, 1);
#endif

  return error;
}


/* ******************************************************************
* Estimate dT increase factor by comparing the 1st and 2nd order
* time approximations.
****************************************************************** */
double Darcy_PK::ErrorEstimate_(double* dt_factor)
{
  double tol, atol(1.0), rtol(1e-5), error, error_max(0.0), p(101325.0);
  double dt_factor_cell;

  *dt_factor = 100.0;
  for (int c = 0; c < ncells_owned; c++) {
    error = fabs((*pdot_cells)[c] - (*pdot_cells_prev)[c]) * dt_ / 2;
    // tol = rtol * fabs(p_cell[0][c]) + atol;
    tol = rtol * p + atol;

    dt_factor_cell = sqrt(tol / std::max(error, FLOW_DT_ADAPTIVE_ERROR_TOLERANCE));
    *dt_factor = std::min(*dt_factor, dt_factor_cell);

    error_max = std::max(error_max, error - tol);
  }

  *dt_factor *= FLOW_DT_ADAPTIVE_SAFETY_FACTOR;
  *dt_factor = std::min(*dt_factor, FLOW_DT_ADAPTIVE_INCREASE);
  *dt_factor = std::max(*dt_factor, FLOW_DT_ADAPTIVE_REDUCTION);

#ifdef HAVE_MPI
    double dt_tmp = *dt_factor;
    solution->Comm()->MinAll(&dt_tmp, dt_factor, 1);  // find the global minimum
 
    double error_tmp = error_max;
    solution->Comm()->MaxAll(&error_tmp, &error_max, 1);  // find the global maximum
#endif

  return error_max;
}

}  // namespace Flow
}  // namespace Amanzi


