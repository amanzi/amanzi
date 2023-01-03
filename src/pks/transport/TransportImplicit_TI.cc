/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Transport PK

*/

#include "TransportImplicit_PK.hh"

namespace Amanzi {
namespace Transport {

/* *******************************************************************
* Routine takes single component and returns functional value
****************************************************************** */
void
TransportImplicit_PK::FunctionalResidual(double t_old,
                                         double t_new,
                                         Teuchos::RCP<TreeVector> u_old,
                                         Teuchos::RCP<TreeVector> u_new,
                                         Teuchos::RCP<TreeVector> f)
{
  double dtp(t_new - t_old);

  S_->GetEvaluator(wc_key_).Update(*S_, "transport");
  const auto wc_c = S_->Get<CompositeVector>(wc_key_, Tags::DEFAULT).ViewComponent("cell");
  const auto wcprev_c = S_->Get<CompositeVector>(prev_wc_key_, Tags::DEFAULT).ViewComponent("cell");
  const auto& sat_c =
    *S_->Get<CompositeVector>(saturation_liquid_key_, Tags::DEFAULT).ViewComponent("cell");

  wc_start = wcprev_c;
  wc_end = wc_c;

  // since dA(C)/dt = F(C), explicit and implicit schemes use F(C) with signs.
  f->PutScalar(0.0);
  FunctionalTimeDerivative_MUSCL_(t_new, *u_new->Data(), *f->Data(), false);
  f->Scale(-1.0);

  auto& uold_c = *u_old->Data()->ViewComponent("cell");
  auto& unew_c = *u_new->Data()->ViewComponent("cell");
  auto& f_c = *f->Data()->ViewComponent("cell");

  if (use_dispersion_) {
    int phase;
    double md;
    CalculateDispersionTensor_(*transport_phi, *wc_c);
    FindDiffusionValue(component_names_[current_component_], &md, &phase);
    if (md != 0.0) CalculateDiffusionTensor_(md, phase, *transport_phi, sat_c, *wc_c);

    op_diff_->global_operator()->Init();
    op_diff_->UpdateMatrices(Teuchos::null, Teuchos::null);
    op_diff_->ApplyBCs(true, true, true);

    CompositeVector g(*f->Data());
    op_diff_->global_operator()->ComputeNegativeResidual(*u_new->Data(), g);
    f->Data()->Update(1.0, g, 1.0);
  }

  // storage term
  for (int c = 0; c < ncells_owned; ++c) {
    double factor = mesh_->getCellVolume(c) / dtp;
    f_c[0][c] += ((*wc_c)[0][c] * unew_c[0][c] - (*wcprev_c)[0][c] * uold_c[0][c]) * factor;
  }
}


/* ******************************************************************
* We assume that operators were populated during residual calculation.
****************************************************************** */
void
TransportImplicit_PK::UpdatePreconditioner(double tp, Teuchos::RCP<const TreeVector> u, double dtp)
{
  S_->GetEvaluator(wc_key_).Update(*S_, "transport");
  const auto& wc = S_->Get<CompositeVector>(wc_key_, Tags::DEFAULT);

  op_->Init();

  op_acc_->AddAccumulationTerm(wc, dtp, "cell");

  op_adv_->UpdateMatrices(S_->GetPtr<CompositeVector>(vol_flowrate_key_, Tags::DEFAULT).ptr());
  op_adv_->ApplyBCs(true, true, true);

  if (use_dispersion_) {
    op_diff_->UpdateMatrices(Teuchos::null, Teuchos::null);
    op_diff_->ApplyBCs(true, true, true);
  }

  op_->AssembleMatrix();
}


/* ******************************************************************
* Apply preconditioner inv(B) * X.
****************************************************************** */
int
TransportImplicit_PK::ApplyPreconditioner(Teuchos::RCP<const TreeVector> X,
                                          Teuchos::RCP<TreeVector> Y)
{
  return op_pc_solver_->ApplyInverse(*X->Data(), *Y->Data());
}


/* ******************************************************************
* Check difference du between the predicted and converged solutions.
* This is a wrapper for various error control methods.
****************************************************************** */
double
TransportImplicit_PK::ErrorNorm(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<const TreeVector> du)
{
  const Epetra_MultiVector& uc = *u->Data()->ViewComponent("cell");
  const Epetra_MultiVector& duc = *du->Data()->ViewComponent("cell");

  double error(0.0);
  for (int c = 0; c < ncells_owned; c++) {
    double tmp = fabs(duc[0][c]) / (1.0 + fabs(uc[0][c]));
    error = std::max(error, tmp);
  }

  double buf = error;
  du->Comm()->MaxAll(&buf, &error, 1);

  return error;
}

} // namespace Transport
} // namespace Amanzi
