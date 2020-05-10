/*
  Energy

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov
*/

#include "FieldEvaluator.hh"

#include "EnergyOnePhase_PK.hh"

namespace Amanzi {
namespace Energy {

/* ******************************************************************
* Computes the non-linear functional g = g(t,u,udot)
****************************************************************** */
void EnergyOnePhase_PK::FunctionalResidual(
    double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
    Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> g)
{
  double dt = t_new - t_old;

  // update BCs and conductivity
  temperature_eval_->SetFieldAsChanged(S_.ptr());
  UpdateSourceBoundaryData(t_old, t_new, *u_new->Data());
  S_->GetFieldEvaluator(conductivity_key_)->HasFieldChanged(S_.ptr(), passwd_);

  // assemble residual for diffusion operator
  op_matrix_->Init();
  op_matrix_diff_->UpdateMatrices(Teuchos::null, solution.ptr());
  op_matrix_diff_->ApplyBCs(true, true, true);

  op_matrix_->ComputeNegativeResidual(*u_new->Data(), *g->Data());

  // add accumulation term
  // -- update the energy at the new time.
  S_->GetFieldEvaluator(energy_key_)->HasFieldChanged(S_.ptr(), passwd_);

  const Epetra_MultiVector& e1 = *S_->GetFieldData(energy_key_)->ViewComponent("cell");
  const Epetra_MultiVector& e0 = *S_->GetFieldData(prev_energy_key_)->ViewComponent("cell");
  Epetra_MultiVector& g_c = *g->Data()->ViewComponent("cell");

  for (int c = 0; c < ncells_owned; ++c) {
    double factor = mesh_->cell_volume(c) / dt;
    g_c[0][c] += factor * (e1[0][c] - e0[0][c]);
  }

  // advect tmp = molar_density_liquid * enthalpy 
  S_->GetFieldEvaluator(enthalpy_key_)->HasFieldChanged(S_.ptr(), passwd_);
  const CompositeVector& enthalpy = *S_->GetFieldData(enthalpy_key_);
  const CompositeVector& n_l = *S_->GetFieldData("molar_density_liquid");

  Teuchos::RCP<const CompositeVector> flux = S_->GetFieldData("darcy_flux");
  op_advection_->Init();
  op_matrix_advection_->Setup(*flux);
  op_matrix_advection_->UpdateMatrices(flux.ptr());
  op_matrix_advection_->ApplyBCs(false, true, false);

  CompositeVector tmp(enthalpy);
  tmp.Multiply(1.0, tmp, n_l, 0.0);

  CompositeVector g_adv(g->Data()->Map());
  // op_advection_->Apply(tmp, g_adv);
  op_advection_->ComputeNegativeResidual(tmp, g_adv);
  g->Data()->Update(1.0, g_adv, 1.0);
}


/* ******************************************************************
* Update the preconditioner at time t and u = up
****************************************************************** */
void EnergyOnePhase_PK::UpdatePreconditioner(
    double t, Teuchos::RCP<const TreeVector> up, double dt)
{
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
    *vo_->os() << "updating preconditioner, T=" << t << std::endl;
  }

  // update BCs and conductivity
  UpdateSourceBoundaryData(t, t + dt, *up->Data());
  S_->GetFieldEvaluator(conductivity_key_)->HasFieldChanged(S_.ptr(), passwd_);

  // assemble matrices for diffusion operator
  op_preconditioner_->Init();
  op_preconditioner_diff_->UpdateMatrices(Teuchos::null, up->Data().ptr());
  op_preconditioner_diff_->ApplyBCs(true, true, true);

  // update with accumulation terms
  // update the accumulation derivatives, dE/dT
  S_->GetFieldEvaluator(energy_key_)->HasFieldDerivativeChanged(S_.ptr(), passwd_, "temperature");
  CompositeVector& dEdT = *S_->GetFieldData("denergy_dtemperature", energy_key_);

  if (dt > 0.0) {
    op_acc_->AddAccumulationDelta(*up->Data().ptr(), dEdT, dEdT, dt, "cell");
  }

  // add advection term dHdT
  if (prec_include_enthalpy_) {
    Teuchos::RCP<const CompositeVector> darcy_flux = S_->GetFieldData("darcy_flux");

    S_->GetFieldEvaluator(enthalpy_key_)->HasFieldDerivativeChanged(S_.ptr(), passwd_, "temperature");
    Teuchos::RCP<const CompositeVector> dHdT = S_->GetFieldData("denthalpy_dtemperature");

    op_preconditioner_advection_->Setup(*darcy_flux);
    op_preconditioner_advection_->UpdateMatrices(darcy_flux.ptr(), dHdT.ptr());
    op_preconditioner_advection_->ApplyBCs(false, true, false);
  }

  // finalize preconditioner
  op_preconditioner_->AssembleMatrix();
  op_preconditioner_->UpdatePreconditioner();
}


/* ******************************************************************
* TBW
****************************************************************** */
double EnergyOnePhase_PK::ErrorNorm(Teuchos::RCP<const TreeVector> u,
                                    Teuchos::RCP<const TreeVector> du)
{
  Teuchos::OSTab tab = vo_->getOSTab();

  // Relative error in cell-centered temperature
  const Epetra_MultiVector& uc = *u->Data()->ViewComponent("cell", false);
  const Epetra_MultiVector& duc = *du->Data()->ViewComponent("cell", false);

  int cell_bad;
  double error_t(0.0);
  double ref_temp(273.0);
  for (int c = 0; c < ncells_owned; c++) {
    double tmp = fabs(duc[0][c]) / (fabs(uc[0][c] - ref_temp) + ref_temp);
    if (tmp > error_t) {
      error_t = tmp;
      cell_bad = c;
    } 
  }

  // Cell error is based upon error in energy conservation relative to
  // a characteristic energy
  double error_e(0.0);

  double error = std::max(error_t, error_e);

#ifdef HAVE_MPI
  double buf = error;
  du->Data()->Comm()->MaxAll(&buf, &error, 1);  // find the global maximum
#endif

  return error;
}

}  // namespace Energy
}  // namespace Amanzi

