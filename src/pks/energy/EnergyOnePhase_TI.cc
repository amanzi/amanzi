/*
  Energy

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov
*/

// Amanzi
#include "Key.hh"
#include "PDE_HelperDiscretization.hh"

// Amanzi::Energy
#include "Evaluator.hh"
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
  temperature_eval_->SetChanged();
  UpdateSourceBoundaryData(t_old, t_new, *u_new->Data());

  auto flux = S_->GetPtr<CompositeVector>(darcy_flux_key_, Tags::DEFAULT);

  S_->GetEvaluator(conductivity_key_).Update(*S_, passwd_);
  if (upwind_.get()) {
     const auto& conductivity = S_->Get<CompositeVector>(conductivity_key_);
    *upw_conductivity_->ViewComponent("cell") = *conductivity.ViewComponent("cell");

    const auto& bc_model = op_bc_->bc_model();
    Operators::CellToBoundaryFaces(bc_model, *upw_conductivity_);
    upwind_->Compute(*flux, *u_new->Data(), bc_model, *upw_conductivity_);
  }

  // assemble residual for diffusion operator
  op_matrix_->Init();
  op_matrix_diff_->UpdateMatrices(Teuchos::null, solution.ptr());
  op_matrix_diff_->ApplyBCs(true, true, true);

  // add sources
  CompositeVector& rhs = *op_matrix_->rhs();
  AddSourceTerms(rhs);
 
  op_matrix_->ComputeNegativeResidual(*u_new->Data(), *g->Data());

  // add accumulation term
  // -- update the energy at the new time.
  S_->GetEvaluator(energy_key_).Update(*S_, passwd_);

  const auto& e1 = *S_->Get<CompositeVector>(energy_key_).ViewComponent("cell");
  const auto& e0 = *S_->Get<CompositeVector>(prev_energy_key_).ViewComponent("cell");
  Epetra_MultiVector& g_c = *g->Data()->ViewComponent("cell");

  for (int c = 0; c < ncells_owned; ++c) {
    double factor = mesh_->cell_volume(c) / dt;
    g_c[0][c] += factor * (e1[0][c] - e0[0][c]);
  }

  // advect tmp = molar_density_liquid * enthalpy 
  S_->GetEvaluator(enthalpy_key_).Update(*S_, passwd_);
  const auto& enthalpy = S_->Get<CompositeVector>(enthalpy_key_);
  const auto& n_l = S_->Get<CompositeVector>(mol_density_liquid_key_);

  op_advection_->Init();
  op_matrix_advection_->Setup(*flux);
  op_matrix_advection_->UpdateMatrices(flux.ptr());
  op_matrix_advection_->ApplyBCs(false, true, false);

  CompositeVector tmp(enthalpy);
  tmp.Multiply(1.0, tmp, n_l, 0.0);

  CompositeVector g_adv(g->Data()->Map());
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

  S_->GetEvaluator(conductivity_key_).Update(*S_, passwd_);

  if (upwind_.get()) {
    const auto& conductivity = S_->Get<CompositeVector>(conductivity_key_);
    *upw_conductivity_->ViewComponent("cell") = *conductivity.ViewComponent("cell");

    auto flux = S_->GetPtr<CompositeVector>(darcy_flux_key_, Tags::DEFAULT);
    const auto& bc_model = op_bc_->bc_model();
    Operators::CellToBoundaryFaces(bc_model, *upw_conductivity_);
    upwind_->Compute(*flux, *up->Data(), bc_model, *upw_conductivity_);
  }

  // assemble matrices for diffusion operator
  op_preconditioner_->Init();
  op_preconditioner_diff_->UpdateMatrices(Teuchos::null, up->Data().ptr());
  op_preconditioner_diff_->ApplyBCs(true, true, true);

  // update with accumulation terms
  // update the accumulation derivatives, dE/dT
  S_->GetEvaluator(energy_key_).UpdateDerivative(*S_, passwd_, temperature_key_, Tags::DEFAULT);
  auto& dEdT = S_->GetDerivativeW<CompositeVector>(energy_key_, Tags::DEFAULT, temperature_key_, Tags::DEFAULT, energy_key_);

  if (dt > 0.0) {
    op_acc_->AddAccumulationDelta(*up->Data().ptr(), dEdT, dEdT, dt, "cell");
  }

  // add advection term dHdT
  if (prec_include_enthalpy_) {
    auto flux = S_->GetPtr<CompositeVector>(darcy_flux_key_, Tags::DEFAULT);

    S_->GetEvaluator(enthalpy_key_).UpdateDerivative(*S_, passwd_, temperature_key_, Tags::DEFAULT);
    auto dHdT = Teuchos::rcp(new CompositeVector(S_->GetDerivative<CompositeVector>(
        enthalpy_key_, Tags::DEFAULT, temperature_key_, Tags::DEFAULT)));

    const auto& n_l = S_->Get<CompositeVector>(mol_density_liquid_key_);
    dHdT->Multiply(1.0, *dHdT, n_l, 0.0);

    if (S_->HasEvaluator(mol_density_liquid_key_)) {
      auto& eval = S_->GetEvaluator(mol_density_liquid_key_);
      if (eval.get_type() == EvaluatorType::SECONDARY) {
        eval.UpdateDerivative(*S_, passwd_, temperature_key_, Tags::DEFAULT);
        auto dRdT = Teuchos::rcp(new CompositeVector(S_->GetDerivative<CompositeVector>(
            mol_density_liquid_key_, Tags::DEFAULT, temperature_key_, Tags::DEFAULT)));

        const auto& H_l = S_->Get<CompositeVector>(mol_density_liquid_key_);
        dRdT->Multiply(1.0, *dRdT, H_l, 0.0);
        dHdT->Update(1.0, *dRdT, 1.0);
      }
    }

    op_preconditioner_advection_->Setup(*flux);
    op_preconditioner_advection_->UpdateMatrices(flux.ptr(), dHdT.ptr());
    op_preconditioner_advection_->ApplyBCs(false, true, false);
  }

  // verify and finalize preconditioner
  op_preconditioner_->Verify();
  op_preconditioner_->ComputeInverse();
}


/* ******************************************************************
* TBW
****************************************************************** */
double EnergyOnePhase_PK::ErrorNorm(Teuchos::RCP<const TreeVector> u,
                                    Teuchos::RCP<const TreeVector> du)
{
  Teuchos::OSTab tab = vo_->getOSTab();

  // Relative error in cell-centered temperature
  const Epetra_MultiVector& uc = *u->Data()->ViewComponent("cell");
  const Epetra_MultiVector& duc = *du->Data()->ViewComponent("cell");

  double error(0.0), error_t(0.0);
  double ref_temp(273.0);
  for (int c = 0; c < ncells_owned; c++) {
    double tmp = fabs(duc[0][c]) / (fabs(uc[0][c] - ref_temp) + ref_temp);
    if (tmp > error_t) {
      error_t = tmp;
    } 
  }

  // Cell error is based upon error in energy conservation relative to
  // a characteristic energy
  double error_e(0.0);
  error = std::max(error_t, error_e);

#ifdef HAVE_MPI
  double buf = error;
  du->Data()->Comm()->MaxAll(&buf, &error, 1);  // find the global maximum
#endif

  return error;
}

}  // namespace Energy
}  // namespace Amanzi

