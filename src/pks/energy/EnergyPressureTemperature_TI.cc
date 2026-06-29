/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov
*/

/*
  Energy

*/

// Amanzi
#include "Key.hh"
#include "LScheme_Helpers.hh"
#include "PDE_HelperDiscretization.hh"

// Amanzi::Energy
#include "Evaluator.hh"
#include "EnergyPressureTemperature_PK.hh"

namespace Amanzi {
namespace Energy {

using CV_t = CompositeVector;

/* ******************************************************************
* Computes the non-linear functional g = g(t,u,udot)
****************************************************************** */
void
EnergyPressureTemperature_PK::FunctionalResidual(double t_old,
                                                 double t_new,
                                                 Teuchos::RCP<const TreeVector> u_old,
                                                 Teuchos::RCP<TreeVector> u_new,
                                                 Teuchos::RCP<TreeVector> g)
{
  double dt = t_new - t_old;

  // update BCs, BOUNDARY_FACE component (if any), and conductivity
  temperature_eval_->SetChanged();
  UpdateSourceBoundaryData(t_old, t_new, *u_new->Data());

  S_->GetEvaluator(mol_flowrate_key_).Update(*S_, passwd_);
  auto flux = S_->GetPtr<CompositeVector>(mol_flowrate_key_, Tags::DEFAULT);

  S_->GetEvaluator(conductivity_gen_key_).Update(*S_, passwd_);
  if (upwind_.get()) {
    const auto& conductivity = S_->Get<CompositeVector>(conductivity_gen_key_);
    *upw_conductivity_->ViewComponent("cell") = *conductivity.ViewComponent("cell");

    auto op_bc_temp = S_->GetPtrW<Operators::BCs>(bcs_temperature_key_, Tags::DEFAULT, "state");
    const auto& bc_model = op_bc_temp->bc_model();
    
    Operators::CellToBoundaryFaces(bc_model, *upw_conductivity_);
    upwind_->Compute(*flux, bc_model, *upw_conductivity_);
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

  if (dt > 0) {
    for (int c = 0; c < ncells_owned; ++c) {
      double factor = mesh_->getCellVolume(c) / dt;
      g_c[0][c] += factor * (e1[0][c] - e0[0][c]);
    }
  }

  // advect tmp = molar_density_liquid * enthalpy
  S_->GetEvaluator(enthalpy_key_).Update(*S_, passwd_);
  const auto& enthalpy = S_->Get<CompositeVector>(enthalpy_key_);

  op_advection_->Init();
  op_matrix_advection_->Setup(*flux);
  op_matrix_advection_->UpdateMatrices(flux.ptr());
  op_matrix_advection_->ApplyBCs(false, true, false);

  CompositeVector g_adv(g->Data()->Map());
  op_advection_->ComputeNegativeResidual(enthalpy, g_adv);
  g->Data()->Update(1.0, g_adv, 1.0);
  
    
  //implicit source models
  if (heat_src_) {
    S_->GetEvaluator(heat_src_key_).Update(*S_, passwd_);
    const auto& src_c = *S_->Get<CV_t>(heat_src_key_, Tags::DEFAULT).ViewComponent("cell");

    for (int c = 0; c < ncells_owned; ++c) {
      g_c[0][c] += src_c[0][c] * mesh_->getCellVolume(c);
    }
  }

  //add optional stabilization term
  if (L_scheme_) {
    const auto& stability_c = *S_->Get<CV_t>(L_scheme_stab_key_).ViewComponent("cell");
    const auto& u_prev_c = *S_->Get<CV_t>(L_scheme_prev_key_).ViewComponent("cell");
    const auto& u_new_c = *u_new->Data()->ViewComponent("cell");
    const auto& u_old_c = *u_old->Data()->ViewComponent("cell");

    double delta(0.0), gnorm(0.0), factor, udiff;
    for (int c = 0; c < ncells_owned; ++c) {
      udiff = u_new_c[0][c] - u_prev_c[0][c];
      delta = std::max(delta, std::fabs(udiff));

      factor = mesh_->getCellVolume(c) / dt_;
      gnorm = std::max(gnorm, g_c[0][c] / (factor * e1[0][c])); // true residual

      g_c[0][c] += stability_c[0][c] * udiff * factor;

    }
 
    // save data
    auto& data = S_->GetW<LSchemeData>(L_scheme_data_key_, "state");
    data[temperature_key_].last_step_increment = delta;
    data[temperature_key_].last_step_residual = gnorm;
    data[temperature_key_].ns_itrs[0]++;
  }
}


/* ******************************************************************
* Update the preconditioner at time t and u = up
****************************************************************** */
void
EnergyPressureTemperature_PK::UpdatePreconditioner(double t,
                                                   Teuchos::RCP<const TreeVector> up,
                                                   double dt)
{
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
    *vo_->os() << "updating preconditioner, T=" << t << std::endl;
  }

  // update BCs and conductivity
  UpdateSourceBoundaryData(t, t + dt, *up->Data());

  S_->GetEvaluator(conductivity_gen_key_).Update(*S_, passwd_);

  if (upwind_.get()) {
    const auto& conductivity = S_->Get<CompositeVector>(conductivity_gen_key_);
    *upw_conductivity_->ViewComponent("cell") = *conductivity.ViewComponent("cell");

    auto op_bc_temp = S_->GetPtrW<Operators::BCs>(bcs_temperature_key_, Tags::DEFAULT, "state");
    const auto& bc_model = op_bc_temp->bc_model();

    auto flux = S_->GetPtr<CompositeVector>(mol_flowrate_key_, Tags::DEFAULT);
    Operators::CellToBoundaryFaces(bc_model, *upw_conductivity_);
    upwind_->Compute(*flux, bc_model, *upw_conductivity_);
  }

  // assemble matrices for diffusion operator
  op_preconditioner_->Init();
  op_preconditioner_diff_->UpdateMatrices(Teuchos::null, up->Data().ptr());
  op_preconditioner_diff_->ApplyBCs(true, true, true);

  // update with accumulation terms
  // update the accumulation derivatives, dE/dT
  S_->GetEvaluator(energy_key_).UpdateDerivative(*S_, passwd_, temperature_key_, Tags::DEFAULT);
  auto& dEdT = S_->GetDerivative<CompositeVector>(energy_key_, Tags::DEFAULT, temperature_key_, Tags::DEFAULT);
  auto& dEdT_c = *dEdT.ViewComponent("cell");

  S_->GetEvaluator(ie_rock_key_).UpdateDerivative(*S_, passwd_, temperature_key_, Tags::DEFAULT);
  const auto& dUrdT = S_->GetDerivative<CV_t>(ie_rock_key_, Tags::DEFAULT, temperature_key_, Tags::DEFAULT);
  auto& dUrdT_c = *dUrdT.ViewComponent("cell");

  S_->GetEvaluator(ie_liquid_key_).UpdateDerivative(*S_, passwd_, temperature_key_, Tags::DEFAULT);
  const auto& dUidT = S_->GetDerivative<CV_t>(ie_liquid_key_, Tags::DEFAULT, temperature_key_, Tags::DEFAULT);
  auto& dUidT_c = *dUrdT.ViewComponent("cell");

  S_->GetEvaluator(mol_density_liquid_key_).Update(*S_, passwd_);
  const auto& eta_c = *S_->Get<CV_t>(mol_density_liquid_key_, Tags::DEFAULT).ViewComponent("cell");
  const auto& rho_r = *S_->Get<CV_t>(particle_density_key_, Tags::DEFAULT).ViewComponent("cell");
  const auto& phi_c = *S_->Get<CV_t>(porosity_key_, Tags::DEFAULT).ViewComponent("cell");

  
  if (dt > 0.0) {
    op_acc_->AddAccumulationDelta(*up->Data().ptr(), dEdT, dEdT, dt, "cell");
    for (int c = 0; c < ncells_owned; ++c) {
      if (dEdT_c[0][c] < 0.0) {
        double tmp = phi_c[0][c];
        dEdT_c[0][c] = rho_r[0][c] * dUrdT_c[0][c] * (1.0 - tmp) + eta_c[0][c] * dUidT_c[0][c] * tmp;
      }
    }
  }

  // implicit source models
  if (heat_src_) {
    S_->GetEvaluator(heat_src_key_).UpdateDerivative(*S_, passwd_, temperature_key_, Tags::DEFAULT);
    auto dQdT = S_->GetDerivative<CV_t>(heat_src_key_, Tags::DEFAULT, temperature_key_, Tags::DEFAULT);
    op_acc_->AddAccumulationTerm(dQdT, 1.0, "cell", true);
  }

  // add advection term dHdT
  if (prec_include_enthalpy_) {
    auto flux = S_->GetPtr<CompositeVector>(mol_flowrate_key_, Tags::DEFAULT);

    S_->GetEvaluator(enthalpy_key_).UpdateDerivative(*S_, passwd_, temperature_key_, Tags::DEFAULT);
    auto dHdT = S_->GetDerivativePtr<CompositeVector>(
      enthalpy_key_, Tags::DEFAULT, temperature_key_, Tags::DEFAULT);

    // should we differentiate alpha = (H eta/mu) in div(H q) ?
    op_preconditioner_advection_->Setup(*flux);
    op_preconditioner_advection_->UpdateMatrices(flux.ptr(), dHdT.ptr());
    op_preconditioner_advection_->ApplyBCs(false, true, false);
  }

  if (L_scheme_) {
    const auto& stability = S_->Get<CV_t>(L_scheme_stab_key_);
    op_acc_->AddAccumulationTerm(stability, dt, "cell");
  }

  // verify and finalize preconditioner
  // op_preconditioner_->Verify();
  op_preconditioner_->ComputeInverse();
}


/* ******************************************************************
* TBW
****************************************************************** */
double
EnergyPressureTemperature_PK::ErrorNorm(Teuchos::RCP<const TreeVector> u,
                                        Teuchos::RCP<const TreeVector> du)
{
  Teuchos::OSTab tab = vo_->getOSTab();

  // Relative error in cell-centered temperature
  double error(0.0), error_t(0.0);
  double ref_temp(273.0);

  for (auto comp = u->Data()->begin(); comp != u->Data()->end(); ++comp) {
    const auto& uc = *u->Data()->ViewComponent(*comp);
    const auto& duc = *du->Data()->ViewComponent(*comp);

    int ncomp = uc.MyLength();
    for (int i = 0; i < ncomp; ++i) {
      double tmp = fabs(duc[0][i]) / (fabs(uc[0][i] - ref_temp) + ref_temp);
      if (tmp > error_t) error_t = tmp;
    }
  }

  // Cell error is based upon error in energy conservation relative to
  // a characteristic energy
  double error_e(0.0);
  error = std::max(error_t, error_e);

#ifdef HAVE_MPI
  double buf = error;
  du->Data()->Comm()->MaxAll(&buf, &error, 1); // find the global maximum
#endif

  return error;
}

} // namespace Energy
} // namespace Amanzi
