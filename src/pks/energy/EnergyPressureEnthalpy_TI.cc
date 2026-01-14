/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov
*/

/*
  Energy PK
*/

// Amanzi
#include "Key.hh"
#include "LScheme_Helpers.hh"
#include "PDE_HelperDiscretization.hh"

// Amanzi::Energy
#include "Evaluator.hh"
#include "EnergyPressureEnthalpy_PK.hh"

namespace Amanzi {
namespace Energy {

using CV_t = CompositeVector;

/* ******************************************************************
* Computes the non-linear functional g = g(t,u,udot)
****************************************************************** */
void
EnergyPressureEnthalpy_PK::FunctionalResidual(double t_old,
                                              double t_new,
                                              Teuchos::RCP<const TreeVector> u_old,
                                              Teuchos::RCP<TreeVector> u_new,
                                              Teuchos::RCP<TreeVector> g)
{
  double dt = t_new - t_old;

  // update BCs and  BOUNDARY_FACE component, if any
  enthalpy_eval_->SetChanged();
  UpdateSourceBoundaryData(t_old, t_new, *u_new->Data());

  // we use the chain rule for gradient of T to avoid unphysical impact
  // of gradient of h to the gradient of T in region 4.
  // add residual for enthalpy diffusion operator 

  S_->GetEvaluator(conductivity_key_).Update(*S_, passwd_);
  const auto& conductivity = S_->Get<CV_t>(conductivity_key_, Tags::DEFAULT);
  auto coef = Teuchos::rcp(new CompositeVector(conductivity));

  S_->GetEvaluator(temperature_key_).UpdateDerivative(*S_, passwd_, enthalpy_key_, Tags::DEFAULT);
  const auto& dTdh = S_->GetDerivative<CV_t>(temperature_key_, Tags::DEFAULT, enthalpy_key_, Tags::DEFAULT);
  coef->Multiply(1.0, *coef, dTdh, 0.0);

  op_matrix_->Init();
  op_matrix_diff_->SetScalarCoefficient(coef, Teuchos::null);
  op_matrix_diff_->UpdateMatrices(Teuchos::null, solution.ptr());
  op_matrix_diff_->ApplyBCs(true, true, true);

  // add sources
  CompositeVector& rhs = *op_matrix_->rhs();
  AddSourceTerms(rhs);

  op_matrix_->ComputeNegativeResidual(*u_new->Data(), *g->Data());

  // add diffusion due to pressure
  CompositeVector g_adv(g->Data()->Map());
  const auto& pres = S_->Get<CV_t>(pressure_key_, Tags::DEFAULT);

  S_->GetEvaluator(temperature_key_).UpdateDerivative(*S_, passwd_, pressure_key_, Tags::DEFAULT);
  const auto& dTdp = S_->GetDerivative<CV_t>(temperature_key_, Tags::DEFAULT, pressure_key_, Tags::DEFAULT);
  coef->Multiply(1.0, conductivity, dTdp, 0.0);

  op_matrix_diff_pres_->global_operator()->Init();
  op_matrix_diff_pres_->SetScalarCoefficient(coef, Teuchos::null);
  op_matrix_diff_pres_->UpdateMatrices(Teuchos::null, Teuchos::null);

  op_matrix_diff_pres_->global_operator()->ComputeNegativeResidual(pres, g_adv);
  g->Data()->Update(1.0, g_adv, 1.0);

  // add accumulation term
  S_->GetEvaluator(energy_key_).Update(*S_, passwd_);

  const auto& e1 = *S_->Get<CompositeVector>(energy_key_).ViewComponent("cell");
  const auto& e0 = *S_->Get<CompositeVector>(prev_energy_key_).ViewComponent("cell");
  Epetra_MultiVector& g_c = *g->Data()->ViewComponent("cell");

  for (int c = 0; c < ncells_owned; ++c) {
    double factor = mesh_->getCellVolume(c) / dt;
    g_c[0][c] += factor * (e1[0][c] - e0[0][c]);
  }

  // add enthalpy advection
  S_->GetEvaluator(mol_flowrate_key_).Update(*S_, passwd_);
  auto flux = S_->GetPtr<CompositeVector>(mol_flowrate_key_, Tags::DEFAULT);

  op_advection_->Init();
  op_matrix_advection_->Setup(*flux);
  op_matrix_advection_->UpdateMatrices(flux.ptr());
  op_matrix_advection_->ApplyBCs(false, true, false);

  op_advection_->ComputeNegativeResidual(*u_new->Data(), g_adv);
  g->Data()->Update(1.0, g_adv, 1.0);
}


/* ******************************************************************
* TBW
****************************************************************** */
double
EnergyPressureEnthalpy_PK::ErrorNorm(Teuchos::RCP<const TreeVector> u,
                                     Teuchos::RCP<const TreeVector> du)
{
  double error(0.0), error_h(0.0);
  double ref_enthalpy(18000.0);

  for (auto comp = u->Data()->begin(); comp != u->Data()->end(); ++comp) {
    const auto& uc = *u->Data()->ViewComponent(*comp);
    const auto& duc = *du->Data()->ViewComponent(*comp);

    int ncomp = uc.MyLength();
    for (int i = 0; i < ncomp; ++i) {
      double tmp = fabs(duc[0][i]) / (fabs(uc[0][i]) + ref_enthalpy);
      error_h = std::max(error_h, tmp);
    }
  }

  // Cell error is based upon error in energy conservation relative to
  // a characteristic energy
  double error_e(0.0);
  error = std::max(error_h, error_e);

  double buf = error;
  du->Data()->Comm()->MaxAll(&buf, &error, 1);

  return error;
}


/* ******************************************************************
* Update the preconditioner at time t and u = up
****************************************************************** */
void
EnergyPressureEnthalpy_PK::UpdatePreconditioner(double t,
                                                Teuchos::RCP<const TreeVector> up,
                                                double dt)
{
  // update BCs
  UpdateSourceBoundaryData(t, t + dt, *up->Data());

  // assemble matrices for diffusion operator
  S_->GetEvaluator(conductivity_key_).Update(*S_, passwd_);
  const auto& conductivity = S_->Get<CV_t>(conductivity_key_, Tags::DEFAULT);
  auto coef = Teuchos::rcp(new CompositeVector(conductivity));

  S_->GetEvaluator(temperature_key_).UpdateDerivative(*S_, passwd_, enthalpy_key_, Tags::DEFAULT);
  const auto& dTdh = S_->GetDerivative<CV_t>(temperature_key_, Tags::DEFAULT, enthalpy_key_, Tags::DEFAULT);
  coef->Multiply(1.0, *coef, dTdh, 0.0);

  op_preconditioner_->Init();
  op_preconditioner_diff_->SetScalarCoefficient(coef, Teuchos::null);
  op_preconditioner_diff_->UpdateMatrices(Teuchos::null, up->Data().ptr());
  op_preconditioner_diff_->ApplyBCs(true, true, true);

  // add the accumulation term derivatives, dE/dT
  S_->GetEvaluator(energy_key_).UpdateDerivative(*S_, passwd_, enthalpy_key_, Tags::DEFAULT);
  const auto& dEdh = S_->GetDerivative<CompositeVector>(energy_key_, Tags::DEFAULT, enthalpy_key_, Tags::DEFAULT);
  op_acc_->AddAccumulationTerm(dEdh, dt, "cell");

  // add matrices for advection term 
  auto flux = S_->GetPtr<CompositeVector>(mol_flowrate_key_, Tags::DEFAULT);

  op_preconditioner_advection_->Setup(*flux);
  op_preconditioner_advection_->UpdateMatrices(flux.ptr());
  op_preconditioner_advection_->ApplyBCs(false, true, false);

  // verify and finalize preconditioner
  // op_preconditioner_->Verify();
  op_preconditioner_->ComputeInverse();
}

} // namespace Energy
} // namespace Amanzi

