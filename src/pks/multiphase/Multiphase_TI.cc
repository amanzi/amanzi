/*
  MultiPhase PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Quan Bui (mquanbui@math.umd.edu)
           Konstantin Lipnikov (lipnikov@lanl.gov)

  Multiphase multi-component flow, see Multiphase_PK.cc for more detail.
*/


// TPLs
#include "Teuchos_RCP.hpp"

// Multiphase
#include "Multiphase_PK.hh"
#include "TotalComponentStorage.hh"

namespace Amanzi {
namespace Multiphase {

/* ******************************************************************
* F(u) implementation
****************************************************************** */
void Multiphase_PK::FunctionalResidual(double t_old, double t_new, 
                                       Teuchos::RCP<TreeVector> u_old,
                                       Teuchos::RCP<TreeVector> u_new,
                                       Teuchos::RCP<TreeVector> f) 
{
  double dtp = t_new - t_old;

  // trigger update of primary variables
  saturation_liquid_eval_->SetFieldAsChanged(S_.ptr());

  // -----------------------
  // process water component
  // -----------------------
  // -- init with diffusion operator for water
  const auto& u0 = u_new->SubVector(0)->Data();
  CompositeVector f1(*u0);

  // -- upwind relative permeability
  S_->GetFieldEvaluator(relperm_liquid_key_)->HasFieldChanged(S_.ptr(), passwd_);
  auto krl = S_->GetFieldData(relperm_liquid_key_, relperm_liquid_key_);

  upwind_->Compute(*krl, *krl, op_bcs_[0]->bc_model(), *krl);
  krl->Scale(rho_l_ / mu_l_);

  // -- create operator
  Teuchos::RCP<std::vector<WhetStone::Tensor> > Kptr = Teuchos::rcpFromRef(K_);
  const auto& rho = S_->GetFieldData("molar_density_liquid");

  auto pde = pde_matrix_diff_[0];
  pde->Setup(Kptr, krl, Teuchos::null, rho, gravity_);
  pde->global_operator()->Init();
  pde->ApplyBCs(true, true, true);
  pde->UpdateMatrices(Teuchos::null, Teuchos::null);
  pde->global_operator()->ComputeNegativeResidual(*u0, f1);
  u0->Update(1.0, f1, 0.0);
 
  // ---------------------------
  // process chemical components
  // ---------------------------
  const auto& u1 = u_new->SubVector(1)->Data();

  // -- gas pressure
  S_->GetFieldEvaluator(pressure_gas_key_)->HasFieldChanged(S_.ptr(), passwd_);
  auto pg = S_->GetFieldData(pressure_gas_key_, pressure_gas_key_);

  // -- mass density of gas phase

  // -- upwind relative permeability
  S_->GetFieldEvaluator(relperm_gas_key_)->HasFieldChanged(S_.ptr(), passwd_);
  auto krg = S_->GetFieldData(relperm_gas_key_, relperm_gas_key_);

  upwind_->Compute(*krg, *krg, op_bcs_[0]->bc_model(), *krg);
  krg->Scale(1.0 / mu_g_);

  for (int i = 0; i < num_primary_; ++i) {
    // -- init with advection operator for specie
  }

  // add accumulation term
  for (int i = 0; i < num_primary_ + 1; ++i) {
    std::string name = eval_acc_[i];
    std::string prev_name = "prev_" + name;
    if (name == "") continue;

    auto eval = S_->GetFieldEvaluator(name);
    if (i > 0) Teuchos::rcp_dynamic_cast<TotalComponentStorage>(eval)->set_component_id(i - 1);

    eval->HasFieldChanged(S_.ptr(), "multiphase");
    const auto& total_c = *S_->GetFieldData(name)->ViewComponent("cell");
    const auto& total_prev_c = *S_->GetFieldData(prev_name)->ViewComponent("cell");

    auto& f_c = *f->SubVector(i)->Data()->ViewComponent("cell");

    for (int c = 0; c < ncells_owned_; ++c) {
      double factor = mesh_->cell_volume(c) / dtp;
      f_c[0][c] += (total_c[0][c] - total_prev_c[0][c]) * factor;
    }
  }
}


/* ******************************************************************
* Update new preconditioner on the interval (tp-dtp, tp].
****************************************************************** */
void Multiphase_PK::UpdatePreconditioner(double tp, Teuchos::RCP<const TreeVector> u, double dtp)
{
  double t_old = tp - dtp;
}


/* ******************************************************************
* Apply preconditioner inv(B) * X.                                                 
****************************************************************** */
int Multiphase_PK::ApplyPreconditioner(Teuchos::RCP<const TreeVector> X, 
                                       Teuchos::RCP<TreeVector> Y)
{
  Y->PutScalar(0.0);
  // return op_pc_solver_->ApplyInverse(*X->Data(), *Y->Data());
  *Y = *X;
  return 0;
}


/* ******************************************************************
* Monitor l2 norm of residual
****************************************************************** */
double Multiphase_PK::ErrorNorm(Teuchos::RCP<const TreeVector> u,
                                Teuchos::RCP<const TreeVector> du) 
{
  double du_l2 = 0.0;
  double resnorm_p, resnorm_s, resnorm_r;

  du->SubVector(0)->Data()->Norm2(&resnorm_p);
  du->SubVector(1)->Data()->Norm2(&resnorm_s);
  du->SubVector(2)->Data()->Norm2(&resnorm_r);
  printf("resnorm_p = %4.6e, resnorm_s = %4.6e, resnorm_r = %4.6e \n", resnorm_p, resnorm_s, resnorm_r);

  du->Norm2(&du_l2);
  return du_l2;
}


/********************************************************************
* Modifies nonlinear update du using .. TBW
****************************************************************** */
AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
Multiphase_PK::ModifyCorrection(double h, Teuchos::RCP<const TreeVector> res,
                               Teuchos::RCP<const TreeVector> u,
                               Teuchos::RCP<TreeVector> du)
{
  Teuchos::RCP<CompositeVector> rho_next = Teuchos::rcp(new CompositeVector(*u->SubVector(2)->Data()));
  rho_next->Update(-1.0, *du->SubVector(2)->Data(), 1.0);

  Epetra_MultiVector& rho_c = *rho_next->ViewComponent("cell");
  for (int c = 0; c < rho_c.MyLength(); c++) {
    rho_c[0][c] = std::max(0.0, rho_c[0][c]);
  }

  du->SubVector(2)->Data()->Update(1.0, *u->SubVector(2)->Data(), -1.0, *rho_next, 0.0);

  return AmanziSolvers::FnBaseDefs::CORRECTION_MODIFIED;
}

}  // namespace Multiphase
}  // namespace Amanzi

