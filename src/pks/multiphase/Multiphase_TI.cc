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

// Amanzi
#include "Tensor.hh"

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
  auto f0 = f->SubVector(0)->Data();

  // -- upwind relative permeability
  S_->GetFieldEvaluator(relperm_liquid_key_)->HasFieldChanged(S_.ptr(), passwd_);
  auto relperm_l = S_->GetFieldData(relperm_liquid_key_, relperm_liquid_key_);

  auto kr = CreateCVforUpwind(mesh_);
  *kr->ViewComponent("cell") = *relperm_l->ViewComponent("cell");
  upwind_->Compute(*kr, *kr, op_bcs_[0]->bc_model(), *kr);

  // -- scale the upwinded relperm assuming constant viscosity
  double factor = eta_l_ / mu_l_;
  kr->Scale(factor);

  // -- create operator
  Teuchos::RCP<std::vector<WhetStone::Tensor> > Kptr = Teuchos::rcpFromRef(K_);

  auto& pde = pde_diff_K_;
  pde->Setup(Kptr, kr, Teuchos::null, rho_l_, gravity_);
  pde->global_operator()->Init();
  pde->ApplyBCs(true, true, true);
  pde->UpdateMatrices(Teuchos::null, Teuchos::null);
  pde->global_operator()->ComputeNegativeResidual(*u0, *f0);
 
  // ---------------------------
  // process chemical components
  // ---------------------------
  const auto& u1 = u_new->SubVector(1)->Data();
  auto f1 = f->SubVector(1)->Data();
  CompositeVector fadd(*f0);

  // -- gas pressure
  S_->GetFieldEvaluator(pressure_gas_key_)->HasFieldChanged(S_.ptr(), passwd_);
  auto pg = S_->GetFieldData(pressure_gas_key_, pressure_gas_key_);

  // -- mass density of gas phase (create evaluator ??? FIXME)
  const auto& eta_g = S_->GetFieldData("molar_density_gas");
  auto rho_g = Teuchos::rcp(new CompositeVector(*eta_g));

  // -- update relative premeability of gas phase
  S_->GetFieldEvaluator(relperm_gas_key_)->HasFieldChanged(S_.ptr(), passwd_);
  auto relperm_g = S_->GetFieldData(relperm_gas_key_, relperm_gas_key_);

  const auto& u1_c = *u1->ViewComponent("cell");
  const auto& eta_gc = *eta_g->ViewComponent("cell");
  const auto& relperm_lc = *relperm_l->ViewComponent("cell");
  const auto& relperm_gc = *relperm_g->ViewComponent("cell");
  auto& kr_c = *kr->ViewComponent("cell");

  for (int i = 0; i < num_primary_; ++i) {
    // -- upwind the scaled relative permeability kr xi eta_l / mu_l
    for (int c = 0; c < ncells_owned_; ++c) {
      kr_c[0][c] = u1_c[i][c] * eta_l_ * relperm_lc[0][c] / mu_l_;
    }
    upwind_->Compute(*kr, *kr, op_bcs_[0]->bc_model(), *kr);

    // -- init with transport in liquid phase
    auto& pdeK = pde_diff_K_;
    pdeK->Setup(Kptr, kr, Teuchos::null, rho_l_, gravity_);
    pdeK->SetBCs(op_bcs_[0], op_bcs_[0]);
    pdeK->global_operator()->Init();
    pdeK->ApplyBCs(true, true, true);
    pdeK->UpdateMatrices(Teuchos::null, Teuchos::null);
    pdeK->global_operator()->ComputeNegativeResidual(*u0, *f1);

    // -- upwind the scaled relative permeability for gas phase
    for (int c = 0; c < ncells_owned_; ++c) {
      kr_c[0][c] = (1.0 - u1_c[i][c]) * eta_gc[0][c] * relperm_gc[0][c] / mu_g_;
    }
    upwind_->Compute(*kr, *kr, op_bcs_[0]->bc_model(), *kr);

    // -- add transport in gas phase
    pdeK->Setup(Kptr, kr, Teuchos::null, rho_g, gravity_);
    pdeK->global_operator()->Init();
    pdeK->ApplyBCs(true, true, true);
    pdeK->UpdateMatrices(Teuchos::null, Teuchos::null);
    pdeK->global_operator()->ComputeNegativeResidual(*pg, fadd);
    f1->Update(1.0, fadd, 1.0);

    // -- add molecular diffusion using harmoic-mean formula
    CompositeVector u1i(fadd);
    *u1i.ViewComponent("cell") = *u1->ViewComponent("cell");

    auto D = Teuchos::rcp(new std::vector<WhetStone::Tensor>);
    WhetStone::Tensor Dc(dim_, 1);
    for (int c = 0; c < ncells_owned_; ++c) {
      Dc(0, 0) = eta_l_ * mol_diff_l_[i] - eta_gc[0][c] * mol_diff_g_[i];
      D->push_back(Dc);
    }

    auto& pdeD = pde_diff_D_;
    pdeD->Setup(D, Teuchos::null, Teuchos::null);
    pdeK->SetBCs(op_bcs_[1], op_bcs_[1]);
    pdeD->global_operator()->Init();
    pdeD->ApplyBCs(true, true, true);
    pdeD->UpdateMatrices(Teuchos::null, Teuchos::null);
    pdeD->global_operator()->ComputeNegativeResidual(u1i, fadd);

    auto f1c = (*f1->ViewComponent("cell"))(i);
    const auto& fadd_c = *fadd.ViewComponent("cell");
    f1c->Update(1.0, fadd_c, 1.0);
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

  // ------------------
  // process constraint
  // ------------------
  f->SubVector(2)->Data()->PutScalar(0.0);
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

