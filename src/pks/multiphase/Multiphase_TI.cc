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
#include "PDE_Accumulation.hh"
#include "PDE_AdvectionUpwind.hh"
#include "PDE_DiffusionFV.hh"
#include "PDE_DiffusionFVwithGravity.hh"
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

  PopulateBCs(0);

  // trigger update of primary variables
  saturation_liquid_eval_->SetFieldAsChanged(S_.ptr());

  // -----------------------
  // process water component
  // -----------------------
  const auto& u0 = u_new->SubVector(0)->Data();
  auto f0 = f->SubVector(0)->Data();

  // -- upwind relative permeability
  S_->GetFieldEvaluator(relperm_liquid_key_)->HasFieldChanged(S_.ptr(), passwd_);
  auto relperm_l = S_->GetFieldData(relperm_liquid_key_, relperm_liquid_key_);

  auto kr = CreateCVforUpwind(mesh_);
  *kr->ViewComponent("cell") = *relperm_l->ViewComponent("cell");
  kr->ViewComponent("dirichlet_faces")->PutScalar(1.0);  // FIXME
  upwind_->Compute(*kr, *kr, op_bcs_[0]->bc_model(), *kr);

  // -- scale the upwinded relperm assuming constant viscosity
  double factor = eta_l_ / mu_l_;
  kr->Scale(factor);

  // -- create operator
  Teuchos::RCP<std::vector<WhetStone::Tensor> > Kptr = Teuchos::rcpFromRef(K_);

  auto& pde = pde_diff_K_;
  pde->Setup(Kptr, kr, Teuchos::null, rho_l_, gravity_);
  pde->global_operator()->Init();
  pde->UpdateMatrices(Teuchos::null, Teuchos::null);
  pde->ApplyBCs(true, true, true);
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
    PopulateBCs(i);

    // -- upwind the scaled relative permeability kr xi eta_l / mu_l
    for (int c = 0; c < ncells_owned_; ++c) {
      kr_c[0][c] = u1_c[i][c] * eta_l_ * relperm_lc[0][c] / mu_l_;
    }
    kr->ViewComponent("dirichlet_faces")->PutScalar(1.0);  // FIXME
    upwind_->Compute(*kr, *kr, op_bcs_[0]->bc_model(), *kr);

    // -- init with transport in liquid phase
    auto& pdeK = pde_diff_K_;
    pdeK->Setup(Kptr, kr, Teuchos::null, rho_l_, gravity_);
    pdeK->SetBCs(op_bcs_[0], op_bcs_[0]);
    pdeK->global_operator()->Init();
    pdeK->UpdateMatrices(Teuchos::null, Teuchos::null);
    pdeK->ApplyBCs(true, true, true);
    pdeK->global_operator()->ComputeNegativeResidual(*u0, *f1);

    // -- upwind the scaled relative permeability for gas phase
    for (int c = 0; c < ncells_owned_; ++c) {
      kr_c[0][c] = (1.0 - u1_c[i][c]) * eta_gc[0][c] * relperm_gc[0][c] / mu_g_;
    }
    kr->ViewComponent("dirichlet_faces")->PutScalar(0.0);  // FIXME
    upwind_->Compute(*kr, *kr, op_bcs_[0]->bc_model(), *kr);

    // -- add transport in gas phase
    pdeK->Setup(Kptr, kr, Teuchos::null, rho_g, gravity_);
    pdeK->SetBCs(op_bc_pg_, op_bc_pg_);
    pdeK->global_operator()->Init();
    pdeK->UpdateMatrices(Teuchos::null, Teuchos::null);
    pdeK->ApplyBCs(true, true, true);
    pdeK->global_operator()->ComputeNegativeResidual(*pg, fadd);
    f1->Update(1.0, fadd, 1.0);

    // -- add molecular diffusion using harmonic-mean transmissibility
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
    pdeD->UpdateMatrices(Teuchos::null, Teuchos::null);
    pdeD->ApplyBCs(true, true, true);
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
  auto f2 = f->SubVector(2)->Data();
  auto& f2_c = *f2->ViewComponent("cell");

  const auto& u2 = u_new->SubVector(2)->Data();
  auto& u2_c = *u2->ViewComponent("cell");  // saturation liquid

  if (ncp_ == "min") {
    for (int c = 0; c < ncells_owned_; ++c) {
      double x_sum(0.0);
      for (int i = 0; i < num_primary_; ++i) {
        x_sum += u1_c[i][c];
      } 
      f2_c[0][c] = std::min(1.0 - u2_c[0][c], 1.0 - x_sum);
    }
  } else if (ncp_ == "Fischer-Burmeister") {
    for (int c = 0; c < ncells_owned_; ++c) {
      double x_sum(0.0);
      for (int i = 0; i < num_primary_; ++i) {
        x_sum += u1_c[i][c];
      } 
      double a = 1.0 - u2_c[0][c];
      double b = 1.0 - x_sum;
      f2_c[0][c] = pow(a * a + b * b, 0.5) - (a + b);
    }
  }
}


/* ******************************************************************
* Update new preconditioner on the interval (tp-dtp, tp].
****************************************************************** */
void Multiphase_PK::UpdatePreconditioner(double tp, Teuchos::RCP<const TreeVector> u, double dtp)
{
  double t_old = tp - dtp;

  PopulateBCs(0);

  // trigger update of primary variables
  saturation_liquid_eval_->SetFieldAsChanged(S_.ptr());

  // -----------------------
  // process water component
  // -----------------------
  const auto& u0 = u->SubVector(0)->Data();

  // accumulation
  // -- upwind relative permeability
  S_->GetFieldEvaluator(relperm_liquid_key_)->HasFieldChanged(S_.ptr(), passwd_);
  auto relperm_l = S_->GetFieldData(relperm_liquid_key_, relperm_liquid_key_);

  auto kr = CreateCVforUpwind(mesh_);
  *kr->ViewComponent("cell") = *relperm_l->ViewComponent("cell");
  kr->ViewComponent("dirichlet_faces")->PutScalar(1.0);  // FIXME
  upwind_->Compute(*kr, *kr, op_bcs_[0]->bc_model(), *kr);

  // -- scale the upwinded relperm assuming constant viscosity
  double factor = eta_l_ / mu_l_;
  kr->Scale(factor);

  // -- init with Darcy operator
  auto& pdeK = pde_diff_K_;
  op_preconditioner_->SetOperatorBlock(0, 0, pdeK->global_operator());

  Teuchos::RCP<std::vector<WhetStone::Tensor> > Kptr = Teuchos::rcpFromRef(K_);
  pdeK->Setup(Kptr, kr, Teuchos::null, rho_l_, gravity_);
  pdeK->SetBCs(op_bcs_[0], op_bcs_[0]);
  pdeK->global_operator()->Init();
  pdeK->UpdateMatrices(Teuchos::null, Teuchos::null);
  pdeK->ApplyBCs(true, true, true);

  auto flux_liquid = S_->GetFieldData(darcy_flux_liquid_key_, passwd_);
  pdeK->UpdateFlux(u0.ptr(), flux_liquid.ptr());

  // ---------------------------
  // process chemical components
  // ---------------------------
  std::vector<Teuchos::RCP<Operators::PDE_AdvectionUpwind> > pdeU(num_primary_);
  std::vector<Teuchos::RCP<Operators::PDE_Accumulation> > pdeA(num_primary_);

  auto& adv_list = mp_list_->sublist("operators").sublist("advection operator").sublist("preconditioner");

  for (int i = 0; i < num_primary_; ++i) {
    PopulateBCs(i);

    // -- init with upwind operator 
    pdeU[i] = Teuchos::rcp(new Operators::PDE_AdvectionUpwind(adv_list, mesh_));
    op_preconditioner_->SetOperatorBlock(i + 1, i + 1, pdeU[i]->global_operator());

    pdeU[i]->Setup(*flux_liquid);
    pdeU[i]->SetBCs(op_bcs_[1], op_bcs_[1]);
    pdeU[i]->global_operator()->Init();
    pdeU[i]->UpdateMatrices(flux_liquid.ptr());
    pdeU[i]->ApplyBCs(true, true, true);

    if (dtp > 0.0) {
      pdeA[i] = Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::CELL, pdeU[i]->global_operator())); 

      S_->GetFieldEvaluator(tcs_key_)->HasFieldDerivativeChanged(S_.ptr(), passwd_, x_liquid_key_);
      CompositeVector& dtcs_dpl = *S_->GetFieldData("dtotal_component_storage_dmolar_fraction_liquid", tcs_key_);

      pdeA[i]->AddAccumulationTerm(dtcs_dpl, dtp, "cell", true);
    }
  }
 
  // ------------------
  // process constraint
  // ------------------
  auto pdeC = Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::CELL, mesh_)); 
  op_preconditioner_->SetOperatorBlock(num_primary_ + 1, num_primary_ + 1, pdeC->global_operator());
 
  // -- reuse temporary vector
  auto& kr_c = *kr->ViewComponent("cell");
  for (int c = 0; c < ncells_owned_; ++c) {
    kr_c[0][c] = 1.0;
  }
  pdeC->AddAccumulationTerm(*kr, "cell");

  // finalize preconditioner
  if (!op_pc_assembled_) {
    op_preconditioner_->SymbolicAssembleMatrix();
    op_pc_assembled_ = true;
  }
  op_preconditioner_->AssembleMatrix();
  op_preconditioner_->UpdatePreconditioner();

  // cleaning temporaty operators
  for (int i = 0; i < num_primary_; ++i) {
    pdeU[i] = Teuchos::null;
    pdeA[i] = Teuchos::null;
  }
}


/* ******************************************************************
* Apply preconditioner inv(B) * X.                                                 
****************************************************************** */
int Multiphase_PK::ApplyPreconditioner(Teuchos::RCP<const TreeVector> X, 
                                       Teuchos::RCP<TreeVector> Y)
{
  Y->PutScalar(0.0);
  return op_pc_solver_->ApplyInverse(*X, *Y);
}


/* ******************************************************************
* Monitor l2 norm of residual
****************************************************************** */
double Multiphase_PK::ErrorNorm(Teuchos::RCP<const TreeVector> u,
                                Teuchos::RCP<const TreeVector> du) 
{
  // pressure error
  auto pc = *u->SubVector(0)->Data()->ViewComponent("cell");
  auto dpc = *du->SubVector(0)->Data()->ViewComponent("cell");

  double atm_pressure(1.0e+5);

  double error_p = 0.0;
  for (int c = 0; c < ncells_owned_; c++) {
    double tmp = fabs(dpc[0][c]) / (fabs(pc[0][c] - atm_pressure) + atm_pressure);
    if (tmp > error_p) {
      error_p = tmp;
    } 
  }

  // saturation error
  auto dsc = *du->SubVector(2)->Data()->ViewComponent("cell");

  double error_s = 0.0;
  for (int c = 0; c < ncells_owned_; c++) {
    error_s  = std::max(error_s, fabs(dsc[0][c]));
  }

  return error_p + error_s;
}


/********************************************************************
* Modifies nonlinear update du using .. TBW
****************************************************************** */
AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
Multiphase_PK::ModifyCorrection(double h, Teuchos::RCP<const TreeVector> res,
                                Teuchos::RCP<const TreeVector> u,
                                Teuchos::RCP<TreeVector> du)
{
  auto u2 = u->SubVector(2);
  auto du2 = du->SubVector(2);

  auto sat_next = Teuchos::rcp(new CompositeVector(*u2->Data()));
  sat_next->Update(-1.0, *du2->Data(), 1.0);

  Epetra_MultiVector& sat_c = *sat_next->ViewComponent("cell");
  for (int c = 0; c < ncells_owned_; ++c) {
    sat_c[0][c] = std::max(0.0, sat_c[0][c]);
  }

  du2->Data()->Update(1.0, *u2->Data(), -1.0, *sat_next, 0.0);

  return AmanziSolvers::FnBaseDefs::CORRECTION_MODIFIED;
}

}  // namespace Multiphase
}  // namespace Amanzi

