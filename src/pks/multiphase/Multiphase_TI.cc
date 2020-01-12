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
  pressure_liquid_eval_->SetFieldAsChanged(S_.ptr());
  xl_liquid_eval_->SetFieldAsChanged(S_.ptr());
  saturation_liquid_eval_->SetFieldAsChanged(S_.ptr());

  // extract pointers to subvectors
  const auto& u0 = u_new->SubVector(0)->Data();
  const auto& u1 = u_new->SubVector(1)->Data();
  const auto& u2 = u_new->SubVector(2)->Data();

  auto f0 = f->SubVector(0)->Data();
  auto f1 = f->SubVector(1)->Data();
  auto f2 = f->SubVector(2)->Data();

  auto& f0c = *f0->ViewComponent("cell");
  auto& f1c = *f1->ViewComponent("cell");
  auto& f2c = *f2->ViewComponent("cell");

  // miscalleneous fields
  auto flux_liquid = S_->GetFieldData(darcy_flux_liquid_key_, passwd_);
  auto flux_gas = S_->GetFieldData(darcy_flux_gas_key_, passwd_);

  // -------------------------
  // process pressure equation
  // -------------------------
  // -- upwind relative permeability
  S_->GetFieldEvaluator(relperm_liquid_key_)->HasFieldChanged(S_.ptr(), passwd_);
  auto relperm_l = S_->GetFieldData(relperm_liquid_key_, relperm_liquid_key_);

  auto kr = CreateCVforUpwind(mesh_);
  *kr->ViewComponent("cell") = *relperm_l->ViewComponent("cell");
  kr->ViewComponent("dirichlet_faces")->PutScalar(1.0);  // FIXME
  upwind_->Compute(*flux_liquid, *kr, op_bcs_[0]->bc_model(), *kr);

  // -- scale the upwinded relperm assuming constant viscosity
  double factor = eta_l_ / mu_l_;
  kr->Scale(factor);

  // -- create operator
  Teuchos::RCP<std::vector<WhetStone::Tensor> > Kptr = Teuchos::rcpFromRef(K_);

  auto& pde = pde_diff_K_;
  pde->Setup(Kptr, kr, Teuchos::null, rho_l_, gravity_);
  pde->SetBCs(op_bcs_[0], op_bcs_[0]);
  pde->global_operator()->Init();
  pde->UpdateMatrices(Teuchos::null, Teuchos::null);
  pde->ApplyBCs(true, false, true);
  pde->global_operator()->ComputeNegativeResidual(*u0, *f0);
 
  // -----------------------------------------
  // process equations for chemical components
  // -----------------------------------------
  const auto& u2c = *u2->ViewComponent("cell");
  CompositeVector fone(*f0), fadd(*f0), u2i(*f0);

  // -- gas pressure
  S_->GetFieldEvaluator(pressure_gas_key_)->HasFieldChanged(S_.ptr(), passwd_);
  auto pg = S_->GetFieldData(pressure_gas_key_, pressure_gas_key_);

  // -- mass density of gas phase (create evaluator ??? FIXME)
  const auto& eta_g = S_->GetFieldData("molar_density_gas");
  auto rho_g = Teuchos::rcp(new CompositeVector(*eta_g));

  // -- update relative premeability of gas phase
  S_->GetFieldEvaluator(relperm_gas_key_)->HasFieldChanged(S_.ptr(), passwd_);
  auto relperm_g = S_->GetFieldData(relperm_gas_key_, relperm_gas_key_);

  const auto& eta_gc = *eta_g->ViewComponent("cell");
  const auto& relperm_lc = *relperm_l->ViewComponent("cell");
  const auto& relperm_gc = *relperm_g->ViewComponent("cell");
  auto& kr_c = *kr->ViewComponent("cell");

  for (int i = 0; i < num_primary_; ++i) {
    PopulateBCs(i);

    // -- upwind the scaled relative permeability for liquid phase, xi eta kr / mu
    for (int c = 0; c < ncells_owned_; ++c) {
      kr_c[0][c] = u2c[i][c] * eta_l_ * relperm_lc[0][c] / mu_l_;
    }
    kr->ViewComponent("dirichlet_faces")->PutScalar(eta_l_ / mu_l_);  // FIXME
    upwind_->Compute(*flux_liquid, *kr, op_bcs_[0]->bc_model(), *kr);

    // -- init with transport in liquid phase
    auto& pdeK = pde_diff_K_;
    pdeK->Setup(Kptr, kr, Teuchos::null, rho_l_, gravity_);
    pdeK->SetBCs(op_bcs_[2], op_bcs_[0]);
    pdeK->global_operator()->Init();
    pdeK->UpdateMatrices(Teuchos::null, Teuchos::null);
    pdeK->ApplyBCs(true, false, true);
    pdeK->global_operator()->ComputeNegativeResidual(*u0, fone);

    // -- upwind the scaled relative permeability for gas phase, xi eta kr / mu
    for (int c = 0; c < ncells_owned_; ++c) {
      kr_c[0][c] = (1.0 - u2c[i][c]) * eta_gc[0][c] * relperm_gc[0][c] / mu_g_;
    }
    kr->ViewComponent("dirichlet_faces")->PutScalar(0.0);  // FIXME
    upwind_->Compute(*flux_gas, *kr, op_bcs_[0]->bc_model(), *kr);

    // -- add transport in gas phase
    pdeK->Setup(Kptr, kr, Teuchos::null, rho_g, gravity_);
    pdeK->SetBCs(op_bcs_[2], op_bc_pg_);
    pdeK->global_operator()->Init();
    pdeK->UpdateMatrices(Teuchos::null, Teuchos::null);
    pdeK->ApplyBCs(false, false, false);
    pdeK->global_operator()->ComputeNegativeResidual(*pg, fadd);
    fone.Update(1.0, fadd, 1.0);

    // -- add molecular diffusion using harmonic-mean transmissibility
    *u2i.ViewComponent("cell") = *u2c(i);

    auto D = Teuchos::rcp(new std::vector<WhetStone::Tensor>);
    WhetStone::Tensor Dc(dim_, 1);
    for (int c = 0; c < ncells_owned_; ++c) {
      Dc(0, 0) = eta_l_ * mol_diff_l_[i] - eta_gc[0][c] * mol_diff_g_[i];
      D->push_back(Dc);
    }

    auto& pdeD = pde_diff_D_;
    pdeD->Setup(D, Teuchos::null, Teuchos::null);
    pdeK->SetBCs(op_bcs_[2], op_bcs_[2]);
    pdeD->global_operator()->Init();
    pdeD->UpdateMatrices(Teuchos::null, Teuchos::null);
    pdeD->ApplyBCs(false, false, false);
    pdeD->global_operator()->ComputeNegativeResidual(u2i, fadd);
    fone.Update(1.0, fadd, 1.0);

    // -- add counter-diffusion to the water component
    f0->Update(-1.0, fadd, 1.0);

    // -- save residual in the proper subvector
    if (i == 0)  {
      *f1 = fone;
    } else {
      const auto& fone_c = *fone.ViewComponent("cell");
      for (int c = 0; c < ncells_owned_; ++c) f2c[i - 1][c] = fone_c[0][c];
    }
  }

  // add accumulation terms to equations for pressure and chemical components
  {
    std::string prev_name = "prev_" + tws_key_;
    S_->GetFieldEvaluator(tws_key_)->HasFieldChanged(S_.ptr(), passwd_);
    const auto& total_c = *S_->GetFieldData(tws_key_)->ViewComponent("cell");
    const auto& total_prev_c = *S_->GetFieldData(prev_name)->ViewComponent("cell");

    for (int c = 0; c < ncells_owned_; ++c) {
      double factor = mesh_->cell_volume(c) / dtp;
      f0c[0][c] += (total_c[0][c] - total_prev_c[0][c]) * factor;
    }
  }

  for (int i = 0; i < num_primary_; ++i) {
    std::string prev_name = "prev_" + tcs_key_;

    auto eval = S_->GetFieldEvaluator(tcs_key_);
    Teuchos::rcp_dynamic_cast<TotalComponentStorage>(eval)->set_component_id(i);
    eval->HasFieldChanged(S_.ptr(), passwd_);

    const auto& total_c = *S_->GetFieldData(tcs_key_)->ViewComponent("cell");
    const auto& total_prev_c = *S_->GetFieldData(prev_name)->ViewComponent("cell");

    if (i == 0) {
      for (int c = 0; c < ncells_owned_; ++c) {
        double factor = mesh_->cell_volume(c) / dtp;
        f1c[0][c] += (total_c[0][c] - total_prev_c[0][c]) * factor;
      }
    } else {
      for (int c = 0; c < ncells_owned_; ++c) {
        double factor = mesh_->cell_volume(c) / dtp;
        f2c[i - 1][c] += (total_c[0][c] - total_prev_c[0][c]) * factor;
      }
    }
  }

  // ------------------
  // process constraint
  // ------------------
  auto& u1c = *u1->ViewComponent("cell");  // saturation liquid

  int n = num_primary_ - 1;
  if (ncp_ == "min") {
    for (int c = 0; c < ncells_owned_; ++c) {
      double x_sum(0.0);
      for (int i = 0; i < num_primary_; ++i) x_sum += u2c[i][c];

      f2c[n][c] = std::min(1.0 - u1c[0][c], 1.0 - x_sum);
    }
  } else if (ncp_ == "Fischer-Burmeister") {
    for (int c = 0; c < ncells_owned_; ++c) {
      double x_sum(0.0);
      for (int i = 0; i < num_primary_; ++i) x_sum += u2c[i][c];

      double a = 1.0 - u1c[0][c];
      double b = 1.0 - x_sum;
      f2c[n][c] = pow(a * a + b * b, 0.5) - (a + b);
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

  // extract pointers to subvectors
  const auto& u0 = u->SubVector(0)->Data();
  const auto& u1 = u->SubVector(1)->Data();
  const auto& u2 = u->SubVector(2)->Data();

  // parameter lists
  auto& adv_list = mp_list_->sublist("operators").sublist("advection operator").sublist("preconditioner");

  // miscalleneous fields
  auto flux_liquid = S_->GetFieldData(darcy_flux_liquid_key_, passwd_);
  auto flux_gas = S_->GetFieldData(darcy_flux_gas_key_, passwd_);

  // -----------------------
  // process water component
  // -----------------------
  // accumulation
  // -- upwind relative permeability
  S_->GetFieldEvaluator(relperm_liquid_key_)->HasFieldChanged(S_.ptr(), passwd_);
  auto relperm_l = S_->GetFieldData(relperm_liquid_key_, relperm_liquid_key_);

  auto kr = CreateCVforUpwind(mesh_);
  *kr->ViewComponent("cell") = *relperm_l->ViewComponent("cell");
  kr->ViewComponent("dirichlet_faces")->PutScalar(1.0);  // FIXME
  upwind_->Compute(*flux_liquid, *kr, op_bcs_[0]->bc_model(), *kr);

  // -- scale the upwinded relperm assuming constant viscosity
  double factor = eta_l_ / mu_l_;
  kr->Scale(factor);

  // -- block (0, 0): Darcy operator
  //    kr is used by matrix update, so we can reuse this variable later
  //    Darcy flux is made consistent with kr and pressure
  auto& pdeK = pde_diff_K_;
  op_preconditioner_->SetOperatorBlock(0, 0, pdeK->global_operator());

  Teuchos::RCP<std::vector<WhetStone::Tensor> > Kptr = Teuchos::rcpFromRef(K_);
  pdeK->Setup(Kptr, kr, Teuchos::null, rho_l_, gravity_);
  pdeK->SetBCs(op_bcs_[0], op_bcs_[0]);
  pdeK->global_operator()->Init();
  pdeK->UpdateMatrices(Teuchos::null, Teuchos::null);
  pdeK->ApplyBCs(true, true, true);

  pdeK->UpdateFlux(u0.ptr(), flux_liquid.ptr());
  flux_liquid->Scale(1.0 / eta_l_);

  // -- block (0, 1)
  // --- advection term
  auto pde01a = Teuchos::rcp(new Operators::PDE_AdvectionUpwind(adv_list, mesh_)); 
  op_preconditioner_->SetOperatorBlock(0, 1, pde01a->global_operator());

  S_->GetFieldEvaluator(relperm_liquid_key_)->HasFieldDerivativeChanged(S_.ptr(), passwd_, saturation_liquid_key_);
  auto dkrl_dsl = S_->GetFieldData("drel_permeability_liquid_dsaturation_liquid", relperm_liquid_key_);

  *kr->ViewComponent("cell") = *dkrl_dsl->ViewComponent("cell");
  kr->ViewComponent("dirichlet_faces")->PutScalar(0.0);  // FIXME
  upwind_->Compute(*flux_liquid, *kr, op_bcs_[0]->bc_model(), *kr);

  auto flux_tmp(*flux_liquid);
  auto& flux_f = *flux_tmp.ViewComponent("face");
  auto& kr_f = *kr->ViewComponent("face");
  for (int f = 0; f < nfaces_owned_; ++f) flux_f[0][f] /= kr_f[0][f];

  pde01a->Setup(flux_tmp);
  pde01a->SetBCs(op_bcs_[0], op_bcs_[1]);
  pde01a->global_operator()->Init();
  pde01a->UpdateMatrices(flux_liquid.ptr());
  pde01a->ApplyBCs(false, false, false);

  // --- accumulation term
  auto pde01b = Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::CELL, pde01a->global_operator())); 

  S_->GetFieldEvaluator(tws_key_)->HasFieldDerivativeChanged(S_.ptr(), passwd_, saturation_liquid_key_);
  auto dtws_dsl = S_->GetFieldData("dtotal_water_storage_dsaturation_liquid", tws_key_);

  auto& kr_c = *kr->ViewComponent("cell");
  auto& dtws_dsl_c = *dtws_dsl->ViewComponent("cell");
  for (int c = 0; c < ncells_owned_; ++c) {
    kr_c[0][c] = dtws_dsl_c[0][c] / dtp;
  }
  pde01b->AddAccumulationTerm(*kr, "cell");

  // -- block (0, 2)

  // ---------------------------
  // process chemical components
  // ---------------------------
  std::vector<Teuchos::RCP<Operators::PDE_AdvectionUpwind> > pde12a(num_primary_), pde11a(num_primary_);
  std::vector<Teuchos::RCP<Operators::PDE_Accumulation> > pde11b(num_primary_), pde12b(num_primary_);
  std::vector<Teuchos::RCP<Operators::PDE_Accumulation> > pde10a(num_primary_);

  for (int i = 0; i < num_primary_; ++i) {
    PopulateBCs(i);

    // -- block (1, 1)
    if (i == 0) {
      // --- advection
      pde11a[i] = Teuchos::rcp(new Operators::PDE_AdvectionUpwind(adv_list, mesh_)); 
      op_preconditioner_->SetOperatorBlock(i + 1, 1, pde11a[i]->global_operator());

      pde11a[i]->Setup(flux_tmp);
      pde11a[i]->SetBCs(op_bcs_[2], op_bcs_[2]);
      pde11a[i]->global_operator()->Init();
      pde11a[i]->UpdateMatrices(flux_liquid.ptr());
      pde11a[i]->ApplyBCs(true, false, false);

      // --- accumulation
      auto eval = S_->GetFieldEvaluator(tcs_key_);
      Teuchos::rcp_dynamic_cast<TotalComponentStorage>(eval)->set_component_id(i);
      eval->HasFieldDerivativeChanged(S_.ptr(), passwd_, saturation_liquid_key_);
      CompositeVector& dtcs_dsl = *S_->GetFieldData("dtotal_component_storage_dsaturation_liquid", tcs_key_);

      pde11b[i] = Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::CELL, pde11a[i]->global_operator())); 
      pde11b[i]->AddAccumulationTerm(dtcs_dsl, dtp, "cell", true);
    }

    // -- block (1, 2)
    // --- init with upwind operator 
    pde12a[i] = Teuchos::rcp(new Operators::PDE_AdvectionUpwind(adv_list, mesh_));
    op_preconditioner_->SetOperatorBlock(i + 1, i + 2, pde12a[i]->global_operator());

    pde12a[i]->Setup(*flux_liquid);
    pde12a[i]->SetBCs(op_bcs_[2], op_bcs_[2]);
    pde12a[i]->global_operator()->Init();
    pde12a[i]->UpdateMatrices(flux_liquid.ptr());
    pde12a[i]->ApplyBCs(false, false, false);

    // --- add storage term 
    if (dtp > 0.0) {
      auto eval = S_->GetFieldEvaluator(tcs_key_);
      Teuchos::rcp_dynamic_cast<TotalComponentStorage>(eval)->set_component_id(i);
      eval->HasFieldDerivativeChanged(S_.ptr(), passwd_, x_liquid_key_);
      CompositeVector& dtcs_dxl = *S_->GetFieldData("dtotal_component_storage_dmolar_fraction_liquid", tcs_key_);

      pde12b[i] = Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::CELL, pde12a[i]->global_operator())); 
      pde12b[i]->AddAccumulationTerm(dtcs_dxl, dtp, "cell", true);
    }

    // -- bloack (1, 0)
    if (dtp > 0.0) {
      auto eval = S_->GetFieldEvaluator(tcs_key_);
      Teuchos::rcp_dynamic_cast<TotalComponentStorage>(eval)->set_component_id(i);
      eval->HasFieldDerivativeChanged(S_.ptr(), passwd_, pressure_liquid_key_);
      CompositeVector& dtcs_dpl = *S_->GetFieldData("dtotal_component_storage_dpressure_liquid", tcs_key_);

      pde10a[i] = Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::CELL, mesh_)); 
      pde10a[i]->AddAccumulationTerm(dtcs_dpl, dtp, "cell", true);
    }
  }
 
  // ---------------------------------
  // process constraint: blocks (2, *)
  // ---------------------------------
  auto& u1c = *u1->ViewComponent("cell");  // saturation liquid
  auto& u2c = *u2->ViewComponent("cell");  // chemical componets

  CompositeVector dfdx(*u1), dfds(*u1);
  auto& dfdx_c = *dfdx.ViewComponent("cell");
  auto& dfds_c = *dfds.ViewComponent("cell");

  std::vector<Teuchos::RCP<Operators::PDE_Accumulation> > pde22a(num_primary_);

  for (int i = 0; i < num_primary_; ++i) {
    pde22a[i] = Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::CELL, mesh_)); 
    op_preconditioner_->SetOperatorBlock(num_primary_ + 1, i + 2, pde22a[i]->global_operator());
 
    // -- identify active set for gas phase
    std::vector<int> active_g(ncells_owned_, 0);
    std::vector<int> inactive_g(ncells_owned_, 1);

    for (int c = 0; c < ncells_owned_; c++) {
      double x_sum(0.0);
      for (int i = 0; i < num_primary_; ++i) x_sum += u2c[i][c];

      if (std::fabs(1.0 - u1c[0][c]) - std::fabs(1.0 - x_sum) > 1e-12) {
        active_g[c] = 1;
        inactive_g[c] = 0;
      }
    }

    if (ncp_ == "min") {
      for (int c = 0; c < ncells_owned_; c++) {
        double volume = mesh_->cell_volume(c);
        if (inactive_g[c] == 1) {
          dfdx_c[0][c] = 0.0;
          dfds_c[0][c] =-1.0;
        } else {
          dfdx_c[0][c] =-1.0;
          dfds_c[0][c] = 0.0;
        }
      }
    }
    pde22a[i]->AddAccumulationTerm(dfdx, "cell");
  }

  auto pde21 = Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::CELL, mesh_)); 
  op_preconditioner_->SetOperatorBlock(num_primary_ + 1, 1, pde21->global_operator());
  pde21->AddAccumulationTerm(dfds, "cell");

  // finalize preconditioner
  if (!op_pc_assembled_) {
    op_preconditioner_->SymbolicAssembleMatrix();
    op_pc_assembled_ = true;
  }
  op_preconditioner_->AssembleMatrix();
  op_preconditioner_->UpdatePreconditioner();

  // cleaning temporaty operators
  for (int i = 0; i < num_primary_; ++i) {
    pde10a[i] = Teuchos::null;
    pde11a[i] = Teuchos::null;
    pde11b[i] = Teuchos::null;
    pde12a[i] = Teuchos::null;
    pde12b[i] = Teuchos::null;
    pde22a[i] = Teuchos::null;
  }
}


/* ******************************************************************
* Apply preconditioner inv(B) * X.                                                 
****************************************************************** */
int Multiphase_PK::ApplyPreconditioner(Teuchos::RCP<const TreeVector> X, 
                                       Teuchos::RCP<TreeVector> Y)
{
  Y->PutScalar(0.0);
  // *Y = *X; return 0;
  int ierr = op_pc_solver_->ApplyInverse(*X, *Y);
// { double aaa; X->Norm2(&aaa); std::cout << aaa << std::endl;}
// { double aaa; Y->Norm2(&aaa); std::cout << aaa << std::endl;}
// std::cout << *op_preconditioner_->A() << std::endl; exit(0);
  return ierr;
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
  // clip mole fraction to range [0; 1]
  const auto& u2c = *u->SubVector(2)->Data()->ViewComponent("cell");
  auto& du2c = *du->SubVector(2)->Data()->ViewComponent("cell");

  for (int i = 0; i < u2c.NumVectors(); ++i) {
    for (int c = 0; c < ncells_owned_; ++c) {
      du2c[i][c] = std::min(du2c[i][c], u2c[i][c]);
      du2c[i][c] = std::max(du2c[i][c], u2c[i][c] - 1.0);
    }    
  }

  // clip saturation (residual saturation is missing, FIXME)
  const auto& u1c = *u->SubVector(1)->Data()->ViewComponent("cell");
  auto& du1c = *du->SubVector(1)->Data()->ViewComponent("cell");

  for (int c = 0; c < ncells_owned_; ++c) {
    du1c[0][c] = std::min(du1c[0][c], u1c[0][c]);
    du1c[0][c] = std::max(du1c[0][c], u1c[0][c] - 1.0);
  }

  return AmanziSolvers::FnBaseDefs::CORRECTION_MODIFIED;
}

}  // namespace Multiphase
}  // namespace Amanzi

