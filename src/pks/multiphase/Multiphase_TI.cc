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

  // extract pointers to subvectors
  std::vector<Teuchos::RCP<CompositeVector> > up, fp;
  for (int i = 0; i < 3; ++i) {
    up.push_back(u_new->SubVector(i)->Data());
    fp.push_back(f->SubVector(i)->Data());
  }

  // miscalleneous fields
  // -- saturation
  auto& sat_lc = *S_->GetFieldData(saturation_liquid_key_, passwd_)->ViewComponent("cell");

  // -- gas pressure
  S_->GetFieldEvaluator(pressure_gas_key_)->HasFieldChanged(S_.ptr(), passwd_);
  auto pg = S_->GetFieldData(pressure_gas_key_, pressure_gas_key_);

  // -- molar densities
  S_->GetFieldEvaluator(molar_density_gas_key_)->HasFieldChanged(S_.ptr(), passwd_);
  const auto& eta_g = S_->GetFieldData(molar_density_gas_key_);
  const auto& eta_gc = *eta_g->ViewComponent("cell");

  const auto& eta_l = S_->GetFieldData(molar_density_liquid_key_);
  const auto& eta_lc = *eta_l->ViewComponent("cell");

  // -- storage
  S_->GetFieldEvaluator(tws_key_)->HasFieldChanged(S_.ptr(), passwd_);
  S_->GetFieldEvaluator(tcs_key_)->HasFieldChanged(S_.ptr(), passwd_);

  // -- porosity
  const auto& phi = *S_->GetFieldData(porosity_key_)->ViewComponent("cell");

  // -- wrapper for absolute permeability
  Teuchos::RCP<std::vector<WhetStone::Tensor> > Kptr = Teuchos::rcpFromRef(K_);

  // work memory for miscalleneous operator
  auto kr = CreateCVforUpwind(mesh_);
  auto& kr_c = *kr->ViewComponent("cell");
  auto& kr_f = *kr->ViewComponent("face");

  // primary variables
  auto tmp = up[0];
  CompositeVector fone(*tmp), fadd(*tmp), comp(*tmp);
  auto& fone_c = *fone.ViewComponent("cell");
  auto& comp_c = *comp.ViewComponent("cell");

  // start loop over physical equations
  Key key;
  for (int n = 0; n < num_primary_ + 1; ++n) {
    ModifyEvaluators(n);
    auto eqn = EquationToSolution(n);
    PopulateBCs(eqn.second, true);
  
    // Richards-type operator for all phases
    fone.PutScalar(0.0);
    std::vector<std::string> varp_name{pressure_liquid_key_, pressure_gas_key_};
    std::vector<std::string> flux_name{darcy_flux_liquid_key_, darcy_flux_gas_key_};

    for (int phase = 0; phase < 2; ++phase) {
      bool bcflag = (phase == 0);
      if ((key = eval_eqns_[n][phase].first) != "") {
        S_->GetFieldEvaluator(key)->HasFieldChanged(S_.ptr(), passwd_);

        // -- upwind cell-centered coefficient
        auto flux = S_->GetFieldData(flux_name[phase], passwd_);
        kr_c = *S_->GetFieldData(key)->ViewComponent("cell");
        kr->ViewComponent("dirichlet_faces")->PutScalar(0.0);  // FIXME (BC data)
        upwind_->Compute(*flux, *kr, op_bcs_[eqn.first]->bc_model(), *kr);

        // -- form operator
        auto& pde = pde_diff_K_;
        pde->Setup(Kptr, kr, Teuchos::null, rho_l_, gravity_);  // FIXME (gravity for gas phase)
        pde->SetBCs(op_bcs_[eqn.first], op_bcs_[eqn.first]);
        pde->global_operator()->Init();
        pde->UpdateMatrices(Teuchos::null, Teuchos::null);
        pde->ApplyBCs(bcflag, false, false);

        // -- add advection term to the residual
        S_->GetFieldEvaluator(varp_name[phase])->HasFieldChanged(S_.ptr(), passwd_);
        auto var = S_->GetFieldData(varp_name[phase]);
        pde->global_operator()->ComputeNegativeResidual(*var, fadd);
        fone.Update(1.0, fadd, 1.0);
      }
    }
 
    // molecular diffusion via harmonic-mean transmissibility
    for (int phase = 0; phase < 2; ++phase) {
      if ((key = eval_eqns_[n][2 + phase].first) != "") {
        double coef;
        auto D = Teuchos::rcp(new std::vector<WhetStone::Tensor>);
        WhetStone::Tensor Dc(dim_, 1);
        for (int c = 0; c < ncells_owned_; ++c) {
          // if (phase == 0) coef = sat_lc[0][c] * eta_lc[0][c] * mol_diff_l_[eqn.second];
          if (phase == 0) coef = sat_lc[0][c] * mol_diff_l_[eqn.second];
          if (phase == 1) coef = (1.0 - sat_lc[0][c]) * eta_gc[0][c] * mol_diff_g_[eqn.second];
          Dc(0, 0) = phi[0][c] * coef;
          D->push_back(Dc);
        }

        // -- form operator
        auto& pde = pde_diff_D_;
        pde->Setup(D, Teuchos::null, Teuchos::null);
        pde->SetBCs(op_bcs_[eqn.first], op_bcs_[eqn.first]);
        pde->global_operator()->Init();
        pde->UpdateMatrices(Teuchos::null, Teuchos::null);
        pde->ApplyBCs(false, false, false);

        // -- add diffusion term to the residual
        S_->GetFieldEvaluator(varx_name_[phase])->HasFieldChanged(S_.ptr(), passwd_);
        auto& tmp = *S_->GetFieldData(varx_name_[phase])->ViewComponent("cell");
        int m = std::max(n - 1, tmp.NumVectors() - 1);
        for (int c = 0; c < ncells_owned_; ++c) {
          comp_c[0][c] = tmp[m][c];
        }
        pde->global_operator()->ComputeNegativeResidual(comp, fadd);

        double factor = eval_eqns_[n][2 + phase].second;
        fone.Update(factor, fadd, 1.0);
      }
    }

    // add storage terms 
    if ((key = eval_eqns_[n][4].first) != "") {
      std::string prev_key = "prev_" + key;
      S_->GetFieldEvaluator(key)->HasFieldChanged(S_.ptr(), passwd_);

      const auto& total_c = *S_->GetFieldData(key)->ViewComponent("cell");
      const auto& total_prev_c = *S_->GetFieldData(prev_key)->ViewComponent("cell");

      for (int c = 0; c < ncells_owned_; ++c) {
        double factor = mesh_->cell_volume(c) / dtp;
        fone_c[0][c] += (total_c[0][c] - total_prev_c[0][c]) * factor;
      }
    }

    // copy temporaty vector to residual
    auto& fc = *fp[eqn.first]->ViewComponent("cell");
    for (int c = 0; c < ncells_owned_; ++c)
      fc[eqn.second][c] = fone_c[0][c];
  }

  // process gas constraints
  int n = num_primary_ + 1;
  key = eval_eqns_[n][0].first;
  S_->GetFieldEvaluator(key)->HasFieldChanged(S_.ptr(), passwd_);
  const auto& ncp_fc = *S_->GetFieldData(key)->ViewComponent("cell");

  key = eval_eqns_[n][1].first;
  S_->GetFieldEvaluator(key)->HasFieldChanged(S_.ptr(), passwd_);
  const auto& ncp_gc = *S_->GetFieldData(key)->ViewComponent("cell");

  auto& fci = *fp[2]->ViewComponent("cell");
  if (ncp_ == "min") {
    for (int c = 0; c < ncells_owned_; ++c) {
      fci[0][c] = std::min(ncp_fc[0][c], ncp_gc[0][c]);
    }
  } else if (ncp_ == "Fischer-Burmeister") {
    for (int c = 0; c < ncells_owned_; ++c) {
      double a = ncp_fc[0][c];
      double b = ncp_gc[0][c];
      fci[0][c] = pow(a * a + b * b, 0.5) - (a + b);
    }
  }
}


/* ******************************************************************
* Update new preconditioner on the interval (tp-dtp, tp].
****************************************************************** */
void Multiphase_PK::UpdatePreconditioner(double tp, Teuchos::RCP<const TreeVector> u, double dtp)
{
  double t_old = tp - dtp;

  // extract pointers to subvectors
  std::vector<Teuchos::RCP<const CompositeVector> > up;
  for (int i = 0; i < 3; ++i) {
    up.push_back(u->SubVector(i)->Data());
  }

  // miscalleneous fields
  // -- saturation
  auto& sat_lc = *S_->GetFieldData(saturation_liquid_key_, passwd_)->ViewComponent("cell");

  // -- porosity
  const auto& phi = *S_->GetFieldData(porosity_key_)->ViewComponent("cell");

  // -- molar densities
  S_->GetFieldEvaluator(molar_density_gas_key_)->HasFieldChanged(S_.ptr(), passwd_);
  const auto& eta_g = S_->GetFieldData(molar_density_gas_key_);
  const auto& eta_gc = *eta_g->ViewComponent("cell");

  const auto& eta_l = S_->GetFieldData(molar_density_liquid_key_);
  const auto& eta_lc = *eta_l->ViewComponent("cell");

  // -- mass density of gas phase 
  auto rho_g = Teuchos::rcp(new CompositeVector(*eta_g));

  // -- gas pressure
  S_->GetFieldEvaluator(pressure_gas_key_)->HasFieldChanged(S_.ptr(), passwd_);
  auto pg = S_->GetFieldData(pressure_gas_key_, pressure_gas_key_);

  // -- wrapper for absolute permeability
  Teuchos::RCP<std::vector<WhetStone::Tensor> > Kptr = Teuchos::rcpFromRef(K_);

  // parameter lists
  auto& adv_list = mp_list_->sublist("operators").sublist("advection operator").sublist("preconditioner");
  auto& ddf_list = mp_list_->sublist("operators").sublist("diffusion operator").sublist("preconditioner");
  auto& mdf_list = mp_list_->sublist("operators").sublist("molecular diffusion operator").sublist("preconditioner");

  // work memory for miscalleneous operator
  Key der_key;
  auto flux_tmp = Teuchos::rcp(new CompositeVector(*S_->GetFieldData(darcy_flux_liquid_key_)));

  auto kr = CreateCVforUpwind(mesh_);
  auto& kr_c = *kr->ViewComponent("cell");
  auto& kr_f = *kr->ViewComponent("face");

  CompositeVector fone(*up[0]);
  auto& fone_c = *fone.ViewComponent("cell");

  // for each operator we linearize (a) functions on which it acts and
  // (b) non-linear coefficients
  Key key;
  for (int row = 0; row < num_primary_ + 1; ++row) {
    ModifyEvaluators(row);
    auto eqnr = EquationToSolution(row);
    Key keyr = soln_names_[eqnr.first];
    PopulateBCs(eqnr.second, false);

    for (int col = 0; col < num_primary_ + 2; ++col) {
      auto eqnc = EquationToSolution(col);
      Key keyc = soln_names_[eqnc.first];

      // add empty operator to have a well-defined global operator pointer
      auto pde = Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::CELL, mesh_));
      auto global_op = pde->global_operator();
      op_preconditioner_->SetOperatorBlock(row, col, global_op);
      kr_c.PutScalar(0.0);
      pde->AddAccumulationTerm(*kr, "cell");

      // Richards operator for all phases
      std::vector<std::string> varp_name{pressure_liquid_key_, pressure_gas_key_};
      std::vector<std::string> flux_name{darcy_flux_liquid_key_, darcy_flux_gas_key_};

      for (int phase = 0; phase < 2; ++phase) {
        // -- diffusion operator div[ (K f du/dv) grad dv ] 
        if ((key = eval_eqns_[row][phase].first) != "") {
          if (varp_name[phase] == keyc) {
            der_key = "constant_field";  // DAG does not calculate derivative when u=v
          } else {
            der_key = "d" + varp_name[phase] + "_d" + keyc;
            S_->GetFieldEvaluator(varp_name[phase])->HasFieldDerivativeChanged(S_.ptr(), passwd_, keyc);
          }

          if (S_->HasField(der_key)) {
            auto pde = Teuchos::rcp(new Operators::PDE_DiffusionFV(ddf_list, global_op));

            S_->GetFieldEvaluator(key)->HasFieldChanged(S_.ptr(), passwd_);
            const auto& coef_c = *S_->GetFieldData(key)->ViewComponent("cell");
            const auto& der_c = *S_->GetFieldData(der_key)->ViewComponent("cell");

            for (int c = 0; c < ncells_owned_; ++c) {
              kr_c[0][c] = der_c[0][c] * coef_c[0][c];
            }
            kr->ViewComponent("dirichlet_faces")->PutScalar(0.0);  // FIXME
            auto flux = *S_->GetFieldData(flux_name[phase], passwd_);
            upwind_->Compute(flux, *kr, op_bcs_[eqnr.first]->bc_model(), *kr);

            pde->Setup(Kptr, kr, Teuchos::null);
            pde->SetBCs(op_bcs_[eqnr.first], op_bcs_[eqnc.first]);
            pde->UpdateMatrices(Teuchos::null, Teuchos::null);
            pde->ApplyBCs(false, false, false);
          }
        }

        // -- advection operator div[ (K f grad du/dv) dv ]
        if ((key = eval_eqns_[row][phase].first) != "") {
          Key der_key = "d" + varp_name[phase] + "_d" + keyc;
          S_->GetFieldEvaluator(varp_name[phase])->HasFieldDerivativeChanged(S_.ptr(), passwd_, keyc);

          if (S_->HasField(der_key)) {
            auto pde = Teuchos::rcp(new Operators::PDE_AdvectionUpwind(adv_list, global_op)); 

            // --- upwind gas molar mobility times molar fraction 
            kr_c = *S_->GetFieldData(key)->ViewComponent("cell");
            kr->ViewComponent("dirichlet_faces")->PutScalar(0.0);  // FIXME
            auto flux = *S_->GetFieldData(flux_name[phase], passwd_);
            upwind_->Compute(flux, *kr, op_bcs_[eqnr.first]->bc_model(), *kr);

            // --- calculate advective flux 
            auto der = S_->GetFieldData(der_key);
            pde_diff_K_->Setup(Kptr, kr, Teuchos::null, rho_l_, gravity_);  // FIXME
            pde_diff_K_->SetBCs(op_bcs_[eqnr.first], op_bcs_[eqnc.first]);
            pde_diff_K_->UpdateMatrices(Teuchos::null, Teuchos::null);
            pde_diff_K_->UpdateFlux(der.ptr(), flux_tmp.ptr());

            // -- populated advection operator
            pde->Setup(*flux_tmp);
            pde->SetBCs(op_bcs_[eqnr.first], op_bcs_[eqnc.first]);
            pde->global_operator()->Init();
            pde->UpdateMatrices(flux_tmp.ptr());
            pde->ApplyBCs(false, false, false);
          }
        }

        // -- advection operator div [ (K df/dv grad p) dv ]
        if ((key = eval_eqns_[row][phase].first) != "") {
          Key der_key = "d" + key + "_d" + keyc;
          S_->GetFieldEvaluator(key)->HasFieldDerivativeChanged(S_.ptr(), passwd_, keyc);

          if (S_->HasField(der_key)) {
            auto pde = Teuchos::rcp(new Operators::PDE_AdvectionUpwind(adv_list, global_op)); 

            // --- upwind derivative
            kr_c = *S_->GetFieldData(der_key)->ViewComponent("cell");
            kr->ViewComponent("dirichlet_faces")->PutScalar(0.0);  // FIXME
            auto flux = *S_->GetFieldData(flux_name[phase], passwd_);
            upwind_->Compute(flux, *kr, op_bcs_[eqnr.first]->bc_model(), *kr);

            // --- calculate advective flux 
            auto var = S_->GetFieldData(varp_name[phase]);
            pde_diff_K_->Setup(Kptr, kr, Teuchos::null, rho_g, gravity_);
            pde_diff_K_->SetBCs(op_bcs_[eqnr.first], op_bcs_[eqnc.first]);
            pde_diff_K_->UpdateMatrices(Teuchos::null, Teuchos::null);
            pde_diff_K_->UpdateFlux(var.ptr(), flux_tmp.ptr());

            // -- populated advection operator
            pde->Setup(*flux_tmp);
            pde->SetBCs(op_bcs_[eqnr.first], op_bcs_[eqnc.first]);
            pde->UpdateMatrices(flux_tmp.ptr());
            pde->ApplyBCs(false, false, false);
          }
        }
      }

      // molecular diffusion
      for (int phase = 0; phase < 2; ++phase) {
        // -- diffusion operator div [ D f du/dv grad dv ]
        if ((key = eval_eqns_[row][2 + phase].first) != "") {
          if (varx_name_[phase] == keyc) {
            der_key = "constant_field";  // DAG does not calculate derivative when u=v
          } else {
            der_key = "d" + varx_name_[phase] + "_d" + keyc;
            S_->GetFieldEvaluator(varx_name_[phase])->HasFieldDerivativeChanged(S_.ptr(), passwd_, keyc);
          }

          if (S_->HasField(der_key)) {
            auto pde = Teuchos::rcp(new Operators::PDE_DiffusionFV(mdf_list, global_op));

            double coef;
            auto D = Teuchos::rcp(new std::vector<WhetStone::Tensor>);
            WhetStone::Tensor Dc(dim_, 1);
            for (int c = 0; c < ncells_owned_; ++c) {
              // if (phase == 0) coef = sat_lc[0][c] * eta_lc[0][c] * mol_diff_l_[eqnr.second];
              if (phase == 0) coef = sat_lc[0][c] * mol_diff_l_[eqnr.second];
              if (phase == 1) coef = (1.0 - sat_lc[0][c]) * eta_gc[0][c] * mol_diff_g_[eqnr.second];
              Dc(0, 0) = phi[0][c] * coef;
              D->push_back(Dc);
            }

            pde->Setup(D, Teuchos::null, Teuchos::null);
            pde->SetBCs(op_bcs_[eqnr.first], op_bcs_[eqnc.first]);
            pde->UpdateMatrices(Teuchos::null, Teuchos::null);
            pde->ApplyBCs(false, false, false);

            double factor = eval_eqns_[row][2 + phase].second;
            if (factor != 1.0) pde->local_op()->Rescale(factor);
          }
        }

        // -- diffusion operator div [ (D df/dv grad du ]
        if ((key = eval_eqns_[row][2 + phase].first) != "" && keyc == saturation_liquid_key_) {
          auto pde = Teuchos::rcp(new Operators::PDE_AdvectionUpwind(adv_list, global_op)); 

          // --- calculate diffusion coefficient
          auto D = Teuchos::rcp(new std::vector<WhetStone::Tensor>);
          WhetStone::Tensor Dc(dim_, 1);
          for (int c = 0; c < ncells_owned_; ++c) {
            if (phase == 0) Dc(0, 0) = phi[0][c] * mol_diff_l_[eqnr.second];
            if (phase == 1) Dc(0, 0) =-eta_gc[0][c] * mol_diff_g_[eqnr.second];
            D->push_back(Dc);
          }

          // --- calculate advective flux 
          auto var = S_->GetFieldData(keyc);
          pde_diff_D_->Setup(D, Teuchos::null, Teuchos::null);
          pde_diff_D_->SetBCs(op_bcs_[eqnr.first], op_bcs_[eqnc.first]);
          pde_diff_D_->UpdateMatrices(Teuchos::null, Teuchos::null);
          pde_diff_D_->UpdateFlux(var.ptr(), flux_tmp.ptr());

          // -- populated advection operator
          pde->Setup(*flux_tmp);
          pde->SetBCs(op_bcs_[eqnr.first], op_bcs_[eqnc.first]);
          pde->UpdateMatrices(flux_tmp.ptr());
          pde->ApplyBCs(false, false, false);

          double factor = eval_eqns_[row][2 + phase].second;
          if (factor != 1.0) pde->local_op()->Rescale(factor);
        }
      }

      // storage term
      if ((key = eval_eqns_[row][4].first) != "") {
        Key der_key = "d" + key + "_d" + keyc;
        S_->GetFieldEvaluator(key)->HasFieldDerivativeChanged(S_.ptr(), passwd_, keyc);

        if (S_->HasField(der_key)) {
          auto pde = Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::CELL, global_op)); 
          auto der = S_->GetFieldData(der_key, key);
          pde->AddAccumulationTerm(*der, dtp, "cell");
        }
      }
    }
  }
      
  // process constraint
  int n = num_primary_ + 1;

  for (int i = 0; i < num_primary_ + 2; ++i) {
    auto eqnc = EquationToSolution(i);
    auto pde = Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::CELL, mesh_)); 
    op_preconditioner_->SetOperatorBlock(n, i, pde->global_operator());
 
    Teuchos::RCP<const Epetra_MultiVector> der_fc, der_gc;

    Key key = eval_eqns_[n][0].first;
    const auto& ncp_fc = *S_->GetFieldData(key)->ViewComponent("cell");

    Key keyc = soln_names_[eqnc.first];
    Key derf_key = "d" + key + "_d" + keyc;
    S_->GetFieldEvaluator(key)->HasFieldDerivativeChanged(S_.ptr(), passwd_, keyc);
    if (S_->HasField(derf_key)) der_fc = S_->GetFieldData(derf_key)->ViewComponent("cell");

    key = eval_eqns_[n][1].first;
    const auto& ncp_gc = *S_->GetFieldData(key)->ViewComponent("cell");

    keyc = soln_names_[eqnc.first];
    Key derg_key = "d" + key + "_d" + keyc;
    S_->GetFieldEvaluator(key)->HasFieldDerivativeChanged(S_.ptr(), passwd_, keyc);
    if (S_->HasField(derg_key)) der_gc = S_->GetFieldData(derg_key)->ViewComponent("cell");

    // -- identify active set for gas phase
    fone.PutScalar(0.0);
    for (int c = 0; c < ncells_owned_; c++) {
      if (ncp_fc[0][c] > ncp_gc[0][c] + 1.0e-12) {
        if (der_gc.get()) fone_c[0][c] = (*der_gc)[0][c];
      } else {
        if (der_fc.get()) fone_c[0][c] = (*der_fc)[0][c];
      }
    }

    pde->AddAccumulationTerm(fone, "cell");
  }


  // finalize preconditioner
  if (!op_pc_assembled_) {
    op_preconditioner_->SymbolicAssembleMatrix();
    op_pc_assembled_ = true;
  }
  op_preconditioner_->AssembleMatrix();
// std::cout << *op_preconditioner_->A() << std::endl; exit(0);
  op_preconditioner_->UpdatePreconditioner();
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
// {double aaa; Y->SubVector(0)->Data()->Norm2(&aaa); std::cout << aaa << std::endl; }
// exit(0);
  return ierr;
}


/* ******************************************************************
* This is called when the time integration scheme changes solution
****************************************************************** */
void Multiphase_PK::ChangedSolution()
{
  for (int i = 0; i < 3; ++i ) {
    auto eval = S_->GetFieldEvaluator(soln_names_[i]);
    Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(eval)->SetFieldAsChanged(S_.ptr());
  }
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
  auto dsc = *du->SubVector(1)->Data()->ViewComponent("cell");

  double error_s = 0.0;
  for (int c = 0; c < ncells_owned_; c++) {
    error_s = std::max(error_s, fabs(dsc[0][c]));
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
      // du2c[i][c] = std::max(du2c[i][c], u2c[i][c] - 1.0);
    }    
  }

  // clip saturation (residual saturation is missing, FIXME)
  const auto& u1c = *u->SubVector(1)->Data()->ViewComponent("cell");
  auto& du1c = *du->SubVector(1)->Data()->ViewComponent("cell");

  for (int c = 0; c < ncells_owned_; ++c) {
    // du1c[0][c] = std::min(du1c[0][c], u1c[0][c]);
    // du1c[0][c] = std::max(du1c[0][c], u1c[0][c] - 1.0);
  }

  return AmanziSolvers::FnBaseDefs::CORRECTION_MODIFIED;
}

}  // namespace Multiphase
}  // namespace Amanzi

