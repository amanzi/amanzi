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
#include "EquationStructure.hh"
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
  // -- gas pressure
  S_->GetEvaluator(pressure_gas_key_).Update(*S_, passwd_);
  auto pg = S_->GetPtrW<CompositeVector>(pressure_gas_key_, pressure_gas_key_);

  // -- molar densities
  S_->GetEvaluator(molar_density_gas_key_).Update(*S_, passwd_);
  S_->GetEvaluator(molar_density_liquid_key_).Update(*S_, passwd_);

  // -- storage
  S_->GetEvaluator(tws_key_).Update(*S_, passwd_);
  S_->GetEvaluator(tcs_key_).Update(*S_, passwd_);

  // -- wrapper for absolute permeability
  Teuchos::RCP<std::vector<WhetStone::Tensor> > Kptr = Teuchos::rcpFromRef(K_);

  // work memory for miscalleneous operator
  auto kr = CreateCVforUpwind(mesh_);
  auto& kr_c = *kr->ViewComponent("cell");
  std::vector<int> bcnone(nfaces_wghost_, Operators::OPERATOR_BC_NONE);

  // primary variables
  auto aux = up[0];
  CompositeVector fone(*aux), fadd(*aux), comp(*aux);
  auto& fone_c = *fone.ViewComponent("cell");
  auto& comp_c = *comp.ViewComponent("cell");

  // start loop over physical equations
  Key key;
  for (int n = 0; n < num_primary_ + 1; ++n) {
    ModifyEvaluators(n);
    auto sol = EquationToSolution(n);
    PopulateBCs(sol.comp, true);
  
    // Richards-type operator for all phases
    fone.PutScalar(0.0);

    for (int phase = 0; phase < 2; ++phase) {
      bool bcflag = (phase == 0);
      if ((key = eqns_[n].advection[phase].first) != "") {
        S_->GetEvaluator(key).Update(*S_, passwd_);

        // -- upwind cell-centered coefficient
        auto& flux = S_->GetW<CompositeVector>(flux_names_[phase], passwd_);
        kr_c = *S_->Get<CompositeVector>(key).ViewComponent("cell");
        upwind_->Compute(flux, *kr, bcnone, *kr);

        // -- form operator
        auto& pde = pde_diff_K_;
        pde->Setup(Kptr, kr, Teuchos::null, rho_l_, gravity_);  // FIXME (gravity for gas phase)
        pde->SetBCs(op_bcs_[sol.var], op_bcs_[sol.var]);
        pde->global_operator()->Init();
        pde->UpdateMatrices(Teuchos::null, Teuchos::null);
        pde->ApplyBCs(bcflag, false, false);

        // -- add advection term to the residual
        Key fname = eqns_[n].advection[phase].second;
        S_->GetEvaluator(fname).Update(*S_, passwd_);
        auto var = S_->GetPtr<CompositeVector>(fname);
        pde->global_operator()->ComputeNegativeResidual(*var, fadd);
        fone.Update(1.0, fadd, 1.0);
      }
    }
 
    // molecular diffusion 
    for (int phase = 0; phase < 2; ++phase) {
      if ((key = eqns_[n].diffusion[phase].first) != "") {
        S_->GetEvaluator(key).Update(*S_, passwd_);
        auto& flux = S_->GetW<CompositeVector>(flux_names_[phase], passwd_);
        kr_c = *S_->Get<CompositeVector>(key).ViewComponent("cell");
        upwind_->Compute(flux, *kr, bcnone, *kr);

        // -- form operator
        auto& pde = pde_diff_D_;
        pde->Setup(Teuchos::null, kr, Teuchos::null);
        pde->SetBCs(op_bcs_[sol.var], op_bcs_[sol.var]);
        pde->global_operator()->Init();
        pde->UpdateMatrices(Teuchos::null, Teuchos::null);
        pde->ApplyBCs(false, false, false);

        // -- add diffusion term to the residual
        Key fname = eqns_[n].diffusion[phase].second;
        S_->GetEvaluator(fname).Update(*S_, passwd_);
        const auto& tmp = *S_->Get<CompositeVector>(fname).ViewComponent("cell");
        int m = std::min(sol.comp, tmp.NumVectors() - 1);
        for (int c = 0; c < ncells_owned_; ++c) {
          comp_c[0][c] = tmp[m][c];
        }
        pde->global_operator()->ComputeNegativeResidual(comp, fadd);

        double factor = eqns_[n].diff_factors[phase];
        fone.Update(factor, fadd, 1.0);
      }
    }

    // add storage terms 
    if ((key = eqns_[n].storage) != "") {
      std::string prev_key = "prev_" + key;
      S_->GetEvaluator(key).Update(*S_, passwd_);

      const auto& total_c = *S_->Get<CompositeVector>(key).ViewComponent("cell");
      const auto& total_prev_c = *S_->Get<CompositeVector>(prev_key).ViewComponent("cell");

      for (int c = 0; c < ncells_owned_; ++c) {
        double factor = mesh_->cell_volume(c) / dtp;
        fone_c[0][c] += (total_c[0][c] - total_prev_c[0][c]) * factor;
      }
    }

    // copy temporary vector to residual
    auto& fc = *fp[sol.var]->ViewComponent("cell");
    for (int c = 0; c < ncells_owned_; ++c)
      fc[sol.comp][c] = fone_c[0][c];
  }

  // process gas constraints
  int n = num_primary_ + 1;
  key = eqns_[n].constraint.first;
  S_->GetEvaluator(key).Update(*S_, passwd_);
  const auto& ncp_fc = *S_->Get<CompositeVector>(key).ViewComponent("cell");

  key = eqns_[n].constraint.second;
  S_->GetEvaluator(key).Update(*S_, passwd_);
  const auto& ncp_gc = *S_->Get<CompositeVector>(key).ViewComponent("cell");

  auto& fci = *fp[2]->ViewComponent("cell");
  if (ncp_ == "min") {
    for (int c = 0; c < ncells_owned_; ++c) {
      fci[0][c] = std::min(ncp_fc[0][c], ncp_gc[0][c]);
    }
  } else if (ncp_ == "Fischer-Burmeister") {
    for (int c = 0; c < ncells_owned_; ++c) {
      double a = ncp_fc[0][c];
      double b = ncp_gc[0][c];
      fci[0][c] = std::pow(a * a + b * b, 0.5) - (a + b);
    }
  }
}


/* ******************************************************************
* Update new preconditioner on the interval (tp-dtp, tp].
****************************************************************** */
void Multiphase_PK::UpdatePreconditioner(double tp, Teuchos::RCP<const TreeVector> u, double dtp)
{
  // extract pointers to subvectors
  std::vector<Teuchos::RCP<const CompositeVector> > up;
  for (int i = 0; i < 3; ++i) {
    up.push_back(u->SubVector(i)->Data());
  }

  // miscalleneous fields
  // -- molar densities
  S_->GetEvaluator(molar_density_gas_key_).Update(*S_, passwd_);
  const auto& eta_g = S_->Get<CompositeVector>(molar_density_gas_key_);

  S_->GetEvaluator(molar_density_liquid_key_).Update(*S_, passwd_);

  // -- mass density of gas phase 
  auto rho_g = Teuchos::rcp(new CompositeVector(eta_g));

  // -- gas pressure
  S_->GetEvaluator(pressure_gas_key_).Update(*S_, passwd_);
  auto pg = S_->GetPtrW<CompositeVector>(pressure_gas_key_, Tags::DEFAULT, pressure_gas_key_);

  // -- wrapper for absolute permeability
  Teuchos::RCP<std::vector<WhetStone::Tensor> > Kptr = Teuchos::rcpFromRef(K_);

  // parameter lists
  auto& adv_list = mp_list_->sublist("operators").sublist("advection operator").sublist("preconditioner");
  auto& ddf_list = mp_list_->sublist("operators").sublist("diffusion operator").sublist("preconditioner");
  auto& mdf_list = mp_list_->sublist("operators").sublist("molecular diffusion operator").sublist("preconditioner");

  // work memory for miscalleneous operator
  auto flux_tmp = Teuchos::rcp(new CompositeVector(S_->Get<CompositeVector>(darcy_flux_liquid_key_)));
  auto flux_acc = Teuchos::rcp(new CompositeVector(*flux_tmp));

  auto kr = CreateCVforUpwind(mesh_);
  auto& kr_c = *kr->ViewComponent("cell");
  std::vector<int> bcnone(nfaces_wghost_, Operators::OPERATOR_BC_NONE);

  auto fone = Teuchos::rcp(new CompositeVector(*up[0]));
  auto fone_c = fone->ViewComponent("cell");
  fone->PutScalar(1.0); 

  const Epetra_MultiVector* der_c;

  // for each operator we linearize (a) functions on which it acts and
  // (b) non-linear coefficients
  Key key;
  for (int row = 0; row < num_primary_ + 1; ++row) {
    ModifyEvaluators(row);
    auto solr = EquationToSolution(row);
    Key keyr = soln_names_[solr.var];
    PopulateBCs(solr.comp, false);

    for (int col = 0; col < num_primary_ + 2; ++col) {
      auto solc = EquationToSolution(col);
      Key keyc = soln_names_[solc.var];

      if (solc.matching_eqn >= 0 && row != solc.matching_eqn) continue;

      bool bcflag = (row == col);

      // add empty operator to have a well-defined global operator pointer
      auto pde0 = Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::CELL, mesh_));
      auto global_op = pde0->global_operator();
      op_preconditioner_->set_operator_block(row, col, global_op);
      kr_c.PutScalar(0.0);
      pde0->AddAccumulationTerm(*kr, "cell");

      // initialize accumulated flux
      flux_acc->PutScalar(0.0);

      //
      // Richards-type operator for all phases, div [K f grad(g)]
      //
      for (int phase = 0; phase < 2; ++phase) {
        // -- diffusion operator div[ (K f dg/dv) grad dv ] 
        if ((key = eqns_[row].advection[phase].first) != "") {
          Key fname = eqns_[row].advection[phase].second;

          if (S_->HasDerivative(fname, keyc) || fname == keyc) {
std::cout << "Try1: " << fname << " " << keyc << " " << row << " " << col << std::endl;
            auto pde = Teuchos::rcp(new Operators::PDE_DiffusionFV(ddf_list, global_op));

            S_->GetEvaluator(key).Update(*S_, passwd_);
            const auto& coef_c = *S_->Get<CompositeVector>(key).ViewComponent("cell");

            if (fname == keyc) {
              der_c = &*fone_c;
            } else {
              S_->GetEvaluator(fname).UpdateDerivative(*S_, passwd_, keyc, Tags::DEFAULT);
              der_c = &*S_->GetDerivative<CompositeVector>(fname, keyc).ViewComponent("cell");
            }

            for (int c = 0; c < ncells_owned_; ++c) {
              kr_c[0][c] = (*der_c)[0][c] * coef_c[0][c];
            }
            auto& flux = S_->GetW<CompositeVector>(flux_names_[phase], passwd_);
            upwind_->Compute(flux, *kr, bcnone, *kr);

            pde->Setup(Kptr, kr, Teuchos::null);
            pde->SetBCs(op_bcs_[solr.var], op_bcs_[solc.var]);
            pde->UpdateMatrices(Teuchos::null, Teuchos::null);
            pde->ApplyBCs(bcflag, false, false);

            double factor = eqns_[row].adv_factors[phase];
            if (factor != 1.0) pde->local_op()->Rescale(factor);
          }
        }

        // -- advection operator div[ (K f grad dg/dv) dv ]
        if ((key = eqns_[row].advection[phase].first) != "") {
          Key fname = eqns_[row].advection[phase].second;

          if (S_->HasDerivative(fname, keyc)) {
std::cout << "Try2: " << fname << " " << keyc << " " << row << " " << col << std::endl;
            // --- upwind gas molar mobility times molar fraction 
            S_->GetEvaluator(key).Update(*S_, passwd_);
            kr_c = *S_->Get<CompositeVector>(key).ViewComponent("cell");
            auto& flux = S_->GetW<CompositeVector>(flux_names_[phase], passwd_);
            upwind_->Compute(flux, *kr, bcnone, *kr);

            // --- calculate advective flux 
            S_->GetEvaluator(fname).UpdateDerivative(*S_, passwd_, keyc, Tags::DEFAULT);
            auto tmp = S_->GetDerivativePtr<CompositeVector>(fname, Tags::DEFAULT, keyc, Tags::DEFAULT);

            pde_diff_K_->Setup(Kptr, kr, Teuchos::null, rho_l_, gravity_);  // FIXME
            pde_diff_K_->SetBCs(op_bcs_[solr.var], op_bcs_[solc.var]);
            pde_diff_K_->UpdateMatrices(Teuchos::null, Teuchos::null);
            pde_diff_K_->UpdateFlux(tmp.ptr(), flux_tmp.ptr());

            double factor = eqns_[row].adv_factors[phase];
            flux_acc->Update(factor, *flux_tmp, 1.0);
          }
        }

        // -- advection operator div [ (K df/dv grad g) dv ]
        if ((key = eqns_[row].advection[phase].first) != "") {

          if (S_->HasDerivative(key, keyc)) {
std::cout << "Try3: " << key << " " << keyc << " " << row << " " << col << std::endl;
            // --- upwind derivative
            S_->GetEvaluator(key).UpdateDerivative(*S_, passwd_, keyc, Tags::DEFAULT);
            kr_c = *S_->GetDerivative<CompositeVector>(key, keyc).ViewComponent("cell");
            const auto& flux = S_->Get<CompositeVector>(flux_names_[phase]);
            upwind_->Compute(flux, *kr, bcnone, *kr);

            // --- calculate advective flux 
            Key fname = eqns_[row].advection[phase].second;
            S_->GetEvaluator(fname).Update(*S_, passwd_);
            auto var = S_->GetPtr<CompositeVector>(fname);
            pde_diff_K_->Setup(Kptr, kr, Teuchos::null, rho_g, gravity_);
            pde_diff_K_->SetBCs(op_bcs_[solr.var], op_bcs_[solc.var]);
            pde_diff_K_->UpdateMatrices(Teuchos::null, Teuchos::null);
            pde_diff_K_->UpdateFlux(var.ptr(), flux_tmp.ptr());

            double factor = eqns_[row].adv_factors[phase];
            flux_acc->Update(factor, *flux_tmp, 1.0);
          }
        }
      }

      //
      // Molecular diffusion
      //
      for (int phase = 0; phase < 2; ++phase) {
        // -- diffusion operator div [ (f dg/dv) grad dv ]
        if ((key = eqns_[row].diffusion[phase].first) != "") {
          Key fname = eqns_[row].diffusion[phase].second;

          if (S_->HasDerivative(fname, keyc) || fname == keyc) {
            auto pde = Teuchos::rcp(new Operators::PDE_DiffusionFV(mdf_list, global_op));

            if (fname == keyc) {
              der_c = &*fone_c;
            } else {
std::cout << "Try4: " << key << " " << keyc << " " << row << " " << col << " fname=" << fname << std::endl;
              S_->GetEvaluator(fname).UpdateDerivative(*S_, passwd_, keyc, Tags::DEFAULT);
              der_c = &*S_->GetDerivative<CompositeVector>(fname, keyc).ViewComponent("cell");
            }

            S_->GetEvaluator(key).Update(*S_, passwd_);
            const auto& coef_c = *S_->Get<CompositeVector>(key).ViewComponent("cell");

            for (int c = 0; c < ncells_owned_; ++c) {
              kr_c[0][c] = (*der_c)[0][c] * coef_c[0][c];
            }
            auto& flux = S_->GetW<CompositeVector>(flux_names_[phase], passwd_);
            upwind_->Compute(flux, *kr, op_bcs_[solr.var]->bc_model(), *kr);

            pde->Setup(Teuchos::null, kr, Teuchos::null);
            pde->SetBCs(op_bcs_[solr.var], op_bcs_[solc.var]);
            pde->UpdateMatrices(Teuchos::null, Teuchos::null);
            pde->ApplyBCs(bcflag, false, false);

            double factor = eqns_[row].diff_factors[phase];
            if (factor != 1.0) pde->local_op()->Rescale(factor);
          }
        }

        // -- advection operator div[ (f grad dg/dv) dv ]
        if ((key = eqns_[row].diffusion[phase].first) != "") {
          Key fname = eqns_[row].diffusion[phase].second;

          if (S_->HasDerivative(fname, keyc)) {
std::cout << "Try5: " << fname << " " << keyc << " " << row << " " << col << std::endl;
            // --- calculate diffusion coefficient
            S_->GetEvaluator(key).Update(*S_, passwd_);
            kr_c = *S_->Get<CompositeVector>(key).ViewComponent("cell");
            auto& flux = S_->GetW<CompositeVector>(flux_names_[phase], passwd_);
            upwind_->Compute(flux, *kr, bcnone, *kr);

            S_->GetEvaluator(fname).UpdateDerivative(*S_, passwd_, keyc, Tags::DEFAULT);
            auto tmp = S_->GetDerivativePtr<CompositeVector>(fname, Tags::DEFAULT, keyc, Tags::DEFAULT);

            // --- calculate advective flux 
            pde_diff_D_->Setup(Teuchos::null, kr, Teuchos::null);
            pde_diff_D_->SetBCs(op_bcs_[solr.var], op_bcs_[solc.var]);
            pde_diff_D_->UpdateMatrices(Teuchos::null, Teuchos::null);
            pde_diff_D_->UpdateFlux(tmp.ptr(), flux_tmp.ptr());

            double factor = eqns_[row].diff_factors[phase];
            flux_acc->Update(factor, *flux_tmp, 1.0);
          }
        }

        // -- advection operator div [ (df/dv grad g) dv ]
        if ((key = eqns_[row].diffusion[phase].first) != "" && keyc == saturation_liquid_key_) {

          if (S_->HasDerivative(key, keyc)) {
std::cout << "Try6: " << key << " " << keyc << " " << row << " " << col << std::endl;
            S_->GetEvaluator(key).UpdateDerivative(*S_, passwd_, keyc, Tags::DEFAULT);
            kr_c = *S_->GetDerivative<CompositeVector>(key, keyc).ViewComponent("cell");

            // --- upwind derivative
            const auto& flux = S_->Get<CompositeVector>(flux_names_[phase]);
            upwind_->Compute(flux, *kr, bcnone, *kr);

            // --- calculate advective flux 
            Key fname = eqns_[row].diffusion[phase].second;
            S_->GetEvaluator(fname).Update(*S_, passwd_);
            auto var = S_->GetPtr<CompositeVector>(fname);
            pde_diff_D_->Setup(Teuchos::null, kr, Teuchos::null);
            pde_diff_D_->SetBCs(op_bcs_[solr.var], op_bcs_[solc.var]);
            pde_diff_D_->UpdateMatrices(Teuchos::null, Teuchos::null);
            pde_diff_D_->UpdateFlux(var.ptr(), flux_tmp.ptr());

            double factor = eqns_[row].diff_factors[phase];
            flux_acc->Update(factor, *flux_tmp, 1.0);
          }
        }
      }

      // populate advection operator
      auto pde1 = Teuchos::rcp(new Operators::PDE_AdvectionUpwind(adv_list, global_op)); 
      pde1->Setup(*flux_acc);
      pde1->SetBCs(op_bcs_[solr.var], op_bcs_[solc.var]);
      pde1->UpdateMatrices(flux_acc.ptr());
      pde1->ApplyBCs(false, false, false);

      // storage term
      if ((key = eqns_[row].storage) != "") {
        if (S_->HasDerivative(key, keyc)) {
std::cout << "Try7: " << key << " " << keyc << " " << row << " " << col << std::endl;
          auto pde = Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::CELL, global_op)); 
          S_->GetEvaluator(key).UpdateDerivative(*S_, passwd_, keyc, Tags::DEFAULT);
          auto der = S_->GetDerivativePtr<CompositeVector>(key, Tags::DEFAULT, keyc, Tags::DEFAULT);
          pde->AddAccumulationTerm(*der, dtp, "cell");
        }
      }
    }
  }
      
  // process constraint
  int n = num_primary_ + 1;

  for (int i = 0; i < num_primary_ + 2; ++i) {
    auto solc = EquationToSolution(i);
    auto pde = Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::CELL, mesh_)); 
    op_preconditioner_->set_operator_block(n, i, pde->global_operator());
 
    // -- derivatives
    Teuchos::RCP<const Epetra_MultiVector> der_fc, der_gc;

    key = eqns_[n].constraint.first;
    S_->GetEvaluator(key).Update(*S_, passwd_);
    const auto& ncp_fc = *S_->Get<CompositeVector>(key).ViewComponent("cell");

    Key keyc = soln_names_[solc.var];
    if (S_->HasDerivative(key, keyc)) {
std::cout << "Try8: " << key << " " << keyc << std::endl;
      S_->GetEvaluator(key).UpdateDerivative(*S_, passwd_, keyc, Tags::DEFAULT);
      der_fc = S_->GetDerivative<CompositeVector>(key, Tags::DEFAULT, keyc, Tags::DEFAULT).ViewComponent("cell");
    }

    key = eqns_[n].constraint.second;
    S_->GetEvaluator(key).Update(*S_, passwd_);
    const auto& ncp_gc = *S_->Get<CompositeVector>(key).ViewComponent("cell");

    keyc = soln_names_[solc.var];
    if (S_->HasDerivative(key, keyc)) {
std::cout << "Try9: " << key << " " << keyc << std::endl;
      S_->GetEvaluator(key).UpdateDerivative(*S_, passwd_, keyc, Tags::DEFAULT);
      der_gc = S_->GetDerivative<CompositeVector>(key, Tags::DEFAULT, keyc, Tags::DEFAULT).ViewComponent("cell");
    }

    // -- identify active set for gas phase
    fone->PutScalar(0.0);
    if (ncp_ == "min") {
      for (int c = 0; c < ncells_owned_; c++) {
        if (ncp_fc[0][c] > ncp_gc[0][c]) {
          if (der_gc.get()) (*fone_c)[0][c] = (*der_gc)[0][c];
        } else {
          if (der_fc.get()) (*fone_c)[0][c] = (*der_fc)[0][c];
        }
      }
    } else if (ncp_ == "Fischer-Burmeister") {
      for (int c = 0; c < ncells_owned_; ++c) {
        double a = ncp_fc[0][c];
        double b = ncp_gc[0][c];

        double da = (der_fc.get()) ? (*der_fc)[0][c] : 0.0;
        double db = (der_gc.get()) ? (*der_gc)[0][c] : 0.0;
        (*fone_c)[0][c] = (a * da + b * db) * std::pow(a * a + b * b, -0.5) - (da + db);
      }
    }

    pde->AddAccumulationTerm(*fone, "cell");
  }

  // finalize preconditioner
  if (!op_pc_assembled_) {
    op_pc_solver_->InitializeInverse();
    op_pc_assembled_ = true;
  }
  op_pc_solver_->ComputeInverse();
}


/* ******************************************************************
* Apply preconditioner inv(B) * X.                                                 
****************************************************************** */
int Multiphase_PK::ApplyPreconditioner(Teuchos::RCP<const TreeVector> X, 
                                       Teuchos::RCP<TreeVector> Y)
{
  Y->PutScalar(0.0);
  // *Y = *X; return 0;
  // return op_preconditioner_->ApplyInverse(*X, *Y);
  return op_pc_solver_->ApplyInverse(*X, *Y);
}


/* ******************************************************************
* This is called when the time integration scheme changes solution
****************************************************************** */
void Multiphase_PK::ChangedSolution()
{
  for (int i = 0; i < 3; ++i ) {
    Teuchos::rcp_dynamic_cast<EvaluatorPrimary<CompositeVector, CompositeVectorSpace> >(
        S_->GetEvaluatorPtr(soln_names_[i]))->SetChanged();
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

  // concentration error
  auto dxc = *du->SubVector(1)->Data()->ViewComponent("cell");
  int n = dxc.NumVectors();

  double error_x = 0.0;
  for (int c = 0; c < ncells_owned_; c++) {
    for (int i = 0; i < n; ++i) {
      error_x = std::max(error_x, fabs(dxc[i][c]));
    }
  }

  // saturation error
  auto dsc = *du->SubVector(2)->Data()->ViewComponent("cell");

  double error_s = 0.0;
  for (int c = 0; c < ncells_owned_; c++) {
    error_s = std::max(error_s, fabs(dsc[0][c]));
  }

  double error = error_p + error_x + error_s;
#ifdef HAVE_MPI
  double buf = error;
  du->Comm()->MaxAll(&buf, &error, 1);  // find the global maximum
#endif

  return error;
}

}  // namespace Multiphase
}  // namespace Amanzi

