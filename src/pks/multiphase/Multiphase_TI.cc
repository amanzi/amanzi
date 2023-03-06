/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Quan Bui (mquanbui@math.umd.edu)
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  MultiPhase PK

  Multiphase multi-component flow, see Multiphase_PK.cc for more detail.
*/


// TPLs
#include "Teuchos_RCP.hpp"

// Amanzi
#include "PDE_Accumulation.hh"
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
void
Multiphase_PK::FunctionalResidual(double t_old,
                                  double t_new,
                                  Teuchos::RCP<TreeVector> u_old,
                                  Teuchos::RCP<TreeVector> u_new,
                                  Teuchos::RCP<TreeVector> f)
{
  double dtp = t_new - t_old;

  // extract pointers to subvectors
  std::vector<Teuchos::RCP<CompositeVector>> up, fp;
  for (int i = 0; i < soln_names_.size(); ++i) {
    up.push_back(u_new->SubVector(i)->Data());
    fp.push_back(f->SubVector(i)->Data());
  }

  // miscalleneous fields
  // -- gas pressure
  S_->GetEvaluator(pressure_gas_key_).Update(*S_, passwd_);
  auto pg = S_->GetPtrW<CompositeVector>(pressure_gas_key_, Tags::DEFAULT, pressure_gas_key_);

  // -- mass and molar densities
  S_->GetEvaluator(mass_density_gas_key_).Update(*S_, passwd_);
  S_->GetEvaluator(mol_density_gas_key_).Update(*S_, passwd_);
  S_->GetEvaluator(mol_density_liquid_key_).Update(*S_, passwd_);

  // -- storage
  S_->GetEvaluator(tws_key_).Update(*S_, passwd_);
  S_->GetEvaluator(tcs_key_).Update(*S_, passwd_);

  // -- wrapper for absolute permeability
  Teuchos::RCP<std::vector<WhetStone::Tensor>> Kptr = Teuchos::rcpFromRef(K_);

  // work memory for miscalleneous operator
  auto kr = CreateCVforUpwind_();
  auto& kr_c = *kr->ViewComponent("cell");

  std::vector<int> bcnone(nfaces_wghost_, Operators::OPERATOR_BC_NONE);

  // primary variables
  auto aux = up[0];
  CompositeVector fone(*aux), fadd(*aux), comp(*aux);
  auto& fone_c = *fone.ViewComponent("cell");
  auto& comp_c = *comp.ViewComponent("cell");

  // start loop over physical equations
  int neqns = eqns_.size();
  int nncps = (system_["constraint eqn"]) ? 1 : 0;

  Key fname, gname, key;
  for (int n = 0; n < neqns - nncps; ++n) {
    Key keyr = soln_names_[eqns_flattened_[n][0]];

    ModifyEvaluators(n);
    PopulateBCs(eqns_flattened_[n][1], true);

    // Richards-type operator for all phases, div [f K grad g]
    fone.PutScalar(0.0);

    for (int phase = 0; phase < 2; ++phase) {
      fname = eqns_[n].advection[phase].first;
      gname = eqns_[n].advection[phase].second;

      if (fname != "") {
        CheckCompatibilityBCs(keyr, gname);
        S_->GetEvaluator(fname).Update(*S_, passwd_);

        // -- upwind cell-centered coefficient
        auto& flux = S_->GetW<CompositeVector>(flux_names_[phase], passwd_);
        kr_c = *S_->Get<CompositeVector>(fname).ViewComponent("cell");
        upwind_->Compute(flux, bcnone, *kr);

        // -- form diffusion operator for variable g
        //    Neuman BCs: separate fluxes for each phase OR the total flux but only once
        //    Dirichlet BC: ellimination could be done independently over phases
        auto pde = pde_diff_K_[phase];
        pde->Setup(Kptr, kr, Teuchos::null); // FIXME (gravity for gas phase)
        pde->SetBCs(op_bcs_[gname], op_bcs_[gname]);
        pde->global_operator()->Init();
        pde->UpdateMatrices(Teuchos::null, Teuchos::null);
        pde->ApplyBCs(true, false, false);

        // -- add advection term to the residual
        S_->GetEvaluator(gname).Update(*S_, passwd_);
        auto var = S_->GetPtr<CompositeVector>(gname, Tags::DEFAULT);
        pde->global_operator()->ComputeNegativeResidual(*var, fadd);

        double factor = eqns_[n].adv_factors[phase];
        fone.Update(factor, fadd, 1.0);
      }
    }

    // Molecular diffusion for all phases, div [f M grad g]
    for (int phase = 0; phase < 2; ++phase) {
      fname = eqns_[n].diffusion[phase].first;
      gname = eqns_[n].diffusion[phase].second;

      if (fname != "") {
        S_->GetEvaluator(fname).Update(*S_, passwd_);
        auto& flux = S_->GetW<CompositeVector>(flux_names_[phase], passwd_);
        kr_c = *S_->Get<CompositeVector>(fname).ViewComponent("cell");
        upwind_->Compute(flux, bcnone, *kr);

        // -- form diffusion operator
        AMANZI_ASSERT(op_bcs_.find(gname) != op_bcs_.end());

        auto pde = pde_diff_D_;
        pde->Setup(Teuchos::null, kr, Teuchos::null);
        pde->SetBCs(op_bcs_[gname], op_bcs_[gname]);
        pde->global_operator()->Init();
        pde->UpdateMatrices(Teuchos::null, Teuchos::null);
        pde->ApplyBCs(true, false, false);

        // -- add diffusion term to the residual
        S_->GetEvaluator(gname).Update(*S_, passwd_);
        const auto& tmp = *S_->Get<CompositeVector>(gname).ViewComponent("cell");
        int m = std::min(eqns_flattened_[n][1], tmp.NumVectors() - 1);
        for (int c = 0; c < ncells_owned_; ++c) { comp_c[0][c] = tmp[m][c]; }
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
    auto& fc = *fp[eqns_flattened_[n][0]]->ViewComponent("cell");
    for (int c = 0; c < ncells_owned_; ++c) fc[eqns_flattened_[n][1]][c] = fone_c[0][c];
  }

  // process gas constraints
  if (nncps > 0) {
    int n = neqns - 1;
    key = eqns_[n].constraint.first;
    S_->GetEvaluator(key).Update(*S_, passwd_);
    const auto& ncp_fc = *S_->Get<CompositeVector>(key).ViewComponent("cell");

    key = eqns_[n].constraint.second;
    S_->GetEvaluator(key).Update(*S_, passwd_);
    const auto& ncp_gc = *S_->Get<CompositeVector>(key).ViewComponent("cell");

    auto& fci = *fp[eqns_flattened_[n][0]]->ViewComponent("cell");
    if (ncp_ == "min") {
      for (int c = 0; c < ncells_owned_; ++c) { fci[0][c] = std::min(ncp_fc[0][c], ncp_gc[0][c]); }
    } else if (ncp_ == "Fischer-Burmeister") {
      for (int c = 0; c < ncells_owned_; ++c) {
        double a = ncp_fc[0][c];
        double b = ncp_gc[0][c];
        fci[0][c] = std::pow(a * a + b * b, 0.5) - (a + b);
      }
    }
  }
}


/* ******************************************************************
* Update new preconditioner on the interval (tp-dtp, tp].
****************************************************************** */
void
Multiphase_PK::UpdatePreconditioner(double tp, Teuchos::RCP<const TreeVector> u, double dtp)
{
  // extract pointers to subvectors
  std::vector<Teuchos::RCP<const CompositeVector>> up;
  for (int i = 0; i < soln_names_.size(); ++i) { up.push_back(u->SubVector(i)->Data()); }

  // miscalleneous fields
  // -- molar densities
  S_->GetEvaluator(mol_density_gas_key_).Update(*S_, passwd_);
  S_->GetEvaluator(mol_density_liquid_key_).Update(*S_, passwd_);

  // -- gas pressure
  S_->GetEvaluator(pressure_gas_key_).Update(*S_, passwd_);
  auto pg = S_->GetPtrW<CompositeVector>(pressure_gas_key_, Tags::DEFAULT, pressure_gas_key_);

  // -- wrapper for absolute permeability
  Teuchos::RCP<std::vector<WhetStone::Tensor>> Kptr = Teuchos::rcpFromRef(K_);

  // work memory for miscalleneous operator
  auto flux_tmp =
    Teuchos::rcp(new CompositeVector(S_->Get<CompositeVector>(vol_flowrate_liquid_key_)));
  auto flux_acc = Teuchos::rcp(new CompositeVector(*flux_tmp));

  auto kr = CreateCVforUpwind_();
  auto& kr_c = *kr->ViewComponent("cell");

  std::vector<int> bcnone(nfaces_wghost_, Operators::OPERATOR_BC_NONE);

  auto fone = Teuchos::rcp(new CompositeVector(*up[0]));
  auto fone_c = fone->ViewComponent("cell");
  fone->PutScalar(1.0);

  const Epetra_MultiVector* der_c;

  // for each operator we linearize (a) functions on which it acts and
  // (b) non-linear coefficients
  int neqns = eqns_.size();
  int nncps = (system_["constraint eqn"]) ? 1 : 0;

  Key fname, gname, key;
  for (int row = 0; row < neqns - nncps; ++row) {
    ModifyEvaluators(row);
    Key keyr = soln_names_[eqns_flattened_[row][0]];
    PopulateBCs(eqns_flattened_[row][1], false);

    for (int col = 0; col < neqns; ++col) {
      Key keyc = soln_names_[eqns_flattened_[col][0]];

      if (eqns_flattened_[col][2] >= 0 && row != eqns_flattened_[col][2]) continue;

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
      // The pair of evaluators or fields describing this equation is (f, g).
      //
      for (int phase = 0; phase < 2; ++phase) {
        fname = eqns_[row].advection[phase].first;
        gname = eqns_[row].advection[phase].second;

        // -- diffusion operator div[ K f grad((dg/dv) dv) ]
        //    BCs are defined by the equation and must be imposed only once
        //    using either the total flux or Dirichlet
        if (fname != "") {
          if (S_->HasDerivative(gname, keyc) || gname == keyc) {
            auto pde = fac_diffK_->Create(global_op);

            S_->GetEvaluator(fname).Update(*S_, passwd_);
            kr_c = *S_->Get<CompositeVector>(fname).ViewComponent("cell");

            auto& flux = S_->GetW<CompositeVector>(flux_names_[phase], passwd_);
            upwind_->Compute(flux, bcnone, *kr);

            pde->Setup(Kptr, kr, Teuchos::null);
            pde->SetBCs(op_bcs_[keyc], op_bcs_[keyc]);
            pde->UpdateMatrices(Teuchos::null, Teuchos::null);

            if (gname != keyc) {
              S_->GetEvaluator(gname).UpdateDerivative(*S_, passwd_, keyc, Tags::DEFAULT);
              const auto& der_c = S_->GetDerivative<CompositeVector>(
                gname, Tags::DEFAULT, keyc, Tags::DEFAULT);

              pde->ScaleMatricesColumns(der_c);
            }

            pde->ApplyBCs(false, false, false);

            double factor = eqns_[row].adv_factors[phase];
            if (factor != 1.0) pde->local_op()->Rescale(factor);
          }
        }

        // -- advection operator div [ (K df/dv grad g) dv ]
        if (fname != "") {
          if (S_->HasDerivative(fname, keyc)) {
            // --- upwind derivative
            S_->GetEvaluator(fname).UpdateDerivative(*S_, passwd_, keyc, Tags::DEFAULT);
            kr_c = *S_->GetDerivative<CompositeVector>(fname, Tags::DEFAULT, keyc, Tags::DEFAULT)
                      .ViewComponent("cell");
            const auto& flux = S_->Get<CompositeVector>(flux_names_[phase]);
            upwind_->Compute(flux, bcnone, *kr);

            // --- calculate advective flux
            S_->GetEvaluator(gname).Update(*S_, passwd_);
            auto var = S_->GetPtr<CompositeVector>(gname, Tags::DEFAULT);

            auto pde = pde_diff_K_[phase];
            pde->Setup(Kptr, kr, Teuchos::null);
            pde->SetBCs(op_bcs_[gname], op_bcs_[gname]);
            pde->UpdateMatrices(Teuchos::null, Teuchos::null);
            pde->UpdateFlux(var.ptr(), flux_tmp.ptr());

            double factor = eqns_[row].adv_factors[phase];
            flux_acc->Update(factor, *flux_tmp, 1.0);
          }
        }
      }

      //
      // Molecular diffusion
      //
      for (int phase = 0; phase < 2; ++phase) {
        fname = eqns_[row].diffusion[phase].first;
        gname = eqns_[row].diffusion[phase].second;

        // -- diffusion operator div [ f grad((dg/dv) dv) ]
        if (fname != "") {
          if (S_->HasDerivative(gname, keyc) || gname == keyc) {
            auto pde = fac_diffD_->Create(global_op);

            S_->GetEvaluator(fname).Update(*S_, passwd_);
            kr_c = *S_->Get<CompositeVector>(fname).ViewComponent("cell");

            auto& flux = S_->GetW<CompositeVector>(flux_names_[phase], passwd_);
            upwind_->Compute(flux, bcnone, *kr);

            pde->Setup(Teuchos::null, kr, Teuchos::null);
            pde->SetBCs(op_bcs_[keyc], op_bcs_[keyc]);
            pde->UpdateMatrices(Teuchos::null, Teuchos::null);

            if (gname != keyc) {
              S_->GetEvaluator(gname).UpdateDerivative(*S_, passwd_, keyc, Tags::DEFAULT);
              const auto& der_c = S_->GetDerivative<CompositeVector>(
                gname, Tags::DEFAULT, keyc, Tags::DEFAULT);

              pde->ScaleMatricesColumns(der_c);
            }

            pde->ApplyBCs(false, false, false);

            double factor = eqns_[row].diff_factors[phase];
            if (factor != 1.0) pde->local_op()->Rescale(factor);
          }
        }

        // -- advection operator div [ (df/dv grad g) dv ]
        if (fname != "" && keyc == saturation_liquid_key_) {
          if (S_->HasDerivative(fname, keyc)) {
            S_->GetEvaluator(fname).UpdateDerivative(*S_, passwd_, keyc, Tags::DEFAULT);
            kr_c = *S_->GetDerivative<CompositeVector>(fname, Tags::DEFAULT, keyc, Tags::DEFAULT)
                      .ViewComponent("cell");

            // --- upwind derivative
            const auto& flux = S_->Get<CompositeVector>(flux_names_[phase]);
            upwind_->Compute(flux, bcnone, *kr);

            // --- calculate advective flux
            S_->GetEvaluator(gname).Update(*S_, passwd_);
            auto var = S_->GetPtr<CompositeVector>(gname, Tags::DEFAULT);
            pde_diff_D_->Setup(Teuchos::null, kr, Teuchos::null);
            pde_diff_D_->SetBCs(op_bcs_[gname], op_bcs_[gname]);
            pde_diff_D_->UpdateMatrices(Teuchos::null, Teuchos::null);
            pde_diff_D_->UpdateFlux(var.ptr(), flux_tmp.ptr());

            double factor = eqns_[row].diff_factors[phase];
            flux_acc->Update(factor, *flux_tmp, 1.0);
          }
        }
      }

      // populate advection operator
      auto pde1 = fac_adv_->Create(global_op);
      pde1->Setup(*flux_acc);
      pde1->SetBCs(op_bcs_[keyr], op_bcs_[keyr]);
      pde1->UpdateMatrices(flux_acc.ptr());
      pde1->ApplyBCs(false, false, false);

      // storage term
      if ((key = eqns_[row].storage) != "") {
        if (S_->HasDerivative(key, keyc)) {
          auto pde = Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::CELL, global_op));
          S_->GetEvaluator(key).UpdateDerivative(*S_, passwd_, keyc, Tags::DEFAULT);
          auto der = S_->GetDerivativePtr<CompositeVector>(key, Tags::DEFAULT, keyc, Tags::DEFAULT);
          pde->AddAccumulationTerm(*der, dtp, "cell");
        }
      }
    }
  }

  // process constraint
  if (nncps > 0) {
    int n = neqns - 1;

    for (int i = 0; i < neqns; ++i) {
      auto pde = Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::CELL, mesh_));
      op_preconditioner_->set_operator_block(n, i, pde->global_operator());

      // -- derivatives
      Teuchos::RCP<const Epetra_MultiVector> der_fc, der_gc;

      key = eqns_[n].constraint.first;
      S_->GetEvaluator(key).Update(*S_, passwd_);
      const auto& ncp_fc = *S_->Get<CompositeVector>(key).ViewComponent("cell");

      Key keyc = soln_names_[eqns_flattened_[i][0]];
      if (S_->HasDerivative(key, keyc)) {
        S_->GetEvaluator(key).UpdateDerivative(*S_, passwd_, keyc, Tags::DEFAULT);
        der_fc = S_->GetDerivative<CompositeVector>(key, Tags::DEFAULT, keyc, Tags::DEFAULT)
                   .ViewComponent("cell");
      }

      key = eqns_[n].constraint.second;
      S_->GetEvaluator(key).Update(*S_, passwd_);
      const auto& ncp_gc = *S_->Get<CompositeVector>(key).ViewComponent("cell");

      keyc = soln_names_[eqns_flattened_[i][0]];
      if (S_->HasDerivative(key, keyc)) {
        S_->GetEvaluator(key).UpdateDerivative(*S_, passwd_, keyc, Tags::DEFAULT);
        der_gc = S_->GetDerivative<CompositeVector>(key, Tags::DEFAULT, keyc, Tags::DEFAULT)
                   .ViewComponent("cell");
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
  }
// op_preconditioner_->SymbolicAssembleMatrix();
// op_preconditioner_->AssembleMatrix();
// auto J = FiniteDifferenceJacobian_(tp - dtp, tp, u, u, 1e-6);
// std::cout << J << std::endl;
// std::cout << *op_preconditioner_->A() << std::endl; exit(0); 

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
int
Multiphase_PK::ApplyPreconditioner(Teuchos::RCP<const TreeVector> X, Teuchos::RCP<TreeVector> Y)
{
  Y->PutScalar(0.0);
  // return op_preconditioner_->ApplyInverse(*X, *Y);
  int ok = op_pc_solver_->ApplyInverse(*X, *Y);
  return ok;
}


/* ******************************************************************
* This is called when the time integration scheme changes solution
****************************************************************** */
void
Multiphase_PK::ChangedSolution()
{
  for (int i = 0; i < soln_names_.size(); ++i) {
    Teuchos::rcp_dynamic_cast<EvaluatorPrimary<CompositeVector, CompositeVectorSpace>>(
      S_->GetEvaluatorPtr(soln_names_[i], Tags::DEFAULT))
      ->SetChanged();
  }
}


/* ******************************************************************
* Monitor l2 norm of residual
****************************************************************** */
double
Multiphase_PK::ErrorNorm(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<const TreeVector> du)
{
  double error(0.0);
  for (int i = 0; i < soln_names_.size(); ++i) {
    double error_tmp(0.0);

    if (soln_names_[i] == pressure_liquid_key_) {
      auto& pc = *u->SubVector(i)->Data()->ViewComponent("cell");
      auto& dpc = *du->SubVector(i)->Data()->ViewComponent("cell");

      double atm_pressure(1.0e+5);
      for (int c = 0; c < ncells_owned_; c++) {
        double tmp = fabs(dpc[0][c]) / (fabs(pc[0][c] - atm_pressure) + atm_pressure);
        error_tmp = std::max(error_tmp, tmp);
      }
    }

    else if (soln_names_[i] == tcc_liquid_key_ || soln_names_[i] == tcc_gas_key_ ||
             soln_names_[i] == x_gas_key_) {
      auto& xc = *u->SubVector(i)->Data()->ViewComponent("cell");
      auto& dxc = *du->SubVector(i)->Data()->ViewComponent("cell");
      int n = dxc.NumVectors();

      double floor(1e-10);
      for (int c = 0; c < ncells_owned_; c++) {
        for (int k = 0; k < n; ++k) {
          error_tmp = std::max(error_tmp, fabs(dxc[k][c]) / (xc[k][c] + floor));
        }
      }
    }

    else if (soln_names_[i] == saturation_liquid_key_) {
      auto& dsc = *du->SubVector(i)->Data()->ViewComponent("cell");

      for (int c = 0; c < ncells_owned_; c++) { error_tmp = std::max(error_tmp, fabs(dsc[0][c])); }
    }
    error += error_tmp;
  }

#ifdef HAVE_MPI
  double buf = error;
  du->Comm()->MaxAll(&buf, &error, 1); // find the global maximum
#endif

  return error;
}


/* ******************************************************************
* Debug tools
****************************************************************** */
WhetStone::DenseMatrix
Multiphase_PK::FiniteDifferenceJacobian_(double t_old,
                                         double t_new,
                                         Teuchos::RCP<const TreeVector> u_old,
                                         Teuchos::RCP<const TreeVector> u_new,
                                         double eps)
{
  auto f0 = Teuchos::rcp(new TreeVector(*u_old));
  auto f1 = Teuchos::rcp(new TreeVector(*u_old));
  auto u0 = Teuchos::rcp(new TreeVector(*u_old));
  auto u1 = Teuchos::rcp_const_cast<TreeVector>(u_new);

  int nJ = 3 * ncells_owned_;
  WhetStone::DenseMatrix J(nJ, nJ);

  J.PutScalar(0.0);
  FunctionalResidual(t_old, t_new, u0, u1, f0);

  for (int ncol = 0; ncol < nJ; ++ncol) {
  // for (int ncol = 102; ncol < 103; ++ncol) {
    int n = ncol / ncells_owned_;
    int c = ncol % ncells_owned_;

    ChangedSolution();

    auto& u1_c = *u1->SubVector(n)->Data()->ViewComponent("cell");
    double factor = eps * u1_c[0][c];
    if (n == 2) factor = -eps;
    if (factor == 0.0) continue;

    u1_c[0][c] += factor;
    FunctionalResidual(t_old, t_new, u0, u1, f1);
    f1->Update(-1.0 / factor, *f0, 1.0 / factor);
    u1_c[0][c] -= factor;

    for (int nrow = 0; nrow < nJ; ++nrow) {
      int m = nrow / ncells_owned_;
      int i = nrow % ncells_owned_;

      auto& f1_c = *f1->SubVector(m)->Data()->ViewComponent("cell");
      J(nrow, ncol) = f1_c[0][i];
    }
  }

  return J;
}

} // namespace Multiphase
} // namespace Amanzi
