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
                                  Teuchos::RCP<const TreeVector> u_old,
                                  Teuchos::RCP<TreeVector> u_new,
                                  Teuchos::RCP<TreeVector> f)
{
  double dtp = t_new - t_old;
  // update to handle time dependent BCs
  for (int i = 0; i < bcs_.size(); ++i) {
    bcs_[i]->Compute(t_new, t_new);
  }
  // update source terms
  for (int i = 0; i < srcs_.size(); ++i) {
    srcs_[i]->Compute(t_new, t_new);
  }
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
    PopulateBCs(eqns_[n].component_id, true);

    // Richards-type "advection" and molecular diffusion operators of
    // the form  div [f K grad g]  where f and g are evaluators
    fone.PutScalar(0.0);

    int nterms = eqns_[n].terms.size();
    for (int i = 0; i < nterms; ++i) {
      fname = eqns_[n].terms[i].first;
      gname = eqns_[n].terms[i].second;

      int type = eqns_[n].types[i];
      int phase = eqns_[n].phases[i];

      if (type == 0) CheckCompatibilityBCs(keyr, gname);
      S_->GetEvaluator(fname).Update(*S_, passwd_);
      S_->GetEvaluator(gname).Update(*S_, passwd_);

      // -- upwind cell-centered coefficient
      auto var = S_->GetPtr<CompositeVector>(gname, Tags::DEFAULT);
      auto flux = S_->GetPtrW<CompositeVector>(flux_names_[phase], Tags::DEFAULT, passwd_);
      kr_c = *S_->Get<CompositeVector>(fname).ViewComponent("cell");
      upwind_->Compute(*flux, bcnone, *kr);

      // -- form diffusion operator for variable g
      //    Neuman BCs: separate fluxes for each phase OR the total flux but only once
      //    Dirichlet BC: ellimination could be done independently over phases
      AMANZI_ASSERT(op_bcs_.find(gname) != op_bcs_.end());

      Teuchos::RCP<Operators::PDE_Diffusion> pde;
      pde = (type == 0) ? pde_diff_K_[phase] : pde_diff_D_;

      pde->SetScalarCoefficient(kr, Teuchos::null); // FIXME (gravity for gas phase)
      pde->SetBCs(op_bcs_[gname], op_bcs_[gname]);
      pde->global_operator()->Init();
      pde->UpdateMatrices(flux.ptr(), var.ptr());
      pde->ApplyBCs(true, false, false);

      // -- add advection term to the residual
      if (var->ViewComponent("cell")->NumVectors() == 1) {
        pde->global_operator()->ComputeNegativeResidual(*var, fadd);
      } else {
        const auto& tmp = *var->ViewComponent("cell");
        int m = std::min(eqns_flattened_[n][1], tmp.NumVectors() - 1);
        for (int c = 0; c < ncells_owned_; ++c) {
          comp_c[0][c] = tmp[m][c];
        }
        pde->global_operator()->ComputeNegativeResidual(comp, fadd);
      }

      double factor = eqns_[n].factors[i];
      fone.Update(factor, fadd, 1.0);
    }

    // add storage terms
    if ((key = eqns_[n].storage) != "") {
      std::string prev_key = "prev_" + key;
      S_->GetEvaluator(key).Update(*S_, passwd_);

      const auto& total_c = *S_->Get<CompositeVector>(key).ViewComponent("cell");
      const auto& total_prev_c = *S_->Get<CompositeVector>(prev_key).ViewComponent("cell");

      for (int c = 0; c < ncells_owned_; ++c) {
        double factor = mesh_->getCellVolume(c) / dtp;
        fone_c[0][c] += (total_c[0][c] - total_prev_c[0][c]) * factor;
      }
    }

    // add source term
    if (srcs_.size() > 0) {
      int id = eqns_[n].component_id;
      for (int j = 0; j < srcs_.size(); ++j) {
        if (srcs_[j]->component_id() == id) {
          for (auto it = srcs_[j]->begin(); it != srcs_[j]->end(); ++it) {
            int c = it->first;
            double factor = mesh_->getCellVolume(c);
            fone_c[0][c] -= it->second[0] * factor;
          }
        }
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
}


/* ******************************************************************
* Update new preconditioner on the interval (tp-dtp, tp].
****************************************************************** */
void
Multiphase_PK::UpdatePreconditioner(double tp, Teuchos::RCP<const TreeVector> u, double dtp)
{
  op_preconditioner_->Init();

  // extract pointers to subvectors
  std::vector<Teuchos::RCP<const CompositeVector>> up;
  for (int i = 0; i < soln_names_.size(); ++i) {
    up.push_back(u->SubVector(i)->Data());
  }

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
    Teuchos::rcp(new CompositeVector(S_->Get<CompositeVector>(mol_flowrate_liquid_key_)));
  auto flux_acc = Teuchos::rcp(new CompositeVector(*flux_tmp));

  auto kr = CreateCVforUpwind_();
  auto& kr_c = *kr->ViewComponent("cell");

  std::vector<int> bcnone(nfaces_wghost_, Operators::OPERATOR_BC_NONE);

  auto fone = Teuchos::rcp(new CompositeVector(*up[0]));
  auto fone_c = fone->ViewComponent("cell");
  fone->PutScalar(1.0);

  // update fluxes
  /*
  for (int phase = 0; phase < 2; ++phase) {
    S_->GetEvaluator(adv_names_[phase]).Update(*S_, passwd_);

    auto var = S_->GetPtr<CompositeVector>(pressure_names_[phase], Tags::DEFAULT);
    auto flux = S_->GetPtrW<CompositeVector>(flux_names_[phase], Tags::DEFAULT, passwd_);

    kr_c = *S_->Get<CompositeVector>(adv_names_[phase]).ViewComponent("cell");
    upwind_->Compute(*flux, bcnone, *kr);

    auto pdeK = pde_diff_K_[phase];
    pdeK->SetScalarCoefficient(kr, Teuchos::null);
    pdeK->SetBCs(op_bcs_[pressure_names_[phase]], op_bcs_[pressure_names_[phase]]);
    pdeK->global_operator()->Init();
    pdeK->UpdateMatrices(flux.ptr(), var.ptr());
    pdeK->UpdateFlux(var.ptr(), flux.ptr());
  }
  */

  // for each operator we linearize (a) functions on which it acts and
  // (b) non-linear coefficients
  int neqns = eqns_.size();
  int nncps = (system_["constraint eqn"]) ? 1 : 0;

  Key fname, gname, key;
  for (int row = 0; row < neqns - nncps; ++row) {
    ModifyEvaluators(row);
    Key keyr = soln_names_[eqns_flattened_[row][0]];
    PopulateBCs(eqns_[row].component_id, false);

    for (int col = 0; col < neqns; ++col) {
      Key keyc = soln_names_[eqns_flattened_[col][0]];

      if (eqns_flattened_[col][2] >= 0 && row != eqns_flattened_[col][2]) continue;

      // add empty operator to have a well-defined global operator pointer
      auto pde0 =
        Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::Entity_kind::CELL, mesh_));
      auto global_op = pde0->global_operator();
      op_preconditioner_->set_operator_block(row, col, global_op);
      kr_c.PutScalar(0.0);
      pde0->AddAccumulationTerm(*kr, "cell");

      // --------------------------------------------------------------------
      // Richards-type operator for all phases, div [K f grad(g)]
      // The pair of evaluators or fields describing this equation is (f, g).
      // --------------------------------------------------------------------
      int nterms = eqns_[row].terms.size();
      for (int i = 0; i < nterms; ++i) {
        int type = eqns_[row].types[i];
        int phase = eqns_[row].phases[i];

        fname = eqns_[row].terms[i].first;
        gname = eqns_[row].terms[i].second;

        // initialize accumulated flux
        flux_acc->PutScalar(0.0);

        // -- diffusion operator div{ K f grad[(dg/dv) dv] }
        //    BCs are defined by the equation and must be imposed only once
        //    using either the total flux or Dirichlet
        if (S_->HasDerivative(gname, keyc) || gname == keyc) {
          Teuchos::RCP<Operators::PDE_Diffusion> pde;
          pde = (type == 0) ? fac_diffK_->Create(global_op) : fac_diffD_->Create(global_op);

          S_->GetEvaluator(fname).Update(*S_, passwd_);
          kr_c = *S_->Get<CompositeVector>(fname).ViewComponent("cell");

          auto var = S_->GetPtr<CompositeVector>(gname, Tags::DEFAULT);
          auto flux = S_->GetPtr<CompositeVector>(flux_names_[phase], Tags::DEFAULT);
          upwind_->Compute(*flux, bcnone, *kr);

          if (type == 0) pde->Setup(Kptr, kr, Teuchos::null);
          else pde->Setup(Teuchos::null, kr, Teuchos::null);
          pde->SetBCs(op_bcs_[keyc], op_bcs_[keyc]);
          pde->UpdateMatrices(flux.ptr(), var.ptr());

          if (gname != keyc) {
            S_->GetEvaluator(gname).UpdateDerivative(*S_, passwd_, keyc, Tags::DEFAULT);
            const auto& der_c =
              S_->GetDerivative<CompositeVector>(gname, Tags::DEFAULT, keyc, Tags::DEFAULT);
            pde->ScaleMatricesColumns(der_c);
          }

          pde->ApplyBCs(false, false, false);

          double factor = eqns_[row].factors[i];
          if (factor != 1.0) pde->local_op()->Rescale(factor);
        }

        // -- advection operator div [ (K df/dv grad g) dv ]
        if (S_->HasDerivative(fname, keyc)) {
          S_->GetEvaluator(fname).UpdateDerivative(*S_, passwd_, keyc, Tags::DEFAULT);
          kr_c = *S_->GetDerivative<CompositeVector>(fname, Tags::DEFAULT, keyc, Tags::DEFAULT)
                    .ViewComponent("cell");

          // --- upwind derivative
          S_->GetEvaluator(gname).Update(*S_, passwd_);
          auto var = S_->GetPtr<CompositeVector>(gname, Tags::DEFAULT);
          auto flux = S_->GetPtr<CompositeVector>(flux_names_[phase], Tags::DEFAULT);
          upwind_->Compute(*flux, bcnone, *kr);

          // --- calculate advective flux from div(K kr grad g)
          Teuchos::RCP<Operators::PDE_Diffusion> pde;
          pde = (type == 0) ? pde_diff_K_[phase] : pde_diff_D_;
          pde->SetScalarCoefficient(kr, Teuchos::null);
          pde->SetBCs(op_bcs_[gname], op_bcs_[gname]);
          pde->UpdateMatrices(flux.ptr(), var.ptr());
          pde->UpdateFlux(var.ptr(), flux_tmp.ptr());

          double factor = eqns_[row].factors[i];
          flux_acc->Update(factor, *flux_tmp, 1.0);
        }

        // populate advection operator. Note that upwind direction is based on
        // the Darcy flux, but flux value is based on derivatives of coefficients
        if (special_pc_term_) {
          auto pde1 = fac_adv_->Create(global_op);
          const auto& flux = S_->Get<CompositeVector>(flux_names_[phase]);
          pde1->Setup(flux);

          pde1->SetBCs(op_bcs_[keyr], op_bcs_[keyr]);
          pde1->UpdateMatrices(flux_acc.ptr(), [](double x) { return x; });
          pde1->ApplyBCs(false, false, false);
        }
      }

      // storage term
      if ((key = eqns_[row].storage) != "") {
        if (S_->HasDerivative(key, keyc)) {
          auto pde =
            Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::Entity_kind::CELL, global_op));
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
      auto pde =
        Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::Entity_kind::CELL, mesh_));
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
            if (der_gc.get() ) (*fone_c)[0][c] = (*der_gc)[0][c];
          } else {
            if (der_fc.get() ) (*fone_c)[0][c] = (*der_fc)[0][c];
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
  // auto J = FiniteDifferenceJacobian(tp - dtp, tp, u, u, 1e-6);
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
  num_ls_itrs_ += op_pc_solver_->num_itrs();
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

      double floor(1.0);
      for (int c = 0; c < ncells_owned_; c++) {
        for (int k = 0; k < n; ++k) {
          error_tmp = std::max(error_tmp, fabs(dxc[k][c]) / (xc[k][c] + floor));
        }
      }
    }

    else if (soln_names_[i] == saturation_liquid_key_) {
      auto& dsc = *du->SubVector(i)->Data()->ViewComponent("cell");

      for (int c = 0; c < ncells_owned_; c++) {
        error_tmp = std::max(error_tmp, fabs(dsc[0][c]));
      }
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
Multiphase_PK::FiniteDifferenceJacobian(double t_old,
                                        double t_new,
                                        Teuchos::RCP<const TreeVector> u_old,
                                        Teuchos::RCP<const TreeVector> u_new,
                                        double eps)
{
  auto f0 = Teuchos::rcp(new TreeVector(*u_old));
  auto f1 = Teuchos::rcp(new TreeVector(*u_old));
  auto u0 = Teuchos::rcp(new TreeVector(*u_old));
  auto u1 = Teuchos::rcp_const_cast<TreeVector>(u_new);

  double umax;
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
    u1_c.NormInf(&umax);
    double factor = eps * umax;
    if (n == 2) factor = -eps;
    if (factor == 0.0) continue;

    u1_c[0][c] += factor;
    FunctionalResidual(t_old, t_new, u0, u1, f1);
    u1_c[0][c] -= factor;

    for (int nrow = 0; nrow < nJ; ++nrow) {
      int m = nrow / ncells_owned_;
      int i = nrow % ncells_owned_;

      auto& f0_c = *f0->SubVector(m)->Data()->ViewComponent("cell");
      auto& f1_c = *f1->SubVector(m)->Data()->ViewComponent("cell");
      J(nrow, ncol) = (f1_c[0][i] - f0_c[0][i]) / factor;
    }
  }

  return J;
}

} // namespace Multiphase
} // namespace Amanzi
