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

  // trigger update of primary variables
  pressure_liquid_eval_->SetFieldAsChanged(S_.ptr());
  xl_liquid_eval_->SetFieldAsChanged(S_.ptr());
  saturation_liquid_eval_->SetFieldAsChanged(S_.ptr());

  // extract pointers to subvectors
  std::vector<Teuchos::RCP<CompositeVector> > up, fp;
  for (int i = 0; i < 3; ++i) {
    up.push_back(u_new->SubVector(i)->Data());
    fp.push_back(f->SubVector(i)->Data());
  }

  // miscalleneous fields
  // -- fluxes
  auto flux_liquid = S_->GetFieldData(darcy_flux_liquid_key_, passwd_);
  auto flux_gas = S_->GetFieldData(darcy_flux_gas_key_, passwd_);

  // -- molar mobilities
  S_->GetFieldEvaluator(molar_mobility_liquid_key_)->HasFieldChanged(S_.ptr(), molar_mobility_liquid_key_);
  auto mobility_l = S_->GetFieldData(molar_mobility_liquid_key_);
  const auto& mobility_lc = *mobility_l->ViewComponent("cell");

  S_->GetFieldEvaluator(molar_mobility_gas_key_)->HasFieldChanged(S_.ptr(), passwd_);
  auto mobility_g = S_->GetFieldData(molar_mobility_gas_key_, molar_mobility_gas_key_);
  const auto& mobility_gc = *mobility_g->ViewComponent("cell");

  // -- gas pressure
  S_->GetFieldEvaluator(pressure_gas_key_)->HasFieldChanged(S_.ptr(), passwd_);
  auto pg = S_->GetFieldData(pressure_gas_key_, pressure_gas_key_);

  // -- molar densities
  S_->GetFieldEvaluator(molar_density_gas_key_)->HasFieldChanged(S_.ptr(), passwd_);
  const auto& eta_g = S_->GetFieldData(molar_density_gas_key_);
  const auto& eta_gc = *eta_g->ViewComponent("cell");

  const auto& eta_l = S_->GetFieldData(molar_density_liquid_key_);
  const auto& eta_lc = *eta_l->ViewComponent("cell");

  // -- mass density of gas phase 
  auto rho_g = Teuchos::rcp(new CompositeVector(*eta_g));

  // work memory for miscalleneous operator
  auto kr = CreateCVforUpwind(mesh_);
  auto& kr_c = *kr->ViewComponent("cell");

  auto psol = PressureToSolution();
  auto tmp = up[psol.first];
  CompositeVector fone(*tmp), fadd(*tmp), comp(*tmp);
  auto& fone_c = *fone.ViewComponent("cell");

  // start loop over physical equations
  Key key;
  for (int n = 0; n < num_primary_ + 1; ++n) {
    auto csol = ComponentToSolution(n);
    const auto& uci = *up[csol.first]->ViewComponent("cell");

    auto eqn = EquationToSolution(n);
    PopulateBCs(eqn.second);
  
    // Richards operator for liquid phase
    if ((key = eval_mobility_liquid_[n]) != "") {
      // -- upwind product of liquid mobility, molar density, and mole fraction
      *kr->ViewComponent("cell") = *mobility_l->ViewComponent("cell");
      if (csol.second < 0) {
        kr_c = mobility_lc;  // FIXME
      } else {
        for (int c = 0; c < ncells_owned_; ++c)
          kr_c[0][c] = mobility_lc[0][c] * uci[csol.second][c];
      }
      kr->ViewComponent("dirichlet_faces")->PutScalar(eta_l_ / mu_l_);  // FIXME
      upwind_->Compute(*flux_liquid, *kr, op_bcs_[eqn.first]->bc_model(), *kr);

      // -- add transport in liquid phase
      Teuchos::RCP<std::vector<WhetStone::Tensor> > Kptr = Teuchos::rcpFromRef(K_);

      auto& pde = pde_diff_K_;
      pde->Setup(Kptr, kr, Teuchos::null, rho_l_, gravity_);
      pde->SetBCs(op_bcs_[eqn.first], op_bcs_[eqn.first]);
      pde->global_operator()->Init();
      pde->UpdateMatrices(Teuchos::null, Teuchos::null);
      pde->ApplyBCs(true, false, true);
      pde->global_operator()->ComputeNegativeResidual(*up[eqn.first], fone);
    } else {
      // primary boundary conditions should be imposed in each equation for at least
      // one term.
      AMANZI_ASSERT(true);
    }
 
    // Richards operator for gas phase
    if ((key = eval_mobility_gas_[n]) != "") {
      // -- upwind product of gas mobility, molar density, and mole fraction
      for (int c = 0; c < ncells_owned_; ++c) {
        kr_c[0][c] = mobility_gc[0][c] * (1.0 - uci[csol.second][c] / kH_[csol.second]);
      }
      kr->ViewComponent("dirichlet_faces")->PutScalar(0.0);  // FIXME
      upwind_->Compute(*flux_gas, *kr, op_bcs_[psol.first]->bc_model(), *kr);

      // -- add transport in gas phase
      Teuchos::RCP<std::vector<WhetStone::Tensor> > Kptr = Teuchos::rcpFromRef(K_);

      auto& pde = pde_diff_K_;
      pde->Setup(Kptr, kr, Teuchos::null, rho_g, gravity_);
      pde->SetBCs(op_bcs_[eqn.first], op_bc_pg_);
      pde->global_operator()->Init();
      pde->UpdateMatrices(Teuchos::null, Teuchos::null);
      pde->ApplyBCs(false, false, false);
      pde->global_operator()->ComputeNegativeResidual(*pg, fadd);
      fone.Update(1.0, fadd, 1.0);
    }

    // molecular diffusion via harmonic-mean transmissibility
    if ((key = eval_molecular_diff_[n]) != "") {
      int i = csol.second;
      *comp.ViewComponent("cell") = *uci(i);

      auto D = Teuchos::rcp(new std::vector<WhetStone::Tensor>);
      WhetStone::Tensor Dc(dim_, 1);
      for (int c = 0; c < ncells_owned_; ++c) {
        Dc(0, 0) = eta_lc[0][c] * mol_diff_l_[i] - eta_gc[0][c] * mol_diff_g_[i];
        D->push_back(Dc);
      }

      auto& pde = pde_diff_D_;
      pde->Setup(D, Teuchos::null, Teuchos::null);
      pde->SetBCs(op_bcs_[eqn.first], op_bcs_[eqn.first]);
      pde->global_operator()->Init();
      pde->UpdateMatrices(Teuchos::null, Teuchos::null);
      pde->ApplyBCs(false, false, false);
      pde->global_operator()->ComputeNegativeResidual(comp, fadd);
      fone.Update(1.0, fadd, 1.0);
    }


    // add counter-diffusion to the water component
    fp[psol.first]->Update(-1.0, fadd, 1.0);

    // add storage terms 
    if ((key = eval_storage_[n]) != "") {
      std::string prev_key = "prev_" + key;

      auto eval = S_->GetFieldEvaluator(key);
      if (key == tcs_key_) {
        Teuchos::rcp_dynamic_cast<TotalComponentStorage>(eval)->set_component_id(eqn.second);
        Teuchos::rcp_dynamic_cast<TotalComponentStorage>(eval)->set_kH(kH_[eqn.second]);
      }
      eval->HasFieldChanged(S_.ptr(), passwd_);

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

  // process gas constraint
  auto ssol = SaturationToSolution();
  auto& usi = *up[ssol.first]->ViewComponent("cell");  // saturation liquid

  auto csol = ComponentToSolution(num_primary_ + 1);
  auto& uci = *up[csol.first]->ViewComponent("cell");  // chemical components

  int n = num_primary_ - 1;
  if (ncp_ == "min") {
    for (int c = 0; c < ncells_owned_; ++c) {
      double x_sum(0.0);
      for (int i = 0; i < num_primary_; ++i) x_sum += uci[i][c] / kH_[i];

      uci[n][c] = std::min(1.0 - usi[0][c], 1.0 - x_sum);
    }
  } else if (ncp_ == "Fischer-Burmeister") {
    for (int c = 0; c < ncells_owned_; ++c) {
      double x_sum(0.0);
      for (int i = 0; i < num_primary_; ++i) x_sum += uci[i][c] / kH_[i];

      double a = 1.0 - usi[0][c];
      double b = 1.0 - x_sum;
      uci[n][c] = pow(a * a + b * b, 0.5) - (a + b);
    }
  }
}


/* ******************************************************************
* Update new preconditioner on the interval (tp-dtp, tp].
****************************************************************** */
void Multiphase_PK::UpdatePreconditioner(double tp, Teuchos::RCP<const TreeVector> u, double dtp)
{
  double t_old = tp - dtp;

  // trigger update of primary variables
  saturation_liquid_eval_->SetFieldAsChanged(S_.ptr());

  // extract pointers to subvectors
  std::vector<Teuchos::RCP<const CompositeVector> > up;
  for (int i = 0; i < 3; ++i) {
    up.push_back(u->SubVector(i)->Data());
  }

  // -- mobilities
  S_->GetFieldEvaluator(molar_mobility_liquid_key_)->HasFieldChanged(S_.ptr(), passwd_);
  auto mobility_l = S_->GetFieldData(molar_mobility_liquid_key_, molar_mobility_liquid_key_);
  const auto& mobility_lc = *mobility_l->ViewComponent("cell");

  // parameter lists
  auto& adv_list = mp_list_->sublist("operators").sublist("advection operator").sublist("preconditioner");
  auto& ddf_list = mp_list_->sublist("operators").sublist("diffusion operator").sublist("preconditioner");
  auto& mdf_list = mp_list_->sublist("operators").sublist("molecular diffusion operator").sublist("preconditioner");

  // miscalleneous fields
  auto flux_liquid = S_->GetFieldData(darcy_flux_liquid_key_, passwd_);
  auto flux_gas = S_->GetFieldData(darcy_flux_gas_key_, passwd_);

  // work memory for miscalleneous operator
  auto kr = CreateCVforUpwind(mesh_);
  auto& kr_c = *kr->ViewComponent("cell");
  auto& kr_f = *kr->ViewComponent("face");

  auto psol = PressureToSolution();
  auto tmp = up[psol.first];
  CompositeVector fone(*tmp);

  // for each operator we linearize (a) functions on which it acts and
  // (b) non-linear coefficients
  Key key;
  for (int row = 0; row < num_primary_ + 1; ++row) {
    auto csol = ComponentToSolution(row);
    const auto& uci = *up[csol.first]->ViewComponent("cell");

    auto eqnr = EquationToSolution(row);
    Key keyr = soln_names_[eqnr.first];
    PopulateBCs(eqnr.second);
  
    for (int col = 0; col < num_primary_ + 2; ++col) {
      auto eqnc = EquationToSolution(col);
      Key keyc = soln_names_[eqnc.first];

      // add empty operator to have a global operator
      auto pde = Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::CELL, mesh_));
      auto global_op = pde->global_operator();
      op_preconditioner_->SetOperatorBlock(row, col, global_op);
      kr_c.PutScalar(0.0);
      pde->AddAccumulationTerm(*kr, "cell");

      // Richards operator for liquid phase
      // -- (a) since liquid pressure is primary, its linearization is either 0 or 1
      if ((key = eval_mobility_liquid_[row]) != "" && keyc == pressure_liquid_key_) {
        auto pde = Teuchos::rcp(new Operators::PDE_DiffusionFV(ddf_list, global_op));

        // -- upwind product of liquid mobility and liquid mole fraction
        *kr->ViewComponent("cell") = *mobility_l->ViewComponent("cell");
        if (csol.second < 0) {
          kr_c = mobility_lc;  // FIXME
        } else {
          for (int c = 0; c < ncells_owned_; ++c)
            kr_c[0][c] = mobility_lc[0][c] * uci[csol.second][c];
        }
        kr->ViewComponent("dirichlet_faces")->PutScalar(eta_l_ / mu_l_);  // FIXME
        upwind_->Compute(*flux_liquid, *kr, op_bcs_[eqnr.first]->bc_model(), *kr);

        Teuchos::RCP<std::vector<WhetStone::Tensor> > Kptr = Teuchos::rcpFromRef(K_);
        pde->Setup(Kptr, kr, Teuchos::null);
        pde->SetBCs(op_bcs_[eqnr.first], op_bcs_[eqnr.first]);
        pde->global_operator()->Init();
        pde->UpdateMatrices(Teuchos::null, Teuchos::null);
        pde->ApplyBCs(true, true, true);
      } 

      // -- (b) molar mobility could be a function of saturation 
      if ((key = eval_mobility_liquid_[row]) != "") {
        Key der_key = "d" + key + "_d" + keyc;
        S_->GetFieldEvaluator(key)->HasFieldDerivativeChanged(S_.ptr(), passwd_, keyc);

        if (S_->HasField(der_key)) {
          auto pde = Teuchos::rcp(new Operators::PDE_AdvectionUpwind(adv_list, global_op)); 

          // --- upwind derivative
          auto der = S_->GetFieldData(der_key, key);
          *kr->ViewComponent("cell") = *der->ViewComponent("cell");
          kr->ViewComponent("dirichlet_faces")->PutScalar(0.0);  // FIXME
          upwind_->Compute(*flux_liquid, *kr, op_bcs_[0]->bc_model(), *kr);

          auto flux_tmp(*flux_liquid);
          auto& flux_f = *flux_tmp.ViewComponent("face");
          for (int f = 0; f < nfaces_owned_; ++f) flux_f[0][f] /= kr_f[0][f];

          pde->Setup(flux_tmp);
          pde->SetBCs(op_bcs_[eqnr.first], op_bcs_[eqnc.first]);
          pde->global_operator()->Init();
          pde->UpdateMatrices(flux_liquid.ptr());
          pde->ApplyBCs(false, false, false);
        }
      }

      // Richards operator for gas phase
      // -- (a) since gas pressure depends on liquid pressure and liquid saturation, 
      // its linearization adds a factor to molar mobility
      if ((key = eval_mobility_gas_[row]) != "") {
        Key der_key = "dpressure_gas_d" + keyc;
        S_->GetFieldEvaluator("pressure_gas")->HasFieldDerivativeChanged(S_.ptr(), passwd_, keyc);

        if (S_->HasField(der_key)) {
          auto pde = Teuchos::rcp(new Operators::PDE_DiffusionFV(ddf_list, global_op));

          // -- upwind product of liquid mobility, liquid mole fraction, and
          //    derivative gas pressure 
          const auto& der_c = *S_->GetFieldData(der_key)->ViewComponent("cell");
          for (int c = 0; c < ncells_owned_; ++c)
            kr_c[0][c] = mobility_lc[0][c] * der_c[0][c];

          kr->ViewComponent("dirichlet_faces")->PutScalar(0.0);  // FIXME
          upwind_->Compute(*flux_gas, *kr, op_bcs_[eqnr.first]->bc_model(), *kr);

          Teuchos::RCP<std::vector<WhetStone::Tensor> > Kptr = Teuchos::rcpFromRef(K_);
          pde->Setup(Kptr, kr, Teuchos::null);
          pde->SetBCs(op_bcs_[eqnr.first], op_bcs_[eqnr.first]);
          pde->global_operator()->Init();
          pde->UpdateMatrices(Teuchos::null, Teuchos::null);
          pde->ApplyBCs(false, false, false);
        }
      }

      // storage term
      if ((key = eval_storage_[row]) != "") {
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
  auto ssol = SaturationToSolution();
  auto& usi = *up[ssol.first]->ViewComponent("cell");  // saturation liquid

  auto csol = ComponentToSolution(num_primary_ + 1);
  auto& uci = *up[csol.first]->ViewComponent("cell");  // chemical components

  CompositeVector dfdx(fone), dfds(fone);
  auto& dfdx_c = *dfdx.ViewComponent("cell");
  auto& dfds_c = *dfds.ViewComponent("cell");

  for (int i = 0; i < num_primary_; ++i) {
    auto pde = Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::CELL, mesh_)); 
    op_preconditioner_->SetOperatorBlock(num_primary_ + 1, i + 2, pde->global_operator());
 
    // -- identify active set for gas phase
    std::vector<int> active_g(ncells_owned_, 0);
    std::vector<int> inactive_g(ncells_owned_, 1);

    for (int c = 0; c < ncells_owned_; c++) {
      double x_sum(0.0);
      for (int k = 0; k < num_primary_; ++k) x_sum += uci[k][c] / kH_[k];

      if (1.0 - usi[0][c] > 1.0 - x_sum + 1e-12) {
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
          dfdx_c[0][c] =-1.0 / kH_[i];
          dfds_c[0][c] = 0.0;
        }
      }
    }
    pde->AddAccumulationTerm(dfdx, "cell");
  }

  auto pde = Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::CELL, mesh_)); 
  op_preconditioner_->SetOperatorBlock(num_primary_ + 1, 1, pde->global_operator());
  pde->AddAccumulationTerm(dfds, "cell");

  // finalize preconditioner
  if (!op_pc_assembled_) {
    op_preconditioner_->SymbolicAssembleMatrix();
    op_pc_assembled_ = true;
  }
  op_preconditioner_->AssembleMatrix();
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
std::cout << "X: " << (*X->SubVector(1)->Data()->ViewComponent("cell"))[0][0] <<
                 "-> Y: " << (*Y->SubVector(1)->Data()->ViewComponent("cell"))[0][0] << std::endl;
  return ierr;
}


/* ******************************************************************
* Monitor l2 norm of residual
****************************************************************** */
double Multiphase_PK::ErrorNorm(Teuchos::RCP<const TreeVector> u,
                                Teuchos::RCP<const TreeVector> du) 
{
// { double aaa; u->SubVector(0)->Norm2(&aaa); std::cout << aaa << std::endl; }
// { std::cout << *du->SubVector(1)->Data()->ViewComponent("cell") << std::endl; }
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
std::cout << "DU: " << dsc[0][0] << " " << dsc[0][1] << std::endl;

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
double aaa = du1c[0][c];
    du1c[0][c] = std::min(du1c[0][c], u1c[0][c]);
    du1c[0][c] = std::max(du1c[0][c], u1c[0][c] - 1.0);
if (fabs(aaa-du1c[0][c]) > 1e-12 && c < 2) std::cout << " clipping in cell " << c << std::endl; 
  }

  return AmanziSolvers::FnBaseDefs::CORRECTION_MODIFIED;
}

}  // namespace Multiphase
}  // namespace Amanzi

