/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Derived MPC for flow and energy.  This couples using a block-diagonal coupler.
------------------------------------------------------------------------- */

#include "wrm_richards_evaluator.hh"
#include "mpc_frozen_prec_coupled_flow_energy.hh"

namespace Amanzi {

#define DEBUG_FLAG 0

RegisteredPKFactory<MPCFrozenCoupledFlowEnergy> MPCFrozenCoupledFlowEnergy::reg_("frozen energy-flow preconditioner coupled");

bool MPCFrozenCoupledFlowEnergy::modify_predictor(double h, Teuchos::RCP<TreeVector> u) {
  if (predictor_type_ == PREDICTOR_HEURISTIC) {
    return modify_predictor_heuristic(h,u);
  } else if (predictor_type_ == PREDICTOR_EWC) {
    if (S_->time() == 0.) { // this needs to be fixed!
      return false;
    } else {
      return modify_predictor_ewc(h,u);
    }
  } else if (predictor_type_ == PREDICTOR_EWC_HEURISTIC) {
    if (S_->time() == 0.) { // this needs to be fixed!
      return false;
    } else {
      return modify_predictor_ewc_heuristic(h,u);
    }
  }
  return false;
};

bool MPCFrozenCoupledFlowEnergy::modify_predictor_ewc_heuristic(double h, Teuchos::RCP<TreeVector> u) {

  Teuchos::OSTab tab = getOSTab();
  if (S_work_ != Teuchos::null) {
    Teuchos::RCP<CompositeVector> temp_guess = u->SubVector("energy")->data();
    Teuchos::RCP<CompositeVector> pres_guess = u->SubVector("flow")->data();
    Epetra_MultiVector& temp_guess_c = *temp_guess->ViewComponent("cell",false);
    Epetra_MultiVector& pres_guess_c = *pres_guess->ViewComponent("cell",false);

    Epetra_MultiVector& T2 = *S_work_->GetFieldData("temperature", energy_pk->name())
      ->ViewComponent("cell",false);
    Epetra_MultiVector& p2 = *S_work_->GetFieldData("pressure", flow_pk->name())
      ->ViewComponent("cell",false);

    const Epetra_MultiVector& T1 = *S_inter_->GetFieldData("temperature")
      ->ViewComponent("cell",false);
    const Epetra_MultiVector& p1 = *S_inter_->GetFieldData("pressure")
      ->ViewComponent("cell",false);

    // get the needed vals for pressure in saturated zone
    const Epetra_MultiVector& cv = *S_inter_->GetFieldData("cell_volume")
      ->ViewComponent("cell",false);
    const Epetra_MultiVector& phi = *S_next_->GetFieldData("porosity")
      ->ViewComponent("cell",false);
    const Epetra_MultiVector& n_g = *S_next_->GetFieldData("molar_density_gas")
      ->ViewComponent("cell",false);
    const Epetra_MultiVector& omega_g = *S_next_->GetFieldData("mol_frac_gas")
      ->ViewComponent("cell",false);
    const Epetra_MultiVector& n_l = *S_next_->GetFieldData("molar_density_liquid")
      ->ViewComponent("cell",false);
    const Epetra_MultiVector& n_i = *S_next_->GetFieldData("molar_density_ice")
      ->ViewComponent("cell",false);
    const Epetra_MultiVector& one_on_A = *S_next_->GetFieldData("wrm_permafrost_one_on_A")
      ->ViewComponent("cell",false);

    const double& p_atm = *S_next_->GetScalarData("atmospheric_pressure");

    // get the WRMs for the saturated zone
    Teuchos::RCP<FieldEvaluator> wrm_B_eval_base = S_next_->GetFieldEvaluator("wrm_permafrost_one_on_B");
    Teuchos::RCP<Flow::FlowRelations::WRMRichardsEvaluator> wrm_B_eval =
      Teuchos::rcp_dynamic_cast<Flow::FlowRelations::WRMRichardsEvaluator>(wrm_B_eval_base);
    Teuchos::RCP<Flow::FlowRelations::WRMRegionPairList> wrms = wrm_B_eval->get_WRMs();

    // UNSATURATED ZONE
    // work is done in S_inter -- get the T and p evaluators
    Teuchos::RCP<FieldEvaluator> Teval_fe = S_work_->GetFieldEvaluator("temperature");
    Teuchos::RCP<FieldEvaluator> peval_fe = S_work_->GetFieldEvaluator("pressure");
    Teuchos::RCP<PrimaryVariableFieldEvaluator> Teval =
      Teuchos::rcp_static_cast<PrimaryVariableFieldEvaluator>(Teval_fe);
    Teuchos::RCP<PrimaryVariableFieldEvaluator> peval =
      Teuchos::rcp_static_cast<PrimaryVariableFieldEvaluator>(peval_fe);

    // project energy and water content
    double dt_next = S_next_->time() - S_inter_->time();
    double dt_prev = S_inter_->time() - S_work_->time();

    if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_EXTREME, true)) {
      *out_ << "Projecting Energy and WC, dt_prev = " << dt_prev
            << ", dt_next = " << dt_next << std::endl;
    }

    // -- get wc and energy data
    const Epetra_MultiVector& wc0 = *S_work_->GetFieldData("prev_water_content")
      ->ViewComponent("cell",false);
    const Epetra_MultiVector& wc1 = *S_inter_->GetFieldData("water_content")
      ->ViewComponent("cell",false);
    Epetra_MultiVector& wc2 = *S_work_->GetFieldData("water_content", "water_content")
      ->ViewComponent("cell",false);

    const Epetra_MultiVector& e0 = *S_work_->GetFieldData("prev_energy")
      ->ViewComponent("cell",false);
    const Epetra_MultiVector& e1 = *S_inter_->GetFieldData("energy")
      ->ViewComponent("cell",false);
    Epetra_MultiVector& e2 = *S_work_->GetFieldData("energy", "energy")
      ->ViewComponent("cell",false);
    const Epetra_MultiVector& sat_g = *S_next_->GetFieldData("saturation_gas")
      ->ViewComponent("cell",false);

    // -- work vectors for target wc and energy
    Epetra_MultiVector wc_t(wc0), e_t(e0);
    Epetra_MultiVector cor_wc(wc0), cor_e(e0);

    // -- project
    wc_t = wc0;
    e_t = e0;
    double dt_ratio = (dt_next + dt_prev) / dt_prev;
    wc_t.Update(dt_ratio, wc1, 1. - dt_ratio);
    e_t.Update(dt_ratio, e1, 1. - dt_ratio);
    //    std::cout << "S0 wc,e: " << wc0[0][99] << ", " << e0[0][99] << std::endl;
    //    std::cout << "S1 wc,e: " << wc1[0][99] << ", " << e1[0][99] << std::endl;
    //    std::cout << "Desired wc,e: " << wc_t[0][99] << ", " << e_t[0][99] << std::endl;

    // -- copy guess for T,p, and the resulting e,wc into S_work_, which we will use for this
    T2 = T1;
    p2 = p1;

    // -- pull out derivatives from S_work_, and copy the values from S_inter_
    Epetra_MultiVector& de_dT = *S_work_->GetFieldData("denergy_dtemperature", "energy")
      ->ViewComponent("cell",false);
    de_dT = *S_inter_->GetFieldData("denergy_dtemperature")->ViewComponent("cell",false);
    Epetra_MultiVector& de_dp = *S_work_->GetFieldData("denergy_dpressure", "energy")
      ->ViewComponent("cell",false);
    de_dp = *S_inter_->GetFieldData("denergy_dpressure")->ViewComponent("cell",false);
    Epetra_MultiVector& dwc_dT = *S_work_->GetFieldData("dwater_content_dtemperature", "water_content")
      ->ViewComponent("cell",false);
    dwc_dT = *S_inter_->GetFieldData("dwater_content_dtemperature")->ViewComponent("cell",false);
    Epetra_MultiVector& dwc_dp = *S_work_->GetFieldData("dwater_content_dpressure", "water_content")
      ->ViewComponent("cell",false);
    dwc_dp = *S_inter_->GetFieldData("dwater_content_dpressure")->ViewComponent("cell",false);

    // initialize the nonlinear solve, evaluate the (relative) residual (stored in 0)
    int ncells = temp_guess->size("cell",false);
    for (int c=0; c!=ncells; ++c) {
      if (sat_g[0][c] > 1.e-2) {
        e2[0][c] = ( e1[0][c] - e_t[0][c]) / e_t[0][c];
        wc2[0][c] = ( wc1[0][c] - wc_t[0][c]) / wc_t[0][c];
      } else {
        e2[0][c] = 0.;
        wc2[0][c] = 0.;
      }
    }

    // get the initial (relative) residual norm
    double norm_e_prev, norm_wc_prev, norm_prev;
    wc2.NormInf(&norm_wc_prev);
    e2.NormInf(&norm_e_prev);
    norm_prev = std::max(norm_e_prev, norm_wc_prev);
    //    std::cout << "Init Norm = " << norm_prev << std::endl;

    // DEBUG find cell that is most of norm
    /*
    int badc = 0;
    for (int c=0; c!=ncells; ++c) {
      if (std::abs(std::abs(wc2[0][c]) - norm_prev) < 1.e-10 ||
          std::abs(std::abs(e2[0][c]) - norm_prev) < 1.e-10) {
        badc = c;
        std::cout << "Found bad cell: " << badc << std::endl;
      }
    }
    */

    // scale back to res, not relative res
    for (int c=0; c!=ncells; ++c) {
      e2[0][c] = e2[0][c] * e_t[0][c];
      wc2[0][c] = wc2[0][c] * wc_t[0][c];
    }

    //    std::cout << "Init p,T: " << p2[0][badc] << ", " << T2[0][badc] << std::endl;
    //    std::cout << "Init wc,e: " << wc1[0][badc] << ", " << e1[0][badc] << std::endl;
    //    std::cout << "Desired wc,e: " << wc_t[0][badc] << ", " << e_t[0][badc] << std::endl;

    // iterate nonlinear solve
    bool converged = false;
    while (!converged) {
      //    std::cout << "Res p,T: " << wc2[0][badc] << ", " << e2[0][badc] << std::endl;

      // cor <-- J^-1 * res
      for (int c=0; c!=ncells; ++c) {
        double dedp = de_dp[0][c];
        double dedT = de_dT[0][c];
        double dwcdp = dwc_dp[0][c];
        double dwcdT = dwc_dT[0][c];

        /*
        if (c == badc) {
          std::cout << "Derivs: " << dwcdp << ",        " << dwcdT << std::endl;
          std::cout << "        " << dedp << ",        " << dedT << std::endl;
        }
        */

        double detJ = dedT*dwcdp - dedp*dwcdT;
        if (std::abs(detJ) < 1.e-20) {
          std::cout << "Crashing, zero det in predictor cell " << c << std::endl;
          std::cout << "  p = " << p2[0][c] << std::endl;
          std::cout << "  T = " << T2[0][c] << std::endl;
          std::cout << "  de_dp = " << dedp << std::endl;
          std::cout << "  de_dT = " << dedT << std::endl;
          std::cout << "  dwc_dp = " << dwcdp << std::endl;
          std::cout << "  dwc_dT = " << dwcdT << std::endl;
          std::cout << "  detJ = " << detJ << std::endl;
          ASSERT(0);
        }

        double res_e = e2[0][c];
        double res_wc = wc2[0][c];

        cor_wc[0][c] = (dedT*res_wc - dwcdT*res_e) / detJ;
        cor_e[0][c] = (-dedp*res_wc + dwcdp*res_e) / detJ;
      }
      //    std::cout << "Cor p,T: " << cor_wc[0][badc] << ", " << cor_e[0][badc] << std::endl;


      // x_n+1 = x_n - cor
      T2.Update(-1.0, cor_e, 1.0);
      p2.Update(-1.0, cor_wc, 1.0);
      //    std::cout << "Next p,T: " << p2[0][badc] << ", " << T2[0][badc] << std::endl;

      // evaluate e(p,T), wc(p,T)
      Teval->SetFieldAsChanged();
      peval->SetFieldAsChanged();
      S_work_->GetFieldEvaluator("water_content")->HasFieldChanged(S_work_.ptr(),name_);
      S_work_->GetFieldEvaluator("energy")->HasFieldChanged(S_work_.ptr(),name_);

      // e(p,T) - e_t, the desired solution
      // e2/wc2 now have the updated values of e(p,T)
      //    std::cout << "Current wc,e: " << wc2[0][badc] << ", " << e2[0][badc] << std::endl;
      //    std::cout << "Desired wc,e: " << wc_t[0][badc] << ", " << e_t[0][badc] << std::endl;
      for (int c=0; c!=ncells; ++c) {
        if (sat_g[0][c] > 1.e-2) {
          e2[0][c] = ( e2[0][c] - e_t[0][c]) / e_t[0][c];
          wc2[0][c] = ( wc2[0][c] - wc_t[0][c]) / wc_t[0][c];
        } else {
          e2[0][c] = 0.;
          wc2[0][c] = 0.;
        }
      }

      // evaulate convergence
      double norm_e, norm_wc, norm_t;
      wc2.NormInf(&norm_wc);
      e2.NormInf(&norm_e);
      norm_t = std::max(norm_e, norm_wc);

      //    std::cout << "New Norm = " << norm_t << std::endl;

      // DEBUG find cell that is most of norm
      /*
      badc = 0;
      for (int c=0; c!=ncells; ++c) {
        if (std::abs(std::abs(wc2[0][c]) - norm_t) < 1.e-10 ||
            std::abs(std::abs(e2[0][c]) - norm_t) < 1.e-10) {
          badc = c;
          std::cout << "Found bad cell: " << badc << std::endl;
        }
      }
      */

      // scale back to res, not relative res
      for (int c=0; c!=ncells; ++c) {
        e2[0][c] = e2[0][c] * e_t[0][c];
        wc2[0][c] = wc2[0][c] * wc_t[0][c];
      }

      double damp = 1.;
      //    std::cout << "damping required? " << (norm_t > norm_prev) << std::endl;
      while (norm_t > norm_prev) { // backtrack
        damp = damp * .5;
        //    std::cout << std::endl;
        //    std::cout << "Backtracking with damping of " << damp << std::endl;
        // backtrack the correction
        T2.Update(damp, cor_e, 1.0);
        p2.Update(damp, cor_wc, 1.0);
        //    std::cout << "Damped p,T: " << p2[0][badc] << ", " << T2[0][badc] << std::endl;

        // evaluate e(p,T), wc(p,T)
        Teval->SetFieldAsChanged();
        peval->SetFieldAsChanged();
        S_work_->GetFieldEvaluator("water_content")->HasFieldChanged(S_work_.ptr(),name_);
        S_work_->GetFieldEvaluator("energy")->HasFieldChanged(S_work_.ptr(),name_);

        // update res
        //    std::cout << "Current wc,e: " << wc2[0][badc] << ", " << e2[0][badc] << std::endl;
        //    std::cout << "Desired wc,e: " << wc_t[0][badc] << ", " << e_t[0][badc] << std::endl;
        for (int c=0; c!=ncells; ++c) {
          if (sat_g[0][c] > 1.e-2) {
            e2[0][c] = ( e2[0][c] - e_t[0][c]) / e_t[0][c];
            wc2[0][c] = ( wc2[0][c] - wc_t[0][c]) / wc_t[0][c];
          } else {
            e2[0][c] = 0.;
            wc2[0][c] = 0.;
          }
        }

        // reevaluate norm
        wc2.NormInf(&norm_wc);
        e2.NormInf(&norm_e);
        norm_t = std::max(norm_e, norm_wc);
        //    std::cout << "Damped New Norm (target) = " << norm_t << ", " << norm_prev << std::endl;

        // scale back to res, not relative res
        for (int c=0; c!=ncells; ++c) {
          e2[0][c] = e2[0][c] * e_t[0][c];
          wc2[0][c] = wc2[0][c] * wc_t[0][c];
        }

        //    std::cout << "damping continued? " << (norm_t > norm_prev) << std::endl;
      }

      converged = (norm_t < 1.e-6);
      norm_prev = norm_t;

      if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_EXTREME, true)) {
        *out_ << "   Solve for T,p: res = " << norm_t << std::endl;
      }

      if (converged) {
        if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_EXTREME, true)) {
          *out_ << " Converged." << std::endl;
        }

      } else {
        // update derivatives
        S_work_->GetFieldEvaluator("water_content")
          ->HasFieldDerivativeChanged(S_work_.ptr(),name_,"temperature");
        S_work_->GetFieldEvaluator("water_content")
          ->HasFieldDerivativeChanged(S_work_.ptr(),name_,"pressure");
        S_work_->GetFieldEvaluator("energy")
          ->HasFieldDerivativeChanged(S_work_.ptr(),name_,"temperature");
        S_work_->GetFieldEvaluator("energy")
          ->HasFieldDerivativeChanged(S_work_.ptr(),name_,"pressure");
      }
    }

    std::vector<bool> changed(ncells, false);

    for (int c=0; c!=ncells; ++c) {
      if (sat_g[0][c] > 1.e-2) {
        // IF SATURATED: copy solution into the result
        temp_guess_c[0][c] = T2[0][c];
        pres_guess_c[0][c] = p2[0][c];
      } else {
        // IF NOT SATURATED: use the heuristic
        if (T1[0][c] >= 273.15 && temp_guess_c[0][c] < 273.15) {
          if (includesVerbLevel(verbosity_, Teuchos::VERB_EXTREME, true)) {
            std::cout << "Freezing cell " << c << std::endl;
            std::cout << "   T_prev = " << T1[0][c]
                      << ",  T_guess = " << temp_guess_c[0][c]
                      << ",  T_corrected = " << 273.15 - 1.e-3 << std::endl;

          }
          // freezing
          temp_guess_c[0][c] = 273.15 - 1.e-3;
          changed[c] = true;

        } else if (273.15 > T1[0][c] &&
                   T1[0][c] >= 273.1 &&
                   (T1[0][c] - temp_guess_c[0][c]) > 1.e-2) {
          if (includesVerbLevel(verbosity_, Teuchos::VERB_EXTREME, true)) {
            std::cout << "2nd Freezing step cell " << c << std::endl;
            std::cout << "   T_prev = " << T1[0][c]
                      << ",  T_guess = " << temp_guess_c[0][c]
                      << ",  T_corrected = " << T1[0][c] << std::endl;
          }

          // catch the 2nd step in freezing -- after the 2nd step the
          // extrapolation should be ok?
          temp_guess_c[0][c] = T1[0][c];
          changed[c] = true;

        } else if (T1[0][c] <= 273.15 && temp_guess_c[0][c] > 273.15) {
          if (modify_thaw_to_prev_) {
            if (includesVerbLevel(verbosity_, Teuchos::VERB_EXTREME, true)) {
              std::cout << "Thawing cell " << c << std::endl;
              std::cout << "   T_prev = " << T1[0][c]
                        << ",  T_guess = " << temp_guess_c[0][c]
                        << ",  T_corrected = " << T1[0][c] << std::endl;
            }

            // thawing
            temp_guess_c[0][c] = T1[0][c];
            changed[c] = true;

          } else {
            if (includesVerbLevel(verbosity_, Teuchos::VERB_EXTREME, true)) {
              std::cout << "Thawing cell " << c << std::endl;
              std::cout << "   T_prev = " << T1[0][c]
                        << ",  T_guess = " << temp_guess_c[0][c]
                        << ",  T_corrected = " << 273.15 - 1.e-3 << std::endl;
            }

            // thawing
            temp_guess_c[0][c] = std::min((temp_guess_c[0][c] + T1[0][c]) / 2., 273.14);
            changed[c] = true;

          }
        }
      }
    }

    // Update pressure on saturated cells
    // update water content, which will get all the needed vals updated at the new temp.
    energy_pk->changed_solution();
    S_next_->GetFieldEvaluator("water_content")->HasFieldChanged(S_next_.ptr(), name_);

    for (Flow::FlowRelations::WRMRegionPairList::iterator region=wrms->begin();
         region!=wrms->end(); ++region) {
      std::string name = region->first;
      int nregion_cells = pres_guess->mesh()->get_set_size(name, AmanziMesh::CELL, AmanziMesh::OWNED);
      AmanziMesh::Entity_ID_List cells(nregion_cells);
      pres_guess->mesh()->get_set_entities(name, AmanziMesh::CELL, AmanziMesh::OWNED, &cells);

      for (AmanziMesh::Entity_ID_List::iterator c=cells.begin(); c!=cells.end(); ++c) {
        if (changed[*c]) {
          double A_minus_one = (1.0/one_on_A[0][*c] - 1.0);

          double wc = wc0[0][*c] / cv[0][*c];
          double sstar = (wc - n_g[0][*c]*omega_g[0][*c]*phi[0][*c]) /
            (phi[0][*c] * (n_l[0][*c] - n_g[0][*c]*omega_g[0][*c]
                    + n_i[0][*c]*A_minus_one) - A_minus_one*wc);

          if (sstar > 0.) {
            double pc = region->second->capillaryPressure(sstar);

            //          if (includesVerbLevel(verbosity_, Teuchos::VERB_EXTREME, true)) {
            std::cout << "   p_prev = " << p1[0][*c]
                      << "   p_guess = " << pres_guess_c[0][*c]
                      << "   p_corrected = " << p_atm - pc << std::endl;
            std::cout << "  (based upon T = " << T1[0][*c] << ")" << std::endl;
            //          }

            pres_guess_c[0][*c] = p_atm - pc;
          }
        }
      }
    }

    // update faces from cells
    temp_guess->ScatterMasterToGhosted("cell");
    pres_guess->ScatterMasterToGhosted("cell");
    AmanziMesh::Entity_ID_List cells;

    int f_owned = temp_guess->size("face");
    for (int f=0; f!=f_owned; ++f) {
      cells.clear();
      temp_guess->mesh()->face_get_cells(f, AmanziMesh::USED, &cells);
      if (cells.size() == 2) {
        (*temp_guess)("face",f) = ((*temp_guess)("cell",cells[0]) +
                (*temp_guess)("cell",cells[1])) / 2.0;
        (*pres_guess)("face",f) = ((*pres_guess)("cell",cells[0]) +
                (*pres_guess)("cell",cells[1])) / 2.0;
      } else {
        (*temp_guess)("face",f) = (*temp_guess)("cell",cells[0]);
        (*pres_guess)("face",f) = (*pres_guess)("cell",cells[0]);
      }
    }

    return true;
  } else {
    return false;
  }
}


bool MPCFrozenCoupledFlowEnergy::modify_predictor_ewc(double h, Teuchos::RCP<TreeVector> u) {

  Teuchos::OSTab tab = getOSTab();
  if (S_work_ != Teuchos::null) {
    Teuchos::RCP<CompositeVector> temp_guess = u->SubVector("energy")->data();
    Teuchos::RCP<CompositeVector> pres_guess = u->SubVector("flow")->data();
    Epetra_MultiVector& temp_guess_c = *temp_guess->ViewComponent("cell",false);
    Epetra_MultiVector& pres_guess_c = *pres_guess->ViewComponent("cell",false);

    Epetra_MultiVector& T2 = *S_work_->GetFieldData("temperature", energy_pk->name())
      ->ViewComponent("cell",false);
    Epetra_MultiVector& p2 = *S_work_->GetFieldData("pressure", flow_pk->name())
      ->ViewComponent("cell",false);

    const Epetra_MultiVector& T1 = *S_inter_->GetFieldData("temperature")
      ->ViewComponent("cell",false);
    const Epetra_MultiVector& p1 = *S_inter_->GetFieldData("pressure")
      ->ViewComponent("cell",false);

    const Epetra_MultiVector& sat_g = *S_next_->GetFieldData("saturation_gas")
      ->ViewComponent("cell",false);

    // UNSATURATED ZONE
    // work is done in S_inter -- get the T and p evaluators
    Teuchos::RCP<FieldEvaluator> Teval_fe = S_work_->GetFieldEvaluator("temperature");
    Teuchos::RCP<FieldEvaluator> peval_fe = S_work_->GetFieldEvaluator("pressure");
    Teuchos::RCP<PrimaryVariableFieldEvaluator> Teval =
      Teuchos::rcp_static_cast<PrimaryVariableFieldEvaluator>(Teval_fe);
    Teuchos::RCP<PrimaryVariableFieldEvaluator> peval =
      Teuchos::rcp_static_cast<PrimaryVariableFieldEvaluator>(peval_fe);

    // project energy and water content
    double dt_next = S_next_->time() - S_inter_->time();
    double dt_prev = S_inter_->time() - S_work_->time();

    if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_EXTREME, true)) {
      *out_ << "Projecting Energy and WC, dt_prev = " << dt_prev
            << ", dt_next = " << dt_next << std::endl;
    }

    // -- get wc and energy data
    const Epetra_MultiVector& wc0 = *S_work_->GetFieldData("prev_water_content")
      ->ViewComponent("cell",false);
    const Epetra_MultiVector& wc1 = *S_inter_->GetFieldData("water_content")
      ->ViewComponent("cell",false);
    Epetra_MultiVector& wc2 = *S_work_->GetFieldData("water_content", "water_content")
      ->ViewComponent("cell",false);

    const Epetra_MultiVector& e0 = *S_work_->GetFieldData("prev_energy")
      ->ViewComponent("cell",false);
    const Epetra_MultiVector& e1 = *S_inter_->GetFieldData("energy")
      ->ViewComponent("cell",false);
    Epetra_MultiVector& e2 = *S_work_->GetFieldData("energy", "energy")
      ->ViewComponent("cell",false);

    // -- work vectors for target wc and energy
    Epetra_MultiVector wc_t(wc0), e_t(e0);
    Epetra_MultiVector cor_wc(wc0), cor_e(e0);

    // -- project
    wc_t = wc0;
    e_t = e0;
    double dt_ratio = (dt_next + dt_prev) / dt_prev;
    wc_t.Update(dt_ratio, wc1, 1. - dt_ratio);
    e_t.Update(dt_ratio, e1, 1. - dt_ratio);
    std::cout << "S0 wc,e: " << wc0[0][99] << ", " << e0[0][99] << std::endl;
    std::cout << "S1 wc,e: " << wc1[0][99] << ", " << e1[0][99] << std::endl;
    std::cout << "Desired wc,e: " << wc_t[0][99] << ", " << e_t[0][99] << std::endl;

    // -- copy guess for T,p, and the resulting e,wc into S_work_, which we will use for this
    T2 = T1;
    p2 = p1;

    // -- pull out derivatives from S_work_, and copy the values from S_inter_
    Epetra_MultiVector& de_dT = *S_work_->GetFieldData("denergy_dtemperature", "energy")
      ->ViewComponent("cell",false);
    de_dT = *S_inter_->GetFieldData("denergy_dtemperature")->ViewComponent("cell",false);
    Epetra_MultiVector& de_dp = *S_work_->GetFieldData("denergy_dpressure", "energy")
      ->ViewComponent("cell",false);
    de_dp = *S_inter_->GetFieldData("denergy_dpressure")->ViewComponent("cell",false);
    Epetra_MultiVector& dwc_dT = *S_work_->GetFieldData("dwater_content_dtemperature", "water_content")
      ->ViewComponent("cell",false);
    dwc_dT = *S_inter_->GetFieldData("dwater_content_dtemperature")->ViewComponent("cell",false);
    Epetra_MultiVector& dwc_dp = *S_work_->GetFieldData("dwater_content_dpressure", "water_content")
      ->ViewComponent("cell",false);
    dwc_dp = *S_inter_->GetFieldData("dwater_content_dpressure")->ViewComponent("cell",false);

    // initialize the nonlinear solve, evaluate the (relative) residual (stored in 0)
    int ncells = temp_guess->size("cell",false);
    for (int c=0; c!=ncells; ++c) {
      if (sat_g[0][c] > 1.e-2) {
        e2[0][c] = ( e1[0][c] - e_t[0][c]) / e_t[0][c];
        wc2[0][c] = ( wc1[0][c] - wc_t[0][c]) / wc_t[0][c];
      } else {
        e2[0][c] = 0.;
        wc2[0][c] = 0.;
      }
    }

    // get the initial (relative) residual norm
    double norm_e_prev, norm_wc_prev, norm_prev;
    wc2.NormInf(&norm_wc_prev);
    e2.NormInf(&norm_e_prev);
    norm_prev = std::max(norm_e_prev, norm_wc_prev);
    //    std::cout << "Init Norm = " << norm_prev << std::endl;

    // DEBUG find cell that is most of norm
    /*
    int badc = 0;
    for (int c=0; c!=ncells; ++c) {
      if (std::abs(std::abs(wc2[0][c]) - norm_prev) < 1.e-10 ||
          std::abs(std::abs(e2[0][c]) - norm_prev) < 1.e-10) {
        badc = c;
        //    std::cout << "Found bad cell: " << badc << std::endl;
      }
    }
    */

    // scale back to res, not relative res
    for (int c=0; c!=ncells; ++c) {
      e2[0][c] = e2[0][c] * e_t[0][c];
      wc2[0][c] = wc2[0][c] * wc_t[0][c];
    }

    //    std::cout << "Init p,T: " << p2[0][badc] << ", " << T2[0][badc] << std::endl;
    //    std::cout << "Init wc,e: " << wc1[0][badc] << ", " << e1[0][badc] << std::endl;
    //    std::cout << "Desired wc,e: " << wc_t[0][badc] << ", " << e_t[0][badc] << std::endl;

    // iterate nonlinear solve
    bool converged = false;
    while (!converged) {
      //    std::cout << "Res p,T: " << wc2[0][badc] << ", " << e2[0][badc] << std::endl;

      // cor <-- J^-1 * res
      for (int c=0; c!=ncells; ++c) {
        double dedp = de_dp[0][c];
        double dedT = de_dT[0][c];
        double dwcdp = dwc_dp[0][c];
        double dwcdT = dwc_dT[0][c];

        /*
        if (c == badc) {
          std::cout << "Derivs: " << dwcdp << ",        " << dwcdT << std::endl;
          std::cout << "        " << dedp << ",        " << dedT << std::endl;
        }
        */

        double detJ = dedT*dwcdp - dedp*dwcdT;
        if (std::abs(detJ) < 1.e-20) {
          std::cout << "Crashing, zero det in predictor cell " << c << std::endl;
          std::cout << "  p = " << p2[0][c] << std::endl;
          std::cout << "  T = " << T2[0][c] << std::endl;
          std::cout << "  de_dp = " << dedp << std::endl;
          std::cout << "  de_dT = " << dedT << std::endl;
          std::cout << "  dwc_dp = " << dwcdp << std::endl;
          std::cout << "  dwc_dT = " << dwcdT << std::endl;
          std::cout << "  detJ = " << detJ << std::endl;
          ASSERT(0);
        }

        double res_e = e2[0][c];
        double res_wc = wc2[0][c];

        cor_wc[0][c] = (dedT*res_wc - dwcdT*res_e) / detJ;
        cor_e[0][c] = (-dedp*res_wc + dwcdp*res_e) / detJ;
      }
      //    std::cout << "Cor p,T: " << cor_wc[0][badc] << ", " << cor_e[0][badc] << std::endl;


      // x_n+1 = x_n - cor
      T2.Update(-1.0, cor_e, 1.0);
      p2.Update(-1.0, cor_wc, 1.0);
      //    std::cout << "Next p,T: " << p2[0][badc] << ", " << T2[0][badc] << std::endl;

      // evaluate e(p,T), wc(p,T)
      Teval->SetFieldAsChanged();
      peval->SetFieldAsChanged();
      S_work_->GetFieldEvaluator("water_content")->HasFieldChanged(S_work_.ptr(),name_);
      S_work_->GetFieldEvaluator("energy")->HasFieldChanged(S_work_.ptr(),name_);

      // e(p,T) - e_t, the desired solution
      // e2/wc2 now have the updated values of e(p,T)
      //    std::cout << "Current wc,e: " << wc2[0][badc] << ", " << e2[0][badc] << std::endl;
      //    std::cout << "Desired wc,e: " << wc_t[0][badc] << ", " << e_t[0][badc] << std::endl;
      for (int c=0; c!=ncells; ++c) {
        if (sat_g[0][c] > 1.e-2) {
          e2[0][c] = ( e2[0][c] - e_t[0][c]) / e_t[0][c];
          wc2[0][c] = ( wc2[0][c] - wc_t[0][c]) / wc_t[0][c];
        } else {
          e2[0][c] = 0.;
          wc2[0][c] = 0.;
        }
      }

      // evaulate convergence
      double norm_e, norm_wc, norm_t;
      wc2.NormInf(&norm_wc);
      e2.NormInf(&norm_e);
      norm_t = std::max(norm_e, norm_wc);

      //    std::cout << "New Norm = " << norm_t << std::endl;

      // DEBUG find cell that is most of norm
      /*
      badc = 0;
      for (int c=0; c!=ncells; ++c) {
        if (std::abs(std::abs(wc2[0][c]) - norm_t) < 1.e-10 ||
            std::abs(std::abs(e2[0][c]) - norm_t) < 1.e-10) {
          badc = c;
          std::cout << "Found bad cell: " << badc << std::endl;
        }
      }
      */

      // scale back to res, not relative res
      for (int c=0; c!=ncells; ++c) {
        e2[0][c] = e2[0][c] * e_t[0][c];
        wc2[0][c] = wc2[0][c] * wc_t[0][c];
      }

      double damp = 1.;
      //    std::cout << "damping required? " << (norm_t > norm_prev) << std::endl;
      while (norm_t > norm_prev) { // backtrack
        damp = damp * .5;
        //    std::cout << std::endl;
        //    std::cout << "Backtracking with damping of " << damp << std::endl;
        // backtrack the correction
        T2.Update(damp, cor_e, 1.0);
        p2.Update(damp, cor_wc, 1.0);
        //    std::cout << "Damped p,T: " << p2[0][badc] << ", " << T2[0][badc] << std::endl;

        // evaluate e(p,T), wc(p,T)
        Teval->SetFieldAsChanged();
        peval->SetFieldAsChanged();
        S_work_->GetFieldEvaluator("water_content")->HasFieldChanged(S_work_.ptr(),name_);
        S_work_->GetFieldEvaluator("energy")->HasFieldChanged(S_work_.ptr(),name_);

        // update res
        //    std::cout << "Current wc,e: " << wc2[0][badc] << ", " << e2[0][badc] << std::endl;
        //    std::cout << "Desired wc,e: " << wc_t[0][badc] << ", " << e_t[0][badc] << std::endl;
        for (int c=0; c!=ncells; ++c) {
          if (sat_g[0][c] > 1.e-2) {
            e2[0][c] = ( e2[0][c] - e_t[0][c]) / e_t[0][c];
            wc2[0][c] = ( wc2[0][c] - wc_t[0][c]) / wc_t[0][c];
          } else {
            e2[0][c] = 0.;
            wc2[0][c] = 0.;
          }
        }

        // reevaluate norm
        wc2.NormInf(&norm_wc);
        e2.NormInf(&norm_e);
        norm_t = std::max(norm_e, norm_wc);
        //    std::cout << "Damped New Norm (target) = " << norm_t << ", " << norm_prev << std::endl;

        // scale back to res, not relative res
        for (int c=0; c!=ncells; ++c) {
          e2[0][c] = e2[0][c] * e_t[0][c];
          wc2[0][c] = wc2[0][c] * wc_t[0][c];
        }

        //    std::cout << "damping continued? " << (norm_t > norm_prev) << std::endl;
      }

      converged = (norm_t < 1.e-6);
      norm_prev = norm_t;

      if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_EXTREME, true)) {
        *out_ << "   Solve for T,p: res = " << norm_t << std::endl;
      }

      if (converged) {
        if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_EXTREME, true)) {
          *out_ << " Converged." << std::endl;
        }

      } else {
        // update derivatives
        S_work_->GetFieldEvaluator("water_content")
          ->HasFieldDerivativeChanged(S_work_.ptr(),name_,"temperature");
        S_work_->GetFieldEvaluator("water_content")
          ->HasFieldDerivativeChanged(S_work_.ptr(),name_,"pressure");
        S_work_->GetFieldEvaluator("energy")
          ->HasFieldDerivativeChanged(S_work_.ptr(),name_,"temperature");
        S_work_->GetFieldEvaluator("energy")
          ->HasFieldDerivativeChanged(S_work_.ptr(),name_,"pressure");
      }
    }

    std::vector<bool> changed(ncells, false);

    for (int c=0; c!=ncells; ++c) {
      if (sat_g[0][c] > 1.e-2) {
        // IF SATURATED: copy solution into the result
        temp_guess_c[0][c] = T2[0][c];
        pres_guess_c[0][c] = p2[0][c];
      }
    }

    // update faces from cells
    temp_guess->ScatterMasterToGhosted("cell");
    pres_guess->ScatterMasterToGhosted("cell");
    AmanziMesh::Entity_ID_List cells;

    int f_owned = temp_guess->size("face");
    for (int f=0; f!=f_owned; ++f) {
      cells.clear();
      temp_guess->mesh()->face_get_cells(f, AmanziMesh::USED, &cells);
      if (cells.size() == 2) {
        (*temp_guess)("face",f) = ((*temp_guess)("cell",cells[0]) +
                (*temp_guess)("cell",cells[1])) / 2.0;
        (*pres_guess)("face",f) = ((*pres_guess)("cell",cells[0]) +
                (*pres_guess)("cell",cells[1])) / 2.0;
      } else {
        (*temp_guess)("face",f) = (*temp_guess)("cell",cells[0]);
        (*pres_guess)("face",f) = (*pres_guess)("cell",cells[0]);
      }
    }

    return true;
  } else {
    return false;
  }
}


bool MPCFrozenCoupledFlowEnergy::modify_predictor_heuristic(double h, Teuchos::RCP<TreeVector> u) {
  Teuchos::RCP<CompositeVector> temp_guess = u->SubVector("energy")->data();
  Epetra_MultiVector& temp_guess_c = *temp_guess->ViewComponent("cell",false);
  Epetra_MultiVector& guess_faces = *temp_guess->ViewComponent("face",false);

  Teuchos::RCP<const CompositeVector> temp = S_inter_->GetFieldData("temperature");
  Teuchos::RCP<const CompositeVector> pres = S_inter_->GetFieldData("pressure");
  const Epetra_MultiVector& temp_c = *temp->ViewComponent("cell",false);
  const Epetra_MultiVector& temp_faces = *temp->ViewComponent("face",false);

  bool update_faces(false);


  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_EXTREME, true)) {
    double minT,maxT,minTL,maxTL;
    temp_guess_c.MinValue(&minT);
    temp_guess_c.MaxValue(&maxT);
    guess_faces.MinValue(&minTL);
    guess_faces.MaxValue(&maxTL);

    *out_ << "--- Modifying Guess: ---" << std::endl;
    *out_ << "  T_min = " << minT << "     T_max = " << maxT << std::endl;
    *out_ << "  Lambda_min = " << minTL << "     Lambda_max = " << maxTL << std::endl;
  }


  int ncells = temp->size("cell",false);
  std::vector<bool> changed(ncells, false);

  for (int c=0; c!=ncells; ++c) {
    if (temp_c[0][c] >= 273.15 && temp_guess_c[0][c] < 273.15) {
      //      if (includesVerbLevel(verbosity_, Teuchos::VERB_EXTREME, true)) {
        std::cout << "Freezing cell " << c << std::endl;
        std::cout << "   T_prev = " << temp_c[0][c]
              << ",  T_guess = " << temp_guess_c[0][c]
              << ",  T_corrected = " << 273.15 - 1.e-3 << std::endl;

        //      }
      // freezing
      temp_guess_c[0][c] = 273.15 - 1.e-3;
      changed[c] = true;
      update_faces = true;

    } else if (273.15 > temp_c[0][c] &&
               temp_c[0][c] >= 273.1 &&
               (temp_c[0][c] - temp_guess_c[0][c]) > 1.e-2) {
      //      if (includesVerbLevel(verbosity_, Teuchos::VERB_EXTREME, true)) {
        std::cout << "2nd Freezing step cell " << c << std::endl;
        std::cout << "   T_prev = " << temp_c[0][c]
              << ",  T_guess = " << temp_guess_c[0][c]
              << ",  T_corrected = " << temp_c[0][c] << std::endl;
        //      }

      // catch the 2nd step in freezing -- after the 2nd step the
      // extrapolation should be ok?
      temp_guess_c[0][c] = temp_c[0][c];
      changed[c] = true;
      update_faces = true;

    } else if (temp_c[0][c] <= 273.15 && temp_guess_c[0][c] > 273.15) {
      if (modify_thaw_to_prev_) {
        //      if (includesVerbLevel(verbosity_, Teuchos::VERB_EXTREME, true)) {
        std::cout << "Thawing cell " << c << std::endl;
        std::cout << "   T_prev = " << temp_c[0][c]
                  << ",  T_guess = " << temp_guess_c[0][c]
                  << ",  T_corrected = " << temp_c[0][c] << std::endl;
        //      }

        // thawing
        temp_guess_c[0][c] = temp_c[0][c];
        changed[c] = true;
        update_faces = true;
      } else {
        //      if (includesVerbLevel(verbosity_, Teuchos::VERB_EXTREME, true)) {
        std::cout << "Thawing cell " << c << std::endl;
        std::cout << "   T_prev = " << temp_c[0][c]
                  << ",  T_guess = " << temp_guess_c[0][c]
                  << ",  T_corrected = " << 273.15 - 1.e-3 << std::endl;
        //      }

        // thawing
        temp_guess_c[0][c] = std::min((temp_guess_c[0][c] + temp_c[0][c]) / 2., 273.14);
        changed[c] = true;
        update_faces = true;
      }
    }
  }


  // update pressure to be consistent with the new temperature and fixed mass
  // Stuff temperature into state
  //  Teuchos::RCP<CompositeVector> temp_from_state = S_next_->GetFieldData("temperature", energy_pk->name());
  //  S_next_->SetData("temperature",energy_pk->name(),temp_guess);
  energy_pk->changed_solution();

  // update water content, which will get all the needed vals updated at the new temp.
  S_next_->GetFieldEvaluator("water_content")
      ->HasFieldChanged(S_next_.ptr(), "richards_pk");

  // get the needed vals
  Teuchos::RCP<const CompositeVector> wc0 = S_inter_->GetFieldData("water_content");
  Teuchos::RCP<const CompositeVector> cv = S_inter_->GetFieldData("cell_volume");
  Teuchos::RCP<const CompositeVector> phi = S_next_->GetFieldData("porosity");
  Teuchos::RCP<const CompositeVector> n_g = S_next_->GetFieldData("molar_density_gas");
  Teuchos::RCP<const CompositeVector> omega_g = S_next_->GetFieldData("mol_frac_gas");
  Teuchos::RCP<const CompositeVector> n_l = S_next_->GetFieldData("molar_density_liquid");
  Teuchos::RCP<const CompositeVector> n_i = S_next_->GetFieldData("molar_density_ice");
  Teuchos::RCP<const CompositeVector> one_on_A = S_next_->GetFieldData("wrm_permafrost_one_on_A");
  Teuchos::RCP<const double> p_atm = S_next_->GetScalarData("atmospheric_pressure");

  // get the WRMs
  Teuchos::RCP<FieldEvaluator> wrm_B_eval_base = S_next_->GetFieldEvaluator("wrm_permafrost_one_on_B");
  Teuchos::RCP<Flow::FlowRelations::WRMRichardsEvaluator> wrm_B_eval =
      Teuchos::rcp_dynamic_cast<Flow::FlowRelations::WRMRichardsEvaluator>(wrm_B_eval_base);
  Teuchos::RCP<Flow::FlowRelations::WRMRegionPairList> wrms = wrm_B_eval->get_WRMs();

  // get the result pressure
  Teuchos::RCP<CompositeVector> pres_guess = u->SubVector("flow")->data();
  for (Flow::FlowRelations::WRMRegionPairList::iterator region=wrms->begin();
       region!=wrms->end(); ++region) {
    std::string name = region->first;
    int ncells = pres_guess->mesh()->get_set_size(name, AmanziMesh::CELL, AmanziMesh::OWNED);
    AmanziMesh::Entity_ID_List cells(ncells);
    pres_guess->mesh()->get_set_entities(name, AmanziMesh::CELL, AmanziMesh::OWNED, &cells);

    for (AmanziMesh::Entity_ID_List::iterator c=cells.begin(); c!=cells.end(); ++c) {
      if (changed[*c]) {
        double A_minus_one = (1.0/(*one_on_A)("cell",*c) - 1.0);

        double wc = (*wc0)("cell",*c) / (*cv)("cell",*c);
        double sstar = (wc - (*n_g)("cell",*c)*(*omega_g)("cell",*c)*(*phi)("cell",*c)) /
            ((*phi)("cell",*c) * ((*n_l)("cell",*c) - (*n_g)("cell",*c)*(*omega_g)("cell",*c)
                    + (*n_i)("cell",*c)*A_minus_one) - A_minus_one*wc);

        if (sstar > 0.) {
          double pc = region->second->capillaryPressure(sstar);

          //          if (includesVerbLevel(verbosity_, Teuchos::VERB_EXTREME, true)) {
            std::cout << "   p_prev = " << (*pres)("cell",*c)
                  << "   p_guess = " << (*pres_guess)("cell",*c)
                  << "   p_corrected = " << *p_atm - pc << std::endl;
            std::cout << "  (based upon T = " << (*temp)("cell",*c) << ")" << std::endl;
            //          }

          (*pres_guess)("cell",*c) = *p_atm - pc;

        }
      }
    }
  }


  // unclear... could do an AllReduce() on update_faces and not communicate if
  // not necessary, but AllReduces() suck... maybe more than extra Scatters?
  temp_guess->ScatterMasterToGhosted("cell");
  pres_guess->ScatterMasterToGhosted("cell");
  AmanziMesh::Entity_ID_List cells;

  int f_owned = temp_guess->size("face");
  for (int f=0; f!=f_owned; ++f) {
    cells.clear();
    temp_guess->mesh()->face_get_cells(f, AmanziMesh::USED, &cells);
    if (cells.size() == 2) {
      (*temp_guess)("face",f) = ((*temp_guess)("cell",cells[0]) +
              (*temp_guess)("cell",cells[1])) / 2.0;
      (*pres_guess)("face",f) = ((*pres_guess)("cell",cells[0]) +
              (*pres_guess)("cell",cells[1])) / 2.0;
    } else {
      (*temp_guess)("face",f) = (*temp_guess)("cell",cells[0]);
      (*pres_guess)("face",f) = (*pres_guess)("cell",cells[0]);
    }
  }

  // clean up -- undo the change of vectors in state
  //  S_next_->SetData("temperature",energy_pk->name(),temp_from_state);

  return true;
}

bool MPCFrozenCoupledFlowEnergy::modify_predictor_temp(double h, Teuchos::RCP<TreeVector> u) {
  // modification of the initial guess occurs by updating T using this guess,
  // then calculating the p that would be required to keep the same water
  // mass.

  // Stuff temperature into state
  Teuchos::RCP<TreeVector> uT = u->SubVector("energy");
  sub_pks_[1]->changed_solution();

#if DEBUG_FLAG
  std::cout << std::endl;
  std::cout << "Modifying guess:" << std::endl;
  std::cout << "  T0: " << (*uT->data())("cell",0) << std::endl;
  std::cout << "  T1: " << (*uT->data())("cell",99) << std::endl;
#endif

  // update water content, which will get all the needed vals updated at the new temp.
  S_next_->GetFieldEvaluator("water_content")
      ->HasFieldChanged(S_next_.ptr(), "richards_pk");

  // get the needed vals
  Teuchos::RCP<const CompositeVector> wc0 = S_inter_->GetFieldData("water_content");
  Teuchos::RCP<const CompositeVector> cv = S_inter_->GetFieldData("cell_volume");
  Teuchos::RCP<const CompositeVector> phi = S_next_->GetFieldData("porosity");
  Teuchos::RCP<const CompositeVector> n_g = S_next_->GetFieldData("molar_density_gas");
  Teuchos::RCP<const CompositeVector> omega_g = S_next_->GetFieldData("mol_frac_gas");
  Teuchos::RCP<const CompositeVector> n_l = S_next_->GetFieldData("molar_density_liquid");
  Teuchos::RCP<const CompositeVector> n_i = S_next_->GetFieldData("molar_density_ice");
  Teuchos::RCP<const CompositeVector> one_on_A = S_next_->GetFieldData("wrm_permafrost_one_on_A");
  Teuchos::RCP<const double> p_atm = S_next_->GetScalarData("atmospheric_pressure");

  // get the WRMs
  Teuchos::RCP<FieldEvaluator> wrm_B_eval_base = S_next_->GetFieldEvaluator("wrm_permafrost_one_on_B");
  Teuchos::RCP<Flow::FlowRelations::WRMRichardsEvaluator> wrm_B_eval =
      Teuchos::rcp_dynamic_cast<Flow::FlowRelations::WRMRichardsEvaluator>(wrm_B_eval_base);
  Teuchos::RCP<Flow::FlowRelations::WRMRegionPairList> wrms = wrm_B_eval->get_WRMs();


  // get the result pressure
  Teuchos::RCP<CompositeVector> pres = u->SubVector("flow")->data();
  for (Flow::FlowRelations::WRMRegionPairList::iterator region=wrms->begin();
       region!=wrms->end(); ++region) {
    std::string name = region->first;
    int ncells = pres->mesh()->get_set_size(name, AmanziMesh::CELL, AmanziMesh::OWNED);
    AmanziMesh::Entity_ID_List cells(ncells);
    pres->mesh()->get_set_entities(name, AmanziMesh::CELL, AmanziMesh::OWNED, &cells);

    for (AmanziMesh::Entity_ID_List::iterator c=cells.begin(); c!=cells.end(); ++c) {
      double p = (*pres)("cell",*c);
      double A_minus_one = (1.0/(*one_on_A)("cell",*c) - 1.0);

      double wc = (*wc0)("cell",*c) / (*cv)("cell",*c);
      double sstar = (wc - (*n_g)("cell",*c)*(*omega_g)("cell",*c)*(*phi)("cell",*c)) /
          ((*phi)("cell",*c) * ((*n_l)("cell",*c) - (*n_g)("cell",*c)*(*omega_g)("cell",*c)
                  + (*n_i)("cell",*c)*A_minus_one) - A_minus_one*wc);

      if (sstar > 0.) {
        if (*c==0) std::cout << "   A-1(0) = " << A_minus_one << std::endl;
        if (*c==99) std::cout << "   A-1(99) = " << A_minus_one << std::endl;
        if (*c==0) std::cout << "   S*(0) = " << sstar << std::endl;
        if (*c==99) std::cout << "   S*(99) = " << sstar << std::endl;

        double pc = region->second->capillaryPressure(sstar);

        if (*c==0) std::cout << "   pc(0) = " << pc << std::endl;
        if (*c==99) std::cout << "   pc(99) = " << pc << std::endl;

        (*pres)("cell",*c) = *p_atm - pc;

        if (*c==0) std::cout << "   p(0) = " << (*pres)("cell",*c) << std::endl;
        if (*c==99) std::cout << "   p(99) = " <<  (*pres)("cell",*c) << std::endl;

      }
    }
  }

  AmanziMesh::Entity_ID_List cells;
  pres->ScatterMasterToGhosted("cell");

  int f_owned = pres->size("face");
  for (int f=0; f!=f_owned; ++f) {
    cells.clear();
    pres->mesh()->face_get_cells(f, AmanziMesh::USED, &cells);
    int ncells = cells.size();

    double face_value = 0.0;
    for (int n=0; n!=ncells; ++n) {
      face_value += (*pres)("cell",cells[n]);
    }
    (*pres)("face",f) = face_value / ncells;
  }

  return true;
};


void MPCFrozenCoupledFlowEnergy::commit_state(double dt, const Teuchos::RCP<State>& S) {
  MPCCoupledFlowEnergy::commit_state(dt,S);

  if ((predictor_type_ == PREDICTOR_EWC || predictor_type_ == PREDICTOR_EWC_HEURISTIC)
       && S_inter_ != Teuchos::null) {
    // stash water content and energy in S_work.
    *S_work_->GetFieldData("prev_water_content",name_) =
      *S_inter_->GetFieldData("water_content");
    *S_work_->GetFieldData("prev_energy",name_) =
      *S_inter_->GetFieldData("energy");
    S_work_->set_time(S_inter_->time());
  }
}


void MPCFrozenCoupledFlowEnergy::setup(const Teuchos::Ptr<State>& S) {
  MPCCoupledFlowEnergy::setup(S);
  flow_pk = Teuchos::rcp_dynamic_cast<Amanzi::Flow::Richards>(sub_pks_[0]);
  energy_pk = Teuchos::rcp_dynamic_cast<Amanzi::Energy::TwoPhase>(sub_pks_[1]);

  std::string predictor_string = plist_.get<std::string>("predictor type", "heuristic");
  if (predictor_string == "none") {
    predictor_type_ = PREDICTOR_NONE;
  } else if (predictor_string == "heuristic") {
    predictor_type_ = PREDICTOR_HEURISTIC;
  } else if (predictor_string == "ewc") {
    predictor_type_ = PREDICTOR_EWC;
    S->RequireField("prev_water_content",name_)->SetMesh(S->GetMesh())->SetGhosted(false)
      ->SetComponent("cell",AmanziMesh::CELL,1);
    S->RequireField("prev_energy",name_)->SetMesh(S->GetMesh())->SetGhosted(false)
      ->SetComponent("cell",AmanziMesh::CELL,1);
  } else if (predictor_string == "ewc and heuristic") {
    predictor_type_ = PREDICTOR_EWC_HEURISTIC;
    S->RequireField("prev_water_content",name_)->SetMesh(S->GetMesh())->SetGhosted(false)
      ->SetComponent("cell",AmanziMesh::CELL,1);
    S->RequireField("prev_energy",name_)->SetMesh(S->GetMesh())->SetGhosted(false)
      ->SetComponent("cell",AmanziMesh::CELL,1);
  } else {
    Errors::Message message(std::string("Invalid predictor type ")+predictor_string);
    Exceptions::amanzi_throw(message);
  }
}

// -- Initialize owned (dependent) variables.
void MPCFrozenCoupledFlowEnergy::initialize(const Teuchos::Ptr<State>& S) {
  Teuchos::RCP<CompositeVector> pres = S->GetFieldData("pressure", flow_pk->name());
  Teuchos::RCP<CompositeVector> temp = S->GetFieldData("temperature", energy_pk->name());

  modify_thaw_to_prev_ = plist_.get<bool>("modify thawing cells to previous temp", true);


  if (plist_.get<bool>("initialize from frozen column", false)) {
    int npoints = 100;
    double ref_T[] = {262.7091428,262.641515665,262.58610644,262.54456626,262.516834608,262.503140657,262.496062716,262.488999243,262.48192875,262.474849863,262.467762477,262.460666527,262.453561996,262.446448892,262.439327226,262.432196991,262.425058094,262.417910659,262.410754807,262.403590438,262.396417523,262.38923607,262.382046106,262.374847613,262.367640513,262.360424583,262.353200074,262.345967205,262.338725906,262.331476326,262.324218348,262.316951874,262.309676754,262.302392788,262.295100307,262.287799762,262.280490954,262.273173738,262.26584795,262.258513003,262.25116961,262.243818573,262.236459336,262.229091793,262.221715879,262.214331552,262.206938887,262.199537835,262.192127993,262.184709455,262.177282338,262.169845193,262.162392238,262.154903957,262.147215752,262.138831145,262.127574119,262.106832867,262.070124703,262.01912569,261.95878179,261.891743368,261.819681263,261.743550148,261.664006281,261.581541241,261.496529897,261.409277942,261.320038524,261.229029146,261.136443739,261.042456713,260.947225061,260.850890294,260.753580116,260.655409681,260.556482544,260.456891469,260.356719229,260.256039443,260.154917462,260.053411276,259.951572408,259.849446765,259.747075431,259.644495381,259.541740117,259.438840219,259.335823827,259.232717041,259.129544246,259.026328383,258.92309113,258.819853031,258.716633546,258.613451025,258.510322629,258.407264182,258.304290022,258.201412884};
    double ref_p[] = {-91144793.8313,-79198250.2123,-64553378.2237,-53817268.5748,-40376600.31,-9186846.56641,-6245241.77642,-5844385.55992,-5808523.22231,-5797105.81683,-5795148.67225,-5807088.92269,-5810167.68087,-5821627.31323,-5820860.37694,-5836780.45965,-5868415.98829,-5840987.06814,-5838227.17932,-5843736.56489,-5852397.98775,-5856537.87728,-5857887.33366,-5869968.3054,-5899256.57528,-5985661.00795,-5927751.02996,-5956801.49746,-5915396.93489,-5903438.58908,-5900368.18448,-5923657.7528,-5971066.03108,-6058214.67597,-5996882.73853,-5952227.64174,-5956081.90559,-5958021.93473,-6019441.59966,-6202372.07811,-6030838.97887,-5995013.54658,-5986821.9308,-5986523.53408,-6003380.78559,-6018426.01063,-6014029.63473,-6050617.29614,-6174214.01679,-6175627.87735,-6271705.74182,-6657225.50566,-7802907.37546,-9925062.27734,-15631485.3766,-21637300.6575,-32255189.5224,-47929250.3074,-63007803.9699,-72942365.2584,-81113258.8592,-87994540.9417,-94323588.0452,-100411902.799,-106347449.951,-112258287.406,-118190275.792,-124185038.048,-130274367.97,-136457417.992,-142720313.132,-149043914.067,-155403433.988,-161767751.951,-168100565.584,-174362689.812,-180514544.613,-186518257.197,-192339144.137,-197946557.031,-203314195.973,-208420043.01,-213246066.64,-217777821.856,-222004036.177,-225916241.135,-229508484.983,-232777146.221,-235720857.447,-238340542.772,-240639567.369,-242623992.387,-244302919.818,-245688897.123,-246798326.804,-247651788.05,-248274123.31,-248694073.748,-248943176.096,-249053590.406};

    double dz = 0.1;
    int ref_water_table_k = 53;

    double water_table = plist_.get<double>("water table height");

    Epetra_MultiVector& temp_c = *temp->ViewComponent("cell",false);
    Epetra_MultiVector& pres_c = *pres->ViewComponent("cell",false);

    int ncells = temp->size("cell",false);
    for (int c=0; c!=ncells; ++c) {
      double c_z = pres->mesh()->cell_centroid(c)[2];
      double z_over_wt = c_z - water_table;
      int dk = std::floor(z_over_wt / dz);
      int k = dk + ref_water_table_k;

      if (k < 0 || k > npoints - 2) {
        Errors::Message message("Initialization of frozen column is off the charts");
        Exceptions::amanzi_throw(message);
      }

      double p0 = ref_p[k];
      double p1 = ref_p[k+1];
      double T0 = ref_T[k];
      double T1 = ref_T[k+1];

      double param = (z_over_wt - dz*dk) / dz;

      temp_c[0][c] = T0*param + T1*(1.-param);
      pres_c[0][c] = p0*param + p1*(1.-param);

    }

    AmanziMesh::Entity_ID_List cells;
    temp->ScatterMasterToGhosted("cell");
    pres->ScatterMasterToGhosted("cell");

    int f_owned = temp->size("face");
    for (int f=0; f!=f_owned; ++f) {
      cells.clear();
      pres->mesh()->face_get_cells(f, AmanziMesh::USED, &cells);
      int ncells = cells.size();

      double face_pres = 0.0;
      double face_temp = 0.0;
      for (int n=0; n!=ncells; ++n) {
        face_temp += (*temp)("cell",cells[n]);
        face_pres += (*pres)("cell",cells[n]);
      }
      (*temp)("face",f) = face_temp / ncells;
      (*pres)("face",f) = face_pres / ncells;
    }

    S->GetField("temperature",energy_pk->name())->set_initialized();
    S->GetField("pressure",flow_pk->name())->set_initialized();
  }

  MPCCoupledFlowEnergy::initialize(S);

  if (predictor_type_ == PREDICTOR_EWC || predictor_type_ == PREDICTOR_EWC_HEURISTIC) {
    *S->GetFieldData("prev_water_content",name_) = *S->GetFieldData("water_content");
    S->GetField("prev_water_content",name_)->set_initialized();
    *S->GetFieldData("prev_energy",name_) = *S->GetFieldData("energy");
    S->GetField("prev_energy",name_)->set_initialized();

    // work is done in S_inter -- get the T and p evaluators
    Teuchos::RCP<FieldEvaluator> Teval_fe = S->GetFieldEvaluator("temperature");
    Teuchos::RCP<FieldEvaluator> peval_fe = S->GetFieldEvaluator("pressure");
    Teuchos::RCP<PrimaryVariableFieldEvaluator> Teval =
      Teuchos::rcp_static_cast<PrimaryVariableFieldEvaluator>(Teval_fe);
    Teuchos::RCP<PrimaryVariableFieldEvaluator> peval =
      Teuchos::rcp_static_cast<PrimaryVariableFieldEvaluator>(peval_fe);
    Teval->SetFieldAsChanged();
    peval->SetFieldAsChanged();
  }
}


void MPCFrozenCoupledFlowEnergy::fun(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
         Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> g) {

  MPCCoupledFlowEnergy::fun(t_old, t_new, u_old, u_new, g);
  g->Norm2(&the_res_norm_);

  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_EXTREME, true)) {
    *out_ << "Setting backtracking res: " << the_res_norm_ << std::endl;
  }
}

bool MPCFrozenCoupledFlowEnergy::is_admissible(Teuchos::RCP<const TreeVector> up) {
  Teuchos::OSTab tab = getOSTab();
  if (!MPCCoupledFlowEnergy::is_admissible(up)) {
    return false;
  }

  if (backtracking_) {
    if (backtracking_count_ == backtracking_iterations_ - 1) {
      if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_EXTREME, true)) {
        *out_ << "Backtracking exausted, giving in." << std::endl;
      }
      backtracking_count_ = 0;
      return true;
    }

    double res_prev = the_res_norm_;
    Teuchos::RCP<TreeVector> up_nc = Teuchos::rcp_const_cast<TreeVector>(up);
    Teuchos::RCP<TreeVector> g = Teuchos::rcp(new TreeVector(*up));

    fun(S_inter_->time(), S_next_->time(), Teuchos::null, up_nc, g);

    if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_EXTREME, true)) {
      *out_ << "Checking admissibility for backtracking (prev, new): " << res_prev << ",  " << the_res_norm_ << std::endl;
    }

    if (the_res_norm_ <= res_prev) {
      backtracking_count_ = 0;
      return true;
    } else {
      backtracking_count_++;
      the_res_norm_ = res_prev;
      return false;
    }
  }

  return true;
}


void MPCFrozenCoupledFlowEnergy::set_states(const Teuchos::RCP<const State>& S,
          const Teuchos::RCP<State>& S_inter,
        const Teuchos::RCP<State>& S_next) {
  MPCCoupledFlowEnergy::set_states(S,S_inter,S_next);

  // update derivatives -- this ensures the derivatives exist in all states, not just S_next
  S_next_->GetFieldEvaluator("water_content")
    ->HasFieldDerivativeChanged(S_next_.ptr(),name_,"temperature");
  S_next_->GetFieldEvaluator("water_content")
    ->HasFieldDerivativeChanged(S_next_.ptr(),name_,"pressure");
  S_next_->GetFieldEvaluator("energy")
    ->HasFieldDerivativeChanged(S_next_.ptr(),name_,"temperature");
  S_next_->GetFieldEvaluator("energy")
    ->HasFieldDerivativeChanged(S_next_.ptr(),name_,"pressure");

  if (predictor_type_ == PREDICTOR_EWC || predictor_type_ == PREDICTOR_EWC_HEURISTIC) {
    S_work_ = Teuchos::rcp(new State(*S_next_));
    *S_work_ = *S_next_;
  }
}

} // namespace
