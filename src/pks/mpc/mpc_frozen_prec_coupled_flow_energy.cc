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
    return modify_predictor_ewc(h,u);
  }
  return false;
};


bool MPCFrozenCoupledFlowEnergy::modify_predictor_ewc(double h, Teuchos::RCP<TreeVector> u) {
  Teuchos::RCP<CompositeVector> temp_guess = u->SubVector("energy")->data();
  Teuchos::RCP<CompositeVector> pres_guess = u->SubVector("flow")->data();
  *temp_guess = *S_next_->GetFieldData("temperature_prediction");
  *pres_guess = *S_next_->GetFieldData("pressure_prediction");
  return true;
}


bool MPCFrozenCoupledFlowEnergy::modify_predictor_heuristic(double h, Teuchos::RCP<TreeVector> u) {
  Teuchos::RCP<CompositeVector> temp_guess = u->SubVector("energy")->data();
  Epetra_MultiVector& guess_cells = *temp_guess->ViewComponent("cell",false);
  Epetra_MultiVector& guess_faces = *temp_guess->ViewComponent("face",false);

  Teuchos::RCP<const CompositeVector> temp = S_next_->GetFieldData("temperature");
  const Epetra_MultiVector& temp_cells = *temp->ViewComponent("cell",false);

  bool update_faces(false);

#if DEBUG_FLAG
  std::cout << "--- Modifying Guess: ---" << std::cout;
#endif

  int ncells = temp->size("cell",false);
  std::vector<bool> changed(ncells, false);
  for (int c=0; c!=ncells; ++c) {
    if (temp_cells[0][c] >= 273.15 && guess_cells[0][c] < 273.15) {
      // freezing
      guess_cells[0][c] = 273.15 - 1.e-3;
      changed[c] = true;
      update_faces = true;
    } else if (temp_cells[0][c] <= 273.15 && guess_cells[0][c] > 273.15) {
      // thawing
      guess_cells[0][c] = 273.15 - 1.e-3;
      changed[c] = true;
      update_faces = true;
    } else if (273.15 > temp_cells[0][c] &&
               temp_cells[0][c] >= 273.1 &&
               (temp_cells[0][c] - guess_cells[0][c]) > 1.e-2) {
      // catch the 2nd step in freezing -- after the 2nd step the
      // extrapolation should be ok?
      guess_cells[0][c] = temp_cells[0][c];
      changed[c] = true;
      update_faces = true;
    }
  }

  // unclear... could do an AllReduce() on update_faces and not communicate if
  // not necessary, but AllReduces() suck... maybe more than extra Scatters?
  /*
  temp_guess->ScatterMasterToGhosted("cell");
  if (update_faces) {
    AmanziMesh::Entity_ID_List cells;

    int f_owned = temp_guess->size("face");
    for (int f=0; f!=f_owned; ++f) {
      cells.clear();
      temp_guess->mesh()->face_get_cells(f, AmanziMesh::USED, &cells);
      if (cells.size() == 2) {
        if (changed[cells[0]] || changed[cells[1]]) {
          guess_faces[0][f] = (guess_cells[0][cells[0]]
                  + guess_cells[0][cells[1]]) / 2.0;
        }
      }
    }
    return true;
  }
  */
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
  Teuchos::OSTab tab = getOSTab();

  if (predictor_type_ == PREDICTOR_EWC) {
    // work is done in S_inter -- get the T and p evaluators
    Teuchos::RCP<FieldEvaluator> Teval_fe = S_inter_->GetFieldEvaluator("temperature");
    Teuchos::RCP<FieldEvaluator> peval_fe = S_inter_->GetFieldEvaluator("pressure");
    Teuchos::RCP<PrimaryVariableFieldEvaluator> Teval =
      Teuchos::rcp_static_cast<PrimaryVariableFieldEvaluator>(Teval_fe);
    Teuchos::RCP<PrimaryVariableFieldEvaluator> peval =
      Teuchos::rcp_static_cast<PrimaryVariableFieldEvaluator>(peval_fe);

    // project energy and water content
    double dt_next = dt_;
    double dt_prev = S_next_->time() - S_inter_->time();

    if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_EXTREME, true)) {
      *out_ << "Projecting Energy and WC, dt_prev = " << dt_prev
            << ", dt_next = " << dt_next << std::endl;
    }

    // -- get wc and energy data
    Teuchos::RCP<CompositeVector> wc0 = S_inter_->GetFieldData("water_content", "water_content");
    Teuchos::RCP<const CompositeVector> wc1 = S_next_->GetFieldData("water_content", "water_content");
    Teuchos::RCP<CompositeVector> e0 = S_inter_->GetFieldData("energy", "energy");
    Teuchos::RCP<const CompositeVector> e1 = S_next_->GetFieldData("energy", "energy");

    // -- work vectors for wc and energy
    CompositeVector wc2(*wc0), e2(*e0);
    wc2 = *wc0;
    e2 = *e0;
    CompositeVector cor_wc(*wc0), cor_e(*e0);

    // -- project
    double dt_ratio = (dt_next + dt_prev) / dt_prev;
    wc2.Update(1.0 - dt_ratio, *wc1, dt_ratio);
    e2.Update(1.0 - dt_ratio, *e1, dt_ratio);

    // -- copy guess for T,p, and the resulting e,wc into S_inter_, which we will use for this
    Teuchos::RCP<CompositeVector> T0 = S_inter_->GetFieldData("temperature", energy_pk->name());
    Teuchos::RCP<CompositeVector> p0 = S_inter_->GetFieldData("pressure", flow_pk->name());
    Teuchos::RCP<const CompositeVector> T1 = S_next_->GetFieldData("temperature");
    Teuchos::RCP<const CompositeVector> p1 = S_next_->GetFieldData("pressure");
    *T0 = *T1;
    *p0 = *p1;

    // -- pull out derivatives from S_inter_, and copy the values from S_next_
    Teuchos::RCP<CompositeVector> de_dT = S_inter_->GetFieldData("denergy_dtemperature", "energy");
    *de_dT = *S_next_->GetFieldData("denergy_dtemperature");
    Teuchos::RCP<CompositeVector> de_dp = S_inter_->GetFieldData("denergy_dpressure", "energy");
    *de_dp = *S_next_->GetFieldData("denergy_dpressure");
    Teuchos::RCP<CompositeVector> dwc_dT = S_inter_->GetFieldData("dwater_content_dtemperature", "water_content");
    *dwc_dT = *S_next_->GetFieldData("dwater_content_dtemperature");
    Teuchos::RCP<CompositeVector> dwc_dp = S_inter_->GetFieldData("dwater_content_dpressure", "water_content");
    *dwc_dp = *S_next_->GetFieldData("dwater_content_dpressure");

    // initialize the nonlinear solve, evaluate the residual (stored in 0)
    *e0 = *e1;
    e0->Update(-1., e2, 1.);
    *wc0 = *wc1;
    wc0->Update(-1., wc2, 1.);

    int ncells = e0->size("cell",false);
    bool converged = false;

    // iterate nonlinear solve
    while (!converged) {
      // cor <-- J^-1 * res
      for (int c=0; c!=ncells; ++c) {
        double dedp = (*de_dp)("cell",c);
        double dedT = (*de_dT)("cell",c);
        double dwcdp = (*dwc_dp)("cell",c);
        double dwcdT = (*dwc_dT)("cell",c);

        double detJ = dedT*dwcdp - dedp*dwcdT;
        if (detJ < 1.e-14) {
          ASSERT(0);
        }

        double res_e = (*e0)("cell",c);
        double res_wc = (*wc0)("cell",c);

        cor_wc("cell",c) = (dedT*res_wc - dwcdT*res_e) / detJ;
        cor_e("cell",c) = (-dedp*res_wc + dwcdp*res_e) / detJ;
      }

      // x_n+1 = x_n - cor
      T0->Update(-1.0, cor_e, 1.0);
      p0->Update(-1.0, cor_wc, 1.0);

      // evaluate e(p,T), wc(p,T)
      Teval->SetFieldAsChanged();
      peval->SetFieldAsChanged();
      S_inter_->GetFieldEvaluator("water_content")->HasFieldChanged(S_inter_.ptr(),name_);
      S_inter_->GetFieldEvaluator("energy")->HasFieldChanged(S_inter_.ptr(),name_);

      // e(p,T) - e2, the desired solution
      // e0/wc0 now have the updated values
      for (int c=0; c!=ncells; ++c) {
        (*e0)("cell",c) = ( (*e0)("cell",c) - e2("cell",c)) / e2("cell",c);
        (*wc0)("cell",c) = ( (*wc0)("cell",c) - wc2("cell",c)) / wc2("cell",c);
      }

      // evaulate convergence
      double norm_e, norm_wc, norm;
      wc0->NormInf(&norm_wc);
      e0->NormInf(&norm_e);
      norm = std::min(norm_e, norm_wc);
      converged = (norm < 1.e-4);

      if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_EXTREME, true)) {
        *out_ << "   Solve for T,p: res = " << norm << std::endl;
      }

      if (converged) {
        if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_EXTREME, true)) {
          *out_ << " Converged." << std::endl;
        }
        // copy solution into the result
        *S_next_->GetFieldData("pressure_prediction",name_) = *p0;
        *S_next_->GetFieldData("temperature_prediction",name_) = *T0;

      } else {
        // update derivatives
        S_inter_->GetFieldEvaluator("dwater_content_dtemperature")
          ->HasFieldDerivativeChanged(S_inter_.ptr(),name_,"temperature");
        S_inter_->GetFieldEvaluator("dwater_content_dpressure")
          ->HasFieldDerivativeChanged(S_inter_.ptr(),name_,"pressure");
        S_inter_->GetFieldEvaluator("denergy_dtemperature")
          ->HasFieldDerivativeChanged(S_inter_.ptr(),name_,"temperature");
        S_inter_->GetFieldEvaluator("denergy_dpressure")
          ->HasFieldDerivativeChanged(S_inter_.ptr(),name_,"pressure");
      }
    }
  }
};

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
    S->RequireField("pressure_prediction",name_)->SetMesh(S->GetMesh())->SetGhosted(false)
      ->SetComponent("cell",AmanziMesh::CELL,1);
    S->RequireField("temperature_prediction",name_)->SetMesh(S->GetMesh())->SetGhosted(false)
      ->SetComponent("cell",AmanziMesh::CELL,1);
  } else {
    Errors::Message message(std::string("Invalid predictor type ")+predictor_string);
    Exceptions::amanzi_throw(message);
  }
}

// -- Initialize owned (dependent) variables.
void MPCFrozenCoupledFlowEnergy::initialize(const Teuchos::Ptr<State>& S) {

  if (plist_.get<bool>("initialize from frozen column", false)) {
    int npoints = 100;
    double ref_T[] = {262.7091428,262.641515665,262.58610644,262.54456626,262.516834608,262.503140657,262.496062716,262.488999243,262.48192875,262.474849863,262.467762477,262.460666527,262.453561996,262.446448892,262.439327226,262.432196991,262.425058094,262.417910659,262.410754807,262.403590438,262.396417523,262.38923607,262.382046106,262.374847613,262.367640513,262.360424583,262.353200074,262.345967205,262.338725906,262.331476326,262.324218348,262.316951874,262.309676754,262.302392788,262.295100307,262.287799762,262.280490954,262.273173738,262.26584795,262.258513003,262.25116961,262.243818573,262.236459336,262.229091793,262.221715879,262.214331552,262.206938887,262.199537835,262.192127993,262.184709455,262.177282338,262.169845193,262.162392238,262.154903957,262.147215752,262.138831145,262.127574119,262.106832867,262.070124703,262.01912569,261.95878179,261.891743368,261.819681263,261.743550148,261.664006281,261.581541241,261.496529897,261.409277942,261.320038524,261.229029146,261.136443739,261.042456713,260.947225061,260.850890294,260.753580116,260.655409681,260.556482544,260.456891469,260.356719229,260.256039443,260.154917462,260.053411276,259.951572408,259.849446765,259.747075431,259.644495381,259.541740117,259.438840219,259.335823827,259.232717041,259.129544246,259.026328383,258.92309113,258.819853031,258.716633546,258.613451025,258.510322629,258.407264182,258.304290022,258.201412884};
    double ref_p[] = {-91144793.8313,-79198250.2123,-64553378.2237,-53817268.5748,-40376600.31,-9186846.56641,-6245241.77642,-5844385.55992,-5808523.22231,-5797105.81683,-5795148.67225,-5807088.92269,-5810167.68087,-5821627.31323,-5820860.37694,-5836780.45965,-5868415.98829,-5840987.06814,-5838227.17932,-5843736.56489,-5852397.98775,-5856537.87728,-5857887.33366,-5869968.3054,-5899256.57528,-5985661.00795,-5927751.02996,-5956801.49746,-5915396.93489,-5903438.58908,-5900368.18448,-5923657.7528,-5971066.03108,-6058214.67597,-5996882.73853,-5952227.64174,-5956081.90559,-5958021.93473,-6019441.59966,-6202372.07811,-6030838.97887,-5995013.54658,-5986821.9308,-5986523.53408,-6003380.78559,-6018426.01063,-6014029.63473,-6050617.29614,-6174214.01679,-6175627.87735,-6271705.74182,-6657225.50566,-7802907.37546,-9925062.27734,-15631485.3766,-21637300.6575,-32255189.5224,-47929250.3074,-63007803.9699,-72942365.2584,-81113258.8592,-87994540.9417,-94323588.0452,-100411902.799,-106347449.951,-112258287.406,-118190275.792,-124185038.048,-130274367.97,-136457417.992,-142720313.132,-149043914.067,-155403433.988,-161767751.951,-168100565.584,-174362689.812,-180514544.613,-186518257.197,-192339144.137,-197946557.031,-203314195.973,-208420043.01,-213246066.64,-217777821.856,-222004036.177,-225916241.135,-229508484.983,-232777146.221,-235720857.447,-238340542.772,-240639567.369,-242623992.387,-244302919.818,-245688897.123,-246798326.804,-247651788.05,-248274123.31,-248694073.748,-248943176.096,-249053590.406};

    double dz = 0.1;
    int ref_water_table_k = 53;

    double water_table = plist_.get<double>("water table height");

    Teuchos::RCP<CompositeVector> temp = S->GetFieldData("temperature", energy_pk->name());
    Epetra_MultiVector& temp_c = *temp->ViewComponent("cell",false);
    Teuchos::RCP<CompositeVector> pres = S->GetFieldData("pressure", flow_pk->name());
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
}


void MPCFrozenCoupledFlowEnergy::fun(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
         Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> g) {
  MPCCoupledFlowEnergy::fun(t_old, t_new, u_old, u_new, g);
  g->Norm2(&the_res_norm_);
}

bool MPCFrozenCoupledFlowEnergy::is_admissible(Teuchos::RCP<const TreeVector> up) {
  if (!MPCCoupledFlowEnergy::is_admissible(up)) {
    return false;
  }

  if (backtracking_) {
    double res_prev = the_res_norm_;
    Teuchos::RCP<TreeVector> up_nc = Teuchos::rcp_const_cast<TreeVector>(up);
    Teuchos::RCP<TreeVector> g = Teuchos::rcp(new TreeVector(*up));
    fun(S_inter_->time(), S_next_->time(), Teuchos::null, up_nc, g);
    if (the_res_norm_ <= res_prev) {
      return true;
    } else {
      return false;
    }
  }

  return true;
}

} // namespace
