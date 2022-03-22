/* -*-  mode++; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Svetlana Tokareva

Solves:

c rho dT/dt + q dot grad h = div Ke grad T + S
------------------------------------------------------------------------- */

#include "advection.hh"
#include "FieldEvaluator.hh"
#include "lake_thermo_pk.hh"
#include "Op.hh"
#include "pk_helpers.hh"

namespace Amanzi {
namespace LakeThermo {

// -------------------------------------------------------------
// Accumulation of energy term c rho dT/dt
// -------------------------------------------------------------
void Lake_Thermo_PK::AddAccumulation_(const Teuchos::Ptr<CompositeVector>& g) {
  double dt = S_next_->time() - S_inter_->time();

//  // update the temperature at both the old and new times.
//  S_next_->GetFieldEvaluator(temperature_key_)->HasFieldChanged(S_next_.ptr(), name_);
//  S_inter_->GetFieldEvaluator(temperature_key_)->HasFieldChanged(S_inter_.ptr(), name_);

  // get the energy at each time
  Teuchos::RCP<const CompositeVector> T1 = S_next_->GetFieldData(temperature_key_);
  Teuchos::RCP<const CompositeVector> T0 = S_inter_->GetFieldData(temperature_key_);

  S_inter_->GetFieldEvaluator(density_key_)->HasFieldChanged(S_inter_.ptr(), name_);

  // evaluate density
  const Epetra_MultiVector& rho =
      *S_inter_->GetFieldData(density_key_)->ViewComponent("cell",false);

  // Update the residual with the accumulation of energy over the
  // timestep, on cells.
  g->ViewComponent("cell", false)
        ->Update(cp_*rho0/dt, *T1->ViewComponent("cell", false),
            -cp_*rho0/dt, *T0->ViewComponent("cell", false), 1.0);


  const Epetra_MultiVector& T1_c = *S_next_->GetFieldData(temperature_key_)->ViewComponent("cell", false);
  const Epetra_MultiVector& T0_c = *S_inter_->GetFieldData(temperature_key_)->ViewComponent("cell", false);

  const Epetra_MultiVector& g_c = *g->ViewComponent("cell", false);

  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

};


// -------------------------------------------------------------
// Advective term for transport of enthalpy, q dot grad h.
// -------------------------------------------------------------
void Lake_Thermo_PK::AddAdvection_(const Teuchos::Ptr<State>& S,
    const Teuchos::Ptr<CompositeVector>& g, bool negate) {

  // set up the operator

  // NOTE: fluxes are a MOLAR flux by choice of the flow pk, i.e.
  // [flux] =  mol/s

  // NOTE: this will be the eventual way to ensure it is up to date,
  // but there is no FieldEvaluator for darcy flux yet.  When there
  // is, we can take the evaluation out of Flow::commit_state(),
  // but for now we'll leave it there and assume it has been updated. --etc
  //  S->GetFieldEvaluator(flux_key_)->HasFieldChanged(S.ptr(), name_);
  Teuchos::RCP<const CompositeVector> flux = S->GetFieldData(flux_key_);

  S_inter_->GetFieldEvaluator(density_key_)->HasFieldChanged(S_inter_.ptr(), name_);

  // evaluate density
  const Epetra_MultiVector& rho =
      *S_inter_->GetFieldData(density_key_)->ViewComponent("cell",false);

  Teuchos::ParameterList& param_list = plist_->sublist("parameters");
  FunctionFactory fac;
  Teuchos::RCP<Function> r_func_ = Teuchos::rcp(fac.Create(param_list.sublist("precipitation")));
  Teuchos::RCP<Function> E_func_ = Teuchos::rcp(fac.Create(param_list.sublist("evaporation")));

  std::vector<double> args(1);
  args[0] = S->time();
  r_ = (*r_func_)(args);
  E_ = (*E_func_)(args);

  // Compute evaporartion rate
  Teuchos::ParameterList& param_list_eval = plist_->sublist("surface flux evaluator");
  Teuchos::ParameterList& param_list_param = param_list_eval.sublist("parameters");
  // read parameters from the met data
  Teuchos::RCP<Amanzi::Function> SS_func_ = Teuchos::rcp(fac.Create(param_list_param.sublist("solar radiation")));
  Teuchos::RCP<Amanzi::Function> E_a_func_ = Teuchos::rcp(fac.Create(param_list_param.sublist("atmospheric downward radiation")));
  Teuchos::RCP<Amanzi::Function> T_a_func_ = Teuchos::rcp(fac.Create(param_list_param.sublist("air temperature")));
  Teuchos::RCP<Amanzi::Function> H_a_func_ = Teuchos::rcp(fac.Create(param_list_param.sublist("air humidity")));
  Teuchos::RCP<Amanzi::Function> P_a_func_ = Teuchos::rcp(fac.Create(param_list_param.sublist("atmospheric pressure")));

  args[0] = S->time();
  double SS = (*SS_func_)(args);
  double E_a = (*E_a_func_)(args);
  double T_a = (*T_a_func_)(args);
  double q_a = (*H_a_func_)(args);
  double P_a = (*P_a_func_)(args);

  // get temperature
  const Epetra_MultiVector& temp_v = *S->GetFieldData(temperature_key_)
		    ->ViewComponent("cell",false);

  Epetra_MultiVector& g_c = *g->ViewComponent("cell",false);
  unsigned int ncomp = g_c.MyLength();

  double T_s = temp_v[0][ncomp-1];

  double b1_vap   = 610.78;        // Coefficient [N m^{-2} = kg m^{-1} s^{-2}]
  double b3_vap   = 273.16;        // Triple point [K]
  double b2w_vap  = 17.2693882;    // Coefficient (water)
  double b2i_vap  = 21.8745584;    // Coefficient (ice)
  double b4w_vap  = 35.86;         // Coefficient (temperature) [K]
  double b4i_vap  = 7.66;          // Coefficient (temperature) [K]

  // Saturation water vapour pressure [N m^{-2} = kg m^{-1} s^{-2}]
  double wvpres_s;

  double h_ice = 0.;
  double h_Ice_min_flk = 1.e-9;
  if (h_ice < h_Ice_min_flk) {  // Water surface
    wvpres_s = b1_vap*std::exp(b2w_vap*(T_s-b3_vap)/(T_s-b4w_vap));
  } else {                      // Ice surface
    wvpres_s = b1_vap*std::exp(b2i_vap*(T_s-b3_vap)/(T_s-b4i_vap));
  }

  // Saturation specific humidity at T=T_s
  double tpsf_R_dryair    = 2.8705e2;  // Gas constant for dry air [J kg^{-1} K^{-1}]
  double tpsf_R_watvap    = 4.6151e2;  // Gas constant for water vapour [J kg^{-1} K^{-1}]
  double tpsf_Rd_o_Rv  = tpsf_R_dryair/tpsf_R_watvap;
  double q_s = tpsf_Rd_o_Rv*wvpres_s/(P_a-(1.-tpsf_Rd_o_Rv)*wvpres_s);

  double height_tq = 2;
  double tpsf_kappa_q_a   = 2.4e-05; // Molecular diffusivity of air for water vapour [m^{2} s^{-1}]

  double LE = -tpsf_kappa_q_a*(q_a-q_s)/height_tq;

  double rho_a = P_a/tpsf_R_dryair/T_s/(1.+(1./tpsf_Rd_o_Rv-1.)*q_s);

  double tpsf_c_a_p  = 1.005e3; // Specific heat of air at constant pressure [J kg^{-1} K^{-1}]
  double tpsf_L_evap = 2.501e6; // Specific heat of evaporation [J kg^{-1}]
  double tpl_L_f     = 3.3e5;   // Latent heat of fusion [J kg^{-1}]

  double Q_watvap   = LE*rho_a;
  LE = tpsf_L_evap;
  if (h_ice >= h_Ice_min_flk) LE = LE + tpl_L_f;   // Add latent heat of fusion over ice
  LE = Q_watvap*LE;

  double row0 = 1.e+3;
  double evap_rate = LE/(row0*tpsf_L_evap); //*1000.;
  evap_rate = abs(evap_rate);

  E_ = evap_rate;

  // <-- end compute evaporation rate

  std::ofstream outfile;

  outfile.open("Precipitation.txt", std::ios_base::app); // append instead of overwrite
  outfile.precision(5);
  outfile.setf(std::ios::scientific,std::ios::floatfield);
  outfile << S->time() << " " << r_ << "\n";

  double dhdt = r_ - E_ - R_s_ - R_b_;
  double B_w  = r_ - E_;

  const Epetra_MultiVector& flux_f = *flux->ViewComponent("face", false);

  int nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);

  double dt = S_next_->time() - S_inter_->time();

  for (int f = 0; f < nfaces_owned; f++) {
    const AmanziGeometry::Point& xcf = mesh_->face_centroid(f);

    AmanziGeometry::Point normal = mesh_->face_normal(f);
    normal /= norm(normal);

    flux_f[0][f] = -1.*cp_*rho0*(dhdt*(1.-xcf[2]) - B_w)/(h_+1.e-6); // *normal or not?
  }

  db_->WriteVector(" adv flux", flux.ptr(), true);
  matrix_adv_->global_operator()->Init();
  matrix_adv_->Setup(*flux);
  matrix_adv_->SetBCs(bc_adv_, bc_adv_);
  matrix_adv_->UpdateMatrices(flux.ptr());

  // apply to temperature
  Teuchos::RCP<const CompositeVector> temp = S->GetFieldData(temperature_key_);
  ApplyDirichletBCsToTemperature_(S.ptr());
  matrix_adv_->ApplyBCs(false, true, false);

  // apply
  matrix_adv_->global_operator()->ComputeNegativeResidual(*temp, *g, false);
}

// -------------------------------------------------------------
// Diffusion term, div K grad T
// -------------------------------------------------------------
void Lake_Thermo_PK::ApplyDiffusion_(const Teuchos::Ptr<State>& S,
    const Teuchos::Ptr<CompositeVector>& g) {
  // update the thermal conductivity
  UpdateConductivityData_(S_next_.ptr());
  Teuchos::RCP<const CompositeVector> conductivity =
      S_next_->GetFieldData(uw_conductivity_key_);

  const Epetra_MultiVector& uw_cond_c = *S_next_->GetFieldData(uw_conductivity_key_)->ViewComponent("face", false);
  int nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);

  // because the lake model has 1/h^2 term
  for (int f = 0; f < nfaces_owned; f++) {
    uw_cond_c[0][f] /= (h_*h_);
  }

  Teuchos::RCP<const CompositeVector> temp = S->GetFieldData(key_);

  // update the stiffness matrix
  matrix_diff_->global_operator()->Init();
  matrix_diff_->SetScalarCoefficient(conductivity, Teuchos::null);
  matrix_diff_->UpdateMatrices(Teuchos::null, temp.ptr());
  matrix_diff_->ApplyBCs(true, true, true);

  // calculate the residual
  matrix_diff_->global_operator()->ComputeNegativeResidual(*temp, *g);
};


// ---------------------------------------------------------------------
// Add in energy source, which are accumulated by a single evaluator.
// ---------------------------------------------------------------------
void Lake_Thermo_PK::AddSources_(const Teuchos::Ptr<State>& S,
    const Teuchos::Ptr<CompositeVector>& g) {
  Teuchos::OSTab tab = vo_->getOSTab();

  is_source_term_ = true;

  // external sources of energy
  if (is_source_term_) {
    Epetra_MultiVector& g_c = *g->ViewComponent("cell",false);

    // Update the source term
    //    S->GetFieldEvaluator(source_key_)->HasFieldChanged(S, name_);
    //    const Epetra_MultiVector& source1 =
    //        *S->GetFieldData(source_key_)->ViewComponent("cell",false);
    const Epetra_MultiVector& cv =
        *S->GetFieldData(Keys::getKey(domain_,"cell_volume"))->ViewComponent("cell",false);

    Teuchos::ParameterList& param_list = plist_->sublist("parameters");
    FunctionFactory fac;
    Teuchos::RCP<Function> r_func_ = Teuchos::rcp(fac.Create(param_list.sublist("precipitation")));
    Teuchos::RCP<Function> E_func_ = Teuchos::rcp(fac.Create(param_list.sublist("evaporation")));

    std::vector<double> args(1);
    args[0] = S->time();
    r_ = (*r_func_)(args);
    E_ = (*E_func_)(args);


    // Compute evaporartion rate
    Teuchos::ParameterList& param_list_eval = plist_->sublist("surface flux evaluator");
    Teuchos::ParameterList& param_list_param = param_list_eval.sublist("parameters");
    // read parameters from the met data
    Teuchos::RCP<Amanzi::Function> SS_func_ = Teuchos::rcp(fac.Create(param_list_param.sublist("solar radiation")));
    Teuchos::RCP<Amanzi::Function> E_a_func_ = Teuchos::rcp(fac.Create(param_list_param.sublist("atmospheric downward radiation")));
    Teuchos::RCP<Amanzi::Function> T_a_func_ = Teuchos::rcp(fac.Create(param_list_param.sublist("air temperature")));
    Teuchos::RCP<Amanzi::Function> H_a_func_ = Teuchos::rcp(fac.Create(param_list_param.sublist("air humidity")));
    Teuchos::RCP<Amanzi::Function> P_a_func_ = Teuchos::rcp(fac.Create(param_list_param.sublist("atmospheric pressure")));

    args[0] = S->time();
    double SS = (*SS_func_)(args);
    double E_a = (*E_a_func_)(args);
    double T_a = (*T_a_func_)(args);
    double q_a = (*H_a_func_)(args);
    double P_a = (*P_a_func_)(args);

    // get temperature
    const Epetra_MultiVector& temp_v = *S->GetFieldData(temperature_key_)
	      ->ViewComponent("cell",false);
    unsigned int ncomp = g_c.MyLength();

    double T_s = temp_v[0][ncomp-1];

    double b1_vap   = 610.78;        // Coefficient [N m^{-2} = kg m^{-1} s^{-2}]
    double b3_vap   = 273.16;        // Triple point [K]
    double b2w_vap  = 17.2693882;    // Coefficient (water)
    double b2i_vap  = 21.8745584;    // Coefficient (ice)
    double b4w_vap  = 35.86;         // Coefficient (temperature) [K]
    double b4i_vap  = 7.66;          // Coefficient (temperature) [K]

    // Saturation water vapour pressure [N m^{-2} = kg m^{-1} s^{-2}]
    double wvpres_s;

    double h_ice = 0.;
    double h_Ice_min_flk = 1.e-9;
    if (h_ice < h_Ice_min_flk) {  // Water surface
      wvpres_s = b1_vap*std::exp(b2w_vap*(T_s-b3_vap)/(T_s-b4w_vap));
    } else {                      // Ice surface
      wvpres_s = b1_vap*std::exp(b2i_vap*(T_s-b3_vap)/(T_s-b4i_vap));
    }

    // Saturation specific humidity at T=T_s
    double tpsf_R_dryair    = 2.8705e2;  // Gas constant for dry air [J kg^{-1} K^{-1}]
    double tpsf_R_watvap    = 4.6151e2;  // Gas constant for water vapour [J kg^{-1} K^{-1}]
    double tpsf_Rd_o_Rv  = tpsf_R_dryair/tpsf_R_watvap;
    double q_s = tpsf_Rd_o_Rv*wvpres_s/(P_a-(1.-tpsf_Rd_o_Rv)*wvpres_s);

    double height_tq = 2;
    double tpsf_kappa_q_a   = 2.4e-05; // Molecular diffusivity of air for water vapour [m^{2} s^{-1}]

    double LE = -tpsf_kappa_q_a*(q_a-q_s)/height_tq;

    double rho_a = P_a/tpsf_R_dryair/T_s/(1.+(1./tpsf_Rd_o_Rv-1.)*q_s);

    double tpsf_c_a_p  = 1.005e3; // Specific heat of air at constant pressure [J kg^{-1} K^{-1}]
    double tpsf_L_evap = 2.501e6; // Specific heat of evaporation [J kg^{-1}]
    double tpl_L_f     = 3.3e5;   // Latent heat of fusion [J kg^{-1}]

    double Q_watvap   = LE*rho_a;
    LE = tpsf_L_evap;
    if (h_ice >= h_Ice_min_flk) LE = LE + tpl_L_f;   // Add latent heat of fusion over ice
    LE = Q_watvap*LE;

    double row0 = 1.e+3;
    double evap_rate = LE/(row0*tpsf_L_evap); //*1000.;
    evap_rate = abs(evap_rate);

    E_ = evap_rate;

    // <-- end compute evaporation rate

    double dhdt = r_ - E_ - R_s_ - R_b_;
    double B_w  = r_ - E_;

    S->GetFieldEvaluator(density_key_)->HasFieldChanged(S.ptr(), name_);

    // evaluate density
    const Epetra_MultiVector& rho =
        *S->GetFieldData(density_key_)->ViewComponent("cell",false);

    // get temperature
    const Epetra_MultiVector& temp = *S->GetFieldData(temperature_key_)
              ->ViewComponent("cell",false);

    // get conductivity
    S->GetFieldEvaluator(conductivity_key_)->HasFieldChanged(S.ptr(), name_);
    const Epetra_MultiVector& lambda_c =
        *S->GetFieldData(conductivity_key_)->ViewComponent("cell",false);

    // water extinction coefficient
    alpha_e_ = 1.04;

    // Add into residual
    unsigned int ncells = g_c.MyLength();
    for (unsigned int c=0; c!=ncells; ++c) {
      const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);

      if (temp[0][ncells-1] < 273.15) {
        S0_ = 0.;
      } else {
        S0_ = 0.; //0.1*SS;
      }

      // -1.* because I switched to vertical xi coordinate
      g_c[0][c] += -1.*( (S0_*exp(-alpha_e_*h_*(1.-xc[2]))*(alpha_e_*h_)/(h_+1.e-6) - cp_*rho0*temp[0][c]*dhdt/(h_+1.e-6)) * cv[0][c] );

      /* TESTING
//      // manufactured solution 1: linear temperature distribution
//      double T0 = 280., T1 = 400.;
//      g_c[0][c] += cp_*rho[0][c]*(T1-T0)/(h_+1.e-6)*(dhdt*xc[2] - B_w) * cv[0][c]; // - cp_*rho[0][c]*temp[0][c]*dhdt/(h_+1.e-6) * cv[0][c];
      // manufactured solution 2: exponential temperature distribution
      double T0 = 280., T1 = 400.;
      double coef = log(T1/T0);
//      std::cout << "lambda_c[0][c] = " << lambda_c[0][c] << std::endl;
//      std::cout << "coef = " << coef << std::endl;
//      std::cout << "exp(coef*xc[2]) = " << exp(coef*xc[2]) << std::endl;
//      std::cout << "h_ = " << h_ << std::endl;
      g_c[0][c] += T0*coef*exp(coef*xc[2])*( 1./(h_+1.e-6)*lambda_c[0][c]*coef - cp_*rho[0][c]*(dhdt*xc[2] - B_w) )/(h_+1.e-6) * cv[0][c];
//          - cp_*rho[0][c]*temp[0][c]*dhdt/(h_+1.e-6) * cv[0][c];
//      // manufactured solution 3: quadratic temperature distribution
//      double T0 = 280., T1 = 400.;
//      double zm = 0.5, Tm = 300.;
//      double d = T0;
//      double b = (Tm-T0-zm*zm*(T1-T0))/(zm*(1.-zm));
//      double a = T1-T0-b;
////      std::cout << "a = " << a << ", b = " << b << std::endl;
//      g_c[0][c] += (2.*a*lambda_c[0][c]/(h_+1.e-6) - cp_*rho[0][c]*(2.*a*xc[2] + b)*(dhdt*xc[2] - B_w))/(h_+1.e-6) * cv[0][c];
//          - cp_*rho[0][c]*temp[0][c]*dhdt/(h_+1.e-6) * cv[0][c];
//      g_c[0][c] += 2.*a*lambda_c[0][c]/(h_+1.e-6)/(h_+1.e-6) * cv[0][c];
      //          - cp_*rho[0][c]*temp[0][c]*dhdt/(h_+1.e-6) * cv[0][c];
//      // manufactured solution without non-conservative part
//      g_c[0][c] -= (-2.*a*lambda_c[0][c]/(h_+1.e-6) -
//          cp_*rho[0][c]/(h_+1.e-6)*( dhdt*(3.*a*xc[2]*xc[2] + 2.*b*xc[2] + c) - B_w*(2.*a*xc[2] + b) ) ) * cv[0][c];
//          - cp_*rho[0][c]*temp[0][c]*dhdt/(h_+1.e-6) * cv[0][c];
//      g_c[0][c] += ( 2.*a*lambda_c[0][c]/(h_+1.e-6)/(h_+1.e-6) - (2.*a*xc[2] + b) ) * cv[0][c];
//       g_c[0][c] += ( 2.*a*lambda_c[0][c]/(h_+1.e-6)/(h_+1.e-6) - cp_*rho[0][c]*dhdt/(h_+1.e-6)*(3.*a*xc[2]*xc[2] + 2.*b*xc[2] + c) + cp_*rho[0][c]*B_w/(h_+1.e-6)*(2.*a*xc[2] + b) ) * cv[0][c]
//                 + cp_*rho[0][c]*dhdt/(h_+1.e-6)*(a*xc[2]*xc[2] + b*xc[2] + c) * cv[0][c];
//      std::cout << "h_ = " << h_ << std::endl;
//      std::cout << "SRC = " << g_c[0][c] << std::endl;
       *
       */
    }

    if (vo_->os_OK(Teuchos::VERB_EXTREME))
      *vo_->os() << "Adding external source term" << std::endl;
    //    db_->WriteVector("  Q_ext", S->GetFieldData(source_key_).ptr(), false);
    db_->WriteVector("res (src)", g, false);
  }
}


void Lake_Thermo_PK::AddSourcesToPrecon_(const Teuchos::Ptr<State>& S, double h) {
//  // external sources of energy (temperature dependent source)
//  if (is_source_term_ && is_source_term_differentiable_ &&
//      S->GetFieldEvaluator(source_key_)->IsDependency(S, key_)) {
//
//    Teuchos::RCP<const CompositeVector> dsource_dT;
//    if (!is_source_term_finite_differentiable_) {
//      // evaluate the derivative through the dag
//      S->GetFieldEvaluator(source_key_)->HasFieldDerivativeChanged(S, name_, key_);
//      dsource_dT = S->GetFieldData(Keys::getDerivKey(source_key_, key_));
//    } else {
//      // evaluate the derivative through finite differences
//      double eps = 1.e-8;
//      S->GetFieldData(key_, name_)->Shift(eps);
//      ChangedSolution();
//      S->GetFieldEvaluator(source_key_)->HasFieldChanged(S, name_);
//      auto dsource_dT_nc = Teuchos::rcp(new CompositeVector(*S->GetFieldData(source_key_)));
//
//      S->GetFieldData(key_, name_)->Shift(-eps);
//      ChangedSolution();
//      S->GetFieldEvaluator(source_key_)->HasFieldChanged(S, name_);
//
//      dsource_dT_nc->Update(-1/eps, *S->GetFieldData(source_key_), 1/eps);
//      dsource_dT = dsource_dT_nc;
//    }
//    db_->WriteVector("  dQ_ext/dT", dsource_dT.ptr(), false);
//    preconditioner_acc_->AddAccumulationTerm(*dsource_dT, -1.0, "cell", true);
//  }
}

// -------------------------------------------------------------
// Plug enthalpy into the boundary faces manually.
// This will be removed once boundary faces exist.
// -------------------------------------------------------------
void Lake_Thermo_PK::ApplyDirichletBCsToTemperature_(const Teuchos::Ptr<State>& S) {
  // put the boundary fluxes in faces for Dirichlet BCs.
  //  S->GetFieldEvaluator(enthalpy_key_)->HasFieldChanged(S, name_);

  const Epetra_MultiVector& temp_bf =
      *S->GetFieldData(temperature_key_)->ViewComponent("boundary_face",false);
  const Epetra_Map& vandalay_map = mesh_->exterior_face_map(false);
  const Epetra_Map& face_map = mesh_->face_map(false);

  int nbfaces = temp_bf.MyLength();
  for (int bf=0; bf!=nbfaces; ++bf) {
    AmanziMesh::Entity_ID f = face_map.LID(vandalay_map.GID(bf));

    if (bc_adv_->bc_model()[f] == Operators::OPERATOR_BC_DIRICHLET) {
      bc_adv_->bc_value()[f] = temp_bf[0][bf];
    }
  }
}



} //namespace LakeThermo
} //namespace Amanzi
