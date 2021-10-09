/*
  Flow PK 
 
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Neil Carlson (nnc@lanl.gov), 
           Konstantin Lipnikov (lipnikov@lanl.gov)
           Daniil Svyatskiy (dasvyat@lanl.gov)

  Interface to the BDF1 time integrator.  
*/

#include <algorithm>
#include <string>
#include <vector>

#include "Key.hh"
#include "Mesh_Algorithms.hh"
#include "CommonDefs.hh"

#include "Richards_PK.hh"

namespace Amanzi {
namespace Flow {

/* ******************************************************************
* Calculate f(u, du/dt) = a d(s(u))/dt + A*u - rhs.
****************************************************************** */
void Richards_PK::FunctionalResidual(
    double t_old, double t_new, 
    Teuchos::RCP<TreeVector> u_old, Teuchos::RCP<TreeVector> u_new, 
    Teuchos::RCP<TreeVector> f)
{ 
  double dtp(t_new - t_old);

  std::vector<int>& bc_model = op_bc_->bc_model();

  if (S_->HasFieldEvaluator(viscosity_liquid_key_)) {
    S_->GetFieldEvaluator(viscosity_liquid_key_)->HasFieldChanged(S_.ptr(), "flow");
  }
  Teuchos::RCP<const CompositeVector> mu = S_->GetFieldData(viscosity_liquid_key_);

  // compute BCs and source terms, update primary field
  UpdateSourceBoundaryData(t_old, t_new, *u_new->Data());

  // upwind diffusion coefficient and its derivative
  darcy_flux_copy->ScatterMasterToGhosted("face");

  pressure_eval_->SetFieldAsChanged(S_.ptr());
  auto alpha = S_->GetFieldData(alpha_key_, alpha_key_);
  S_->GetFieldEvaluator(alpha_key_)->HasFieldChanged(S_.ptr(), "flow");
  
  if (!flow_on_manifold_) {
    *alpha_upwind_->ViewComponent("cell") = *alpha->ViewComponent("cell");
    Operators::BoundaryFacesToFaces(bc_model, *alpha, *alpha_upwind_);
    upwind_->Compute(*darcy_flux_copy, *u_new->Data(), bc_model, *alpha_upwind_);
  }

  // modify relative permeability coefficient for influx faces
  // UpwindInflowBoundary_New(u_new->Data());

  if (!flow_on_manifold_) {
    Key der_key = Keys::getDerivKey(alpha_key_, pressure_key_);
    S_->GetFieldEvaluator(alpha_key_)->HasFieldDerivativeChanged(S_.ptr(), passwd_, pressure_key_);
    auto alpha_dP = S_->GetFieldData(der_key);

    *alpha_upwind_dP_->ViewComponent("cell") = *alpha_dP->ViewComponent("cell");
    Operators::BoundaryFacesToFaces(bc_model, *alpha_dP, *alpha_upwind_dP_);
    upwind_->Compute(*darcy_flux_copy, *u_new->Data(), bc_model, *alpha_upwind_dP_);
  }

  // assemble residual for diffusion operator
  op_matrix_->Init();
  op_matrix_diff_->UpdateMatrices(darcy_flux_copy.ptr(), solution.ptr());
  op_matrix_diff_->ApplyBCs(true, true, true);

  Teuchos::RCP<CompositeVector> rhs = op_matrix_->rhs();
  AddSourceTerms(*rhs);

  op_matrix_->ComputeNegativeResidual(*u_new->Data(), *f->Data());

  // add accumulation term 
  Epetra_MultiVector& f_cell = *f->Data()->ViewComponent("cell");

  pressure_eval_->SetFieldAsChanged(S_.ptr());
  S_->GetFieldEvaluator(porosity_key_)->HasFieldChanged(S_.ptr(), "flow");
  const Epetra_MultiVector& phi_c = *S_->GetFieldData(porosity_key_)->ViewComponent("cell");

  S_->GetFieldEvaluator(water_content_key_)->HasFieldChanged(S_.ptr(), "flow");
  const Epetra_MultiVector& wc_c = *S_->GetFieldData(water_content_key_)->ViewComponent("cell");
  const Epetra_MultiVector& wc_prev_c = *S_->GetFieldData(prev_water_content_key_)->ViewComponent("cell");

  for (int c = 0; c < ncells_owned; ++c) {
    double wc1 = wc_c[0][c];
    double wc2 = wc_prev_c[0][c];

    double factor = mesh_->cell_volume(c) / dtp;
    f_cell[0][c] += (wc1 - wc2) * factor;
  }

  // add vapor diffusion 
  if (vapor_diffusion_) {
    Functional_AddVaporDiffusion_(f->Data());
  }

  // add water content in matrix
  if (multiscale_porosity_) {
    pressure_matrix_eval_->SetFieldAsChanged(S_.ptr());
    Functional_AddMassTransferMatrix_(dtp, f->Data());
  }

  // calculate normalized residual
  functional_max_norm = 0.0;
  functional_max_cell = 0;

  for (int c = 0; c < ncells_owned; ++c) {
    const auto& dens_c = *S_->GetFieldData(mol_density_liquid_key_)->ViewComponent("cell");
    double factor = mesh_->cell_volume(c) * dens_c[0][c] * phi_c[0][c] / dtp;
    double tmp = fabs(f_cell[0][c]) / factor;
    if (tmp > functional_max_norm) {
      functional_max_norm = tmp;
      functional_max_cell = c;        
    }
  } 
}


/* ******************************************************************
* Calculate additional conribution to Richards functional:
*  f += -div (phi s_g tau_g n_g D_g \grad X_g), where
*    s_g   - gas saturation
*    tau_g - gas tortuosity
*    n_g   - molar density of gas
*    D_g   - diffusion coefficient
*    X_g   - the molar fraction of water in gas (vapor) phase
*
* Accumulation term due to water vapor is included in the volumetric
* water content field.
****************************************************************** */
void Richards_PK::Functional_AddVaporDiffusion_(Teuchos::RCP<CompositeVector> f)
{
  Key temperature_key = Keys::getKey(domain_, "temperature"); 

  const CompositeVector& pres = *S_->GetFieldData(pressure_key_);
  const CompositeVector& temp = *S_->GetFieldData(temperature_key);

  // Compute conductivities
  Teuchos::RCP<CompositeVector> kvapor_pres = Teuchos::rcp(new CompositeVector(f->Map()));
  Teuchos::RCP<CompositeVector> kvapor_temp = Teuchos::rcp(new CompositeVector(f->Map()));
  CalculateVaporDiffusionTensor_(kvapor_pres, kvapor_temp);

  // Calculate vapor contribution due to temperature.
  // We assume the same DOFs for pressure and temperature. 
  // We assume that field temperature has already essential BCs.
  op_vapor_->Init();
  op_vapor_diff_->SetScalarCoefficient(kvapor_temp, Teuchos::null);
  op_vapor_diff_->UpdateMatrices(Teuchos::null, Teuchos::null);
  op_vapor_diff_->ApplyBCs(false, false, false);

  // -- Calculate residual due to temperature
  CompositeVector g(*f);
  op_vapor_->ComputeNegativeResidual(temp, g);
  f->Update(1.0, g, 1.0);

  // Calculate vapor contribution due to capillary pressure.
  // We elliminate essential BCs to re-use the local Op for PC.
  op_vapor_->Init();
  op_vapor_diff_->SetScalarCoefficient(kvapor_pres, Teuchos::null);
  op_vapor_diff_->UpdateMatrices(Teuchos::null, Teuchos::null);
  op_vapor_diff_->ApplyBCs(false, true, false);

  // -- Calculate residual due to pressure
  op_vapor_->ComputeNegativeResidual(pres, g);
  f->Update(1.0, g, 1.0);
}


/* ******************************************************************
* Calculation of diffusion coefficient for vapor diffusion operator.
****************************************************************** */
void Richards_PK::CalculateVaporDiffusionTensor_(Teuchos::RCP<CompositeVector>& kvapor_pres,
                                                 Teuchos::RCP<CompositeVector>& kvapor_temp)
{
  AMANZI_ASSERT(domain_ == "domain");
  Key temperature_key = Keys::getKey(domain_, "temperature"); 

  S_->GetFieldEvaluator("molar_density_gas")->HasFieldChanged(S_.ptr(), passwd_);
  const Epetra_MultiVector& n_g = *S_->GetFieldData("molar_density_gas")->ViewComponent("cell");

  S_->GetFieldEvaluator(porosity_key_)->HasFieldChanged(S_.ptr(), passwd_);
  const Epetra_MultiVector& phi = *S_->GetFieldData(porosity_key_)->ViewComponent("cell");

  S_->GetFieldEvaluator(saturation_liquid_key_)->HasFieldChanged(S_.ptr(), passwd_);
  const Epetra_MultiVector& s_l = *S_->GetFieldData(saturation_liquid_key_)->ViewComponent("cell");

  S_->GetFieldEvaluator("molar_density_liquid")->HasFieldChanged(S_.ptr(), passwd_);
  const Epetra_MultiVector& n_l = *S_->GetFieldData("molar_density_liquid")->ViewComponent("cell");

  S_->GetFieldEvaluator("molar_fraction_gas")->HasFieldChanged(S_.ptr(), passwd_);
  const Epetra_MultiVector& x_g = *S_->GetFieldData("molar_fraction_gas")->ViewComponent("cell");

  std::string der_name = Keys::getDerivKey("molar_fraction_gas", temperature_key);
  S_->GetFieldEvaluator("molar_fraction_gas")->HasFieldDerivativeChanged(S_.ptr(), passwd_, temperature_key);
  const Epetra_MultiVector& dxgdT = *S_->GetFieldData(der_name)->ViewComponent("cell");

  const Epetra_MultiVector& temp = *S_->GetFieldData(temperature_key)->ViewComponent("cell");
  const Epetra_MultiVector& pres = *S_->GetFieldData(pressure_key_)->ViewComponent("cell");

  Epetra_MultiVector& kp_cell = *kvapor_pres->ViewComponent("cell");
  Epetra_MultiVector& kt_cell = *kvapor_temp->ViewComponent("cell");

  double a = 4.0 / 3.0;
  double b = 10.0 / 3.0;
  double Dref = 0.282;
  double Pref = atm_pressure_;
  double Tref = 298.0;  // Kelvins
  double R = CommonDefs::IDEAL_GAS_CONSTANT_R;

  for (int c = 0; c != ncells_owned; ++c) {
    // Millington Quirk fit for tortuosity
    double tau_phi_sat_g = pow(phi[0][c], a) * pow((1.0 - s_l[0][c]), b);
    double D_g = Dref * (Pref / atm_pressure_) * pow(temp[0][c] / Tref, 1.8);
    double tmp = tau_phi_sat_g * n_g[0][c] * D_g;

    double nRT = n_l[0][c] * temp[0][c] * R;
    double pc = atm_pressure_ - pres[0][c];
    tmp *= exp(-pc / nRT);

    kp_cell[0][c] = tmp * x_g[0][c] / nRT;
    kt_cell[0][c] = tmp * (dxgdT[0][c] / atm_pressure_ 
                        +  x_g[0][c] * pc / (nRT * temp[0][c]));
  }
}


/* ******************************************************************
* Calculate additional conribution to Richards functional:
*  f += alpha (p_f - p_m), where
*    p_f   - pressure in the fracture
*    p_m   - pressure in the matrix
*    alpha - piecewise constant mass fransfer coeffiecient
****************************************************************** */
void Richards_PK::Functional_AddMassTransferMatrix_(double dt, Teuchos::RCP<CompositeVector> f)
{
  const Epetra_MultiVector& pcf = *S_->GetFieldData(pressure_key_)->ViewComponent("cell");
  const Epetra_MultiVector& pcm = *S_->GetFieldData("pressure_matrix")->ViewComponent("cell");

  S_->GetFieldEvaluator("porosity_matrix")->HasFieldChanged(S_.ptr(), "flow");
  const Epetra_MultiVector& phi = *S_->GetFieldData("porosity_matrix")->ViewComponent("cell");

  const Epetra_MultiVector& wcm_prev = *S_->GetFieldData("prev_water_content_matrix")->ViewComponent("cell");
  Epetra_MultiVector& wcm = *S_->GetFieldData("water_content_matrix", passwd_)->ViewComponent("cell");

  Epetra_MultiVector& fc = *f->ViewComponent("cell");

  double phi0, wcm0, wcm1, pcf0;
  WhetStone::DenseVector pcm0(1);
  for (int c = 0; c < ncells_owned; ++c) {
    pcf0 = atm_pressure_ - pcf[0][c];
    pcm0(0) = atm_pressure_ - pcm[0][c];
    wcm0 = wcm_prev[0][c];
    phi0 = phi[0][c];

    int max_itrs(100);
    wcm1 = msp_->second[(*msp_->first)[c]]->WaterContentMatrix(
        pcf0, pcm0, wcm0, dt, phi0, molar_rho_, max_itrs);

    fc[0][c] += (wcm1 - wcm0) / dt;
    wcm[0][c] = wcm1;
    pcm[0][c] = atm_pressure_ - pcm0(0);

    ms_itrs_ += max_itrs;
  }
  ms_calls_ += ncells_owned;
}


/* ******************************************************************
* Calculate volumetric water content in matrix
****************************************************************** */
void Richards_PK::CalculateVWContentMatrix_()
{
  const Epetra_MultiVector& pcm = *S_->GetFieldData("pressure_matrix")->ViewComponent("cell");
  S_->GetFieldEvaluator("porosity_matrix")->HasFieldChanged(S_.ptr(), "flow");
  const Epetra_MultiVector& phi = *S_->GetFieldData("porosity_matrix")->ViewComponent("cell");
  Epetra_MultiVector& wcm = *S_->GetFieldData("water_content_matrix", passwd_)->ViewComponent("cell");

  double phi0, pcm0;
  for (int c = 0; c < ncells_owned; ++c) {
    pcm0 = atm_pressure_ - pcm[0][c];
    phi0 = phi[0][c];
    wcm[0][c] = msp_->second[(*msp_->first)[c]]->ComputeField(phi0, molar_rho_, pcm0);
  }
}


/* ******************************************************************
* Apply preconditioner inv(B) * X.                                                 
****************************************************************** */
int Richards_PK::ApplyPreconditioner(Teuchos::RCP<const TreeVector> X, 
                                      Teuchos::RCP<TreeVector> Y)
{
  Y->PutScalar(0.0);
  return op_pc_solver_->ApplyInverse(*X->Data(), *Y->Data());
}


/* ******************************************************************
* Update new preconditioner on the interval (tp-dtp, tp].
****************************************************************** */
void Richards_PK::UpdatePreconditioner(double tp, Teuchos::RCP<const TreeVector> u, double dtp)
{
  double t_old = tp - dtp;

  Teuchos::RCP<const CompositeVector> mu = S_->GetFieldData(viscosity_liquid_key_);

  std::vector<int>& bc_model = op_bc_->bc_model();
  std::vector<double>& bc_value = op_bc_->bc_value();

  // update BCs and source terms
  UpdateSourceBoundaryData(t_old, tp, *u->Data());

  // -- Darcy flux
  if (upwind_frequency_ == FLOW_UPWIND_UPDATE_ITERATION) {
    op_matrix_diff_->UpdateFlux(solution.ptr(), darcy_flux_copy.ptr());
    Epetra_MultiVector& flux = *darcy_flux_copy->ViewComponent("face");
    flux.Scale(1.0 / molar_rho_);  // FIXME
  }
  darcy_flux_copy->ScatterMasterToGhosted("face");

  // diffusion coefficient and its derivative
  pressure_eval_->SetFieldAsChanged(S_.ptr());
  auto alpha = S_->GetFieldData(alpha_key_, alpha_key_);
  S_->GetFieldEvaluator(alpha_key_)->HasFieldChanged(S_.ptr(), "flow");
  
  if (!flow_on_manifold_) {
    *alpha_upwind_->ViewComponent("cell") = *alpha->ViewComponent("cell");
    Operators::BoundaryFacesToFaces(bc_model, *alpha, *alpha_upwind_);
    upwind_->Compute(*darcy_flux_copy, *u->Data(), bc_model, *alpha_upwind_);
  }

  // modify relative permeability coefficient for influx faces
  // UpwindInflowBoundary_New(u->Data());

  if (!flow_on_manifold_) {
    Key der_key = Keys::getDerivKey(alpha_key_, pressure_key_);
    S_->GetFieldEvaluator(alpha_key_)->HasFieldDerivativeChanged(S_.ptr(), passwd_, pressure_key_);
    auto alpha_dP = S_->GetFieldData(der_key);

    *alpha_upwind_dP_->ViewComponent("cell") = *alpha_dP->ViewComponent("cell");
    Operators::BoundaryFacesToFaces(bc_model, *alpha_dP, *alpha_upwind_dP_);
    upwind_->Compute(*darcy_flux_copy, *u->Data(), bc_model, *alpha_upwind_dP_);
  }

  // create diffusion operators
  op_preconditioner_->Init();
  op_preconditioner_diff_->UpdateMatrices(darcy_flux_copy.ptr(), solution.ptr());
  op_preconditioner_diff_->UpdateMatricesNewtonCorrection(darcy_flux_copy.ptr(), solution.ptr(), molar_rho_);
  op_preconditioner_diff_->ApplyBCs(true, true, true);

  // add time derivative
  if (dtp > 0.0) {
    std::string der_name = Keys::getDerivKey(water_content_key_, pressure_key_);
    S_->GetFieldEvaluator(water_content_key_)->HasFieldDerivativeChanged(S_.ptr(), passwd_, pressure_key_);
    CompositeVector& dwc_dp = *S_->GetFieldData(der_name, water_content_key_);

    op_acc_->AddAccumulationDelta(*u->Data(), dwc_dp, dwc_dp, dtp, "cell");
 
    // estimate CNLS limiters
    if (algebraic_water_content_balance_) {
      const CompositeVector& wc = *S_->GetFieldData(water_content_key_);
      CalculateCNLSLimiter_(wc, dwc_dp, bdf1_dae_->tol_solver());
    }
  }

  // Add vapor diffusion. We assume that the corresponding local operator
  // has been already populated during functional evaluation.
  if (vapor_diffusion_) {
    Teuchos::RCP<CompositeVector> kvapor_pres = Teuchos::rcp(new CompositeVector(u->Data()->Map()));
    Teuchos::RCP<CompositeVector> kvapor_temp = Teuchos::rcp(new CompositeVector(u->Data()->Map()));
    CalculateVaporDiffusionTensor_(kvapor_pres, kvapor_temp);

    op_vapor_->Init();
    op_vapor_diff_->SetScalarCoefficient(kvapor_pres, Teuchos::null);
    op_vapor_diff_->UpdateMatrices(Teuchos::null, Teuchos::null);
    op_vapor_diff_->ApplyBCs(false, true, false);
  }
 
  // finalize preconditioner
  op_pc_solver_->ComputeInverse();
}


/* ******************************************************************
* Modify preconditior as needed.
****************************************************************** */
bool Richards_PK::ModifyPredictor(double dt, Teuchos::RCP<const TreeVector> u0,
                                  Teuchos::RCP<TreeVector> u)
{
  Teuchos::RCP<TreeVector> du = Teuchos::rcp(new TreeVector(*u));
  du->Update(-1.0, *u0, 1.0);
 
  ModifyCorrection(dt, Teuchos::null, u0, du);

  *u = *u0;
  u->Update(1.0, *du, 1.0);
  return true;
}


/* ******************************************************************
* Check difference du between the predicted and converged solutions.
* This is a wrapper for various error control methods. 
****************************************************************** */
double Richards_PK::ErrorNorm(Teuchos::RCP<const TreeVector> u, 
                              Teuchos::RCP<const TreeVector> du)
{
  double error;
  error = ErrorNormSTOMP(*u->Data(), *du->Data());

  // exact algebraic relation between saturation and Darcy flux
  // requires to save the last increment.
  if (algebraic_water_content_balance_) {
    *cnls_limiter_->ViewComponent("dpre") = *du->Data()->ViewComponent("cell");
  }
 
  return error;
}


/* ******************************************************************
* Error control a-la STOMP.
****************************************************************** */
double Richards_PK::ErrorNormSTOMP(const CompositeVector& u, const CompositeVector& du)
{
  const Epetra_MultiVector& uc = *u.ViewComponent("cell");
  const Epetra_MultiVector& duc = *du.ViewComponent("cell");

  double error, error_p, error_r;
  int cell_p(0);

  if (error_control_ & FLOW_TI_ERROR_CONTROL_PRESSURE) {
    error_p = 0.0;
    for (int c = 0; c < ncells_owned; c++) {
      double tmp = fabs(duc[0][c]) / (fabs(uc[0][c] - atm_pressure_) + atm_pressure_);
      if (tmp > error_p) {
        error_p = tmp;
        cell_p = c;
      } 
    }
  } else {
    error_p = 0.0;
  }

  if (error_control_ & FLOW_TI_ERROR_CONTROL_RESIDUAL) {
    error_r = functional_max_norm;
  } else {
    error_r = 0.0;
  }

  error = std::max(error_r, error_p);

#ifdef HAVE_MPI
  double buf = error;
  du.Comm()->MaxAll(&buf, &error, 1);  // find the global maximum
#endif

  // maximum error is printed out only on one processor
  if (vo_->getVerbLevel() >= Teuchos::VERB_EXTREME) {
    if (error == buf) {
      int c = functional_max_cell;
      const AmanziGeometry::Point& xp = mesh_->cell_centroid(c);

      Teuchos::OSTab tab = vo_->getOSTab();
      *vo_->os() << "residual=" << functional_max_norm << " at point";
      for (int i = 0; i < dim; i++) *vo_->os() << " " << xp[i];
      *vo_->os() << std::endl;
 
      c = cell_p;
      const AmanziGeometry::Point& yp = mesh_->cell_centroid(c);

      *vo_->os() << "pressure err=" << error_p << " at point";
      for (int i = 0; i < dim; i++) *vo_->os() << " " << yp[i];
      *vo_->os() << std::endl;

      double s = wrm_->second[(*wrm_->first)[c]]->saturation(atm_pressure_ - uc[0][c]);
      *vo_->os() << "saturation=" << s << " pressure=" << uc[0][c] << std::endl;
    }
  }

  return error;
}


/********************************************************************
* Modifies nonlinear update du based on the maximum allowed change
* of saturation.
****************************************************************** */
AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
    Richards_PK::ModifyCorrection(double dt, Teuchos::RCP<const TreeVector> f,
                                  Teuchos::RCP<const TreeVector> u,
                                  Teuchos::RCP<TreeVector> du)
{
  const Epetra_MultiVector& uc = *u->Data()->ViewComponent("cell");
  const Epetra_MultiVector& duc = *du->Data()->ViewComponent("cell");

  AmanziGeometry::Point face_centr, cell_cntr;
  double max_sat_pert(0.25), damping_factor(0.5);
  
  if (fp_list_->isSublist("clipping parameters")) {
    Teuchos::ParameterList& clip_list = fp_list_->sublist("clipping parameters");
    max_sat_pert = clip_list.get<double>("maximum saturation change", 0.25);
    damping_factor = clip_list.get<double>("pressure damping factor", 0.5);
  }

  int nsat_clipped(0), npre_clipped(0);

  for (int c = 0; c < ncells_owned; c++) {
    double pc = atm_pressure_ - uc[0][c];
    double sat = wrm_->second[(*wrm_->first)[c]]->saturation(pc);
    double sat_pert;
    if (sat >= 0.5) sat_pert = sat - max_sat_pert;
    else sat_pert = sat + max_sat_pert;
    
    double press_pert = atm_pressure_ - wrm_->second[(*wrm_->first)[c]]->capillaryPressure(sat_pert);
    double du_pert_max = fabs(uc[0][c] - press_pert); 

    if ((fabs(duc[0][c]) > du_pert_max) && (1 - sat > 1e-5)) {
      // std::cout << "clip saturation: c=" << c << " p=" << uc[0][c]
      //           << " dp: " << duc[0][c] << " -> " << du_pert_max << std::endl;

      if (duc[0][c] >= 0.0) duc[0][c] = du_pert_max;
      else duc[0][c] = -du_pert_max;
      
      nsat_clipped++;
    }    
  }

  for (int c = 0; c < ncells_owned; c++) {
    double unew = uc[0][c] - duc[0][c];
    double tmp = duc[0][c];

    if ((unew > atm_pressure_) && (uc[0][c] < atm_pressure_)) {
      // std::cout << "pressure change: " << uc[0][c] << " -> " << unew << std::endl;
      duc[0][c] = tmp * damping_factor;
      npre_clipped++;
    }
  }

  // output statistics
  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
    int nsat_tmp = nsat_clipped, npre_tmp = npre_clipped;
    mesh_->get_comm()->SumAll(&nsat_tmp, &nsat_clipped, 1);
    mesh_->get_comm()->SumAll(&npre_tmp, &npre_clipped, 1);

    if (nsat_clipped > 0 || npre_clipped > 0) {
      Teuchos::OSTab tab = vo_->getOSTab();
      *vo_->os() << vo_->color("green") << "saturation/pressure clipped in " 
                 << nsat_clipped << "/" << npre_clipped << " cells" << vo_->reset() << std::endl;
    }
  }

  return (nsat_clipped + npre_clipped) > 0 ? AmanziSolvers::FnBaseDefs::CORRECTION_MODIFIED :
      AmanziSolvers::FnBaseDefs::CORRECTION_NOT_MODIFIED;
}

}  // namespace Flow
}  // namespace Amanzi



