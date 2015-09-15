/*
  This is the flow component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Neil Carlson (nnc@lanl.gov), 
           Konstantin Lipnikov (lipnikov@lanl.gov)
           Daniil Svyatskiy (dasvyat@lanl.gov)

  The routine implements interface to the BDF1 time integrator.  
*/

#include <algorithm>
#include <string>
#include <vector>

#include "CommonDefs.hh"
#include "FieldMaps.hh"
#include "LinearOperatorFactory.hh"

#include "Richards_PK.hh"

namespace Amanzi {
namespace Flow {

/* ******************************************************************
* Calculate f(u, du/dt) = a d(s(u))/dt + A*u - rhs.
****************************************************************** */
void Richards_PK::Functional(double t_old, double t_new, 
                             Teuchos::RCP<TreeVector> u_old, Teuchos::RCP<TreeVector> u_new, 
                             Teuchos::RCP<TreeVector> f)
{ 
  double dtp(t_new - t_old);

  if (S_->HasFieldEvaluator("viscosity_liquid")) {
    S_->GetFieldEvaluator("viscosity_liquid")->HasFieldChanged(S_.ptr(), "flow");
  }
  Teuchos::RCP<const CompositeVector> mu = S_->GetFieldData("viscosity_liquid");

  // update coefficients
  darcy_flux_copy->ScatterMasterToGhosted("face");

  relperm_->Compute(u_new->Data(), krel_); 
  RelPermUpwindFn func1 = &RelPerm::Compute;
  upwind_->Compute(*darcy_flux_upwind, *u_new->Data(), bc_model, bc_value, *krel_, *krel_, func1);
  Operators::CellToFace_ScaleInverse(mu, krel_);
  krel_->ScaleMasterAndGhosted(molar_rho_);

  // modify relative permeability coefficient for influx faces
  // UpwindInflowBoundary_New(u_new->Data());

  relperm_->ComputeDerivative(u_new->Data(), dKdP_); 
  RelPermUpwindFn func2 = &RelPerm::ComputeDerivative;
  upwind_->Compute(*darcy_flux_upwind, *u_new->Data(), bc_model, bc_value, *dKdP_, *dKdP_, func2);
  Operators::CellToFace_ScaleInverse(mu, dKdP_);
  dKdP_->ScaleMasterAndGhosted(molar_rho_);

  UpdateSourceBoundaryData(t_old, t_new, *u_new->Data());
  
  // assemble residual for diffusion operator
  op_matrix_->Init();
  op_matrix_diff_->UpdateMatrices(darcy_flux_copy.ptr(), solution.ptr());
  op_matrix_diff_->ApplyBCs(true, true);

  Teuchos::RCP<CompositeVector> rhs = op_matrix_->rhs();
  if (src_sink != NULL) AddSourceTerms(*rhs);

  op_matrix_->ComputeNegativeResidual(*u_new->Data(), *f->Data());

  // add accumulation term 
  Epetra_MultiVector& f_cell = *f->Data()->ViewComponent("cell");

  pressure_eval_->SetFieldAsChanged(S_.ptr());
  S_->GetFieldEvaluator("porosity")->HasFieldChanged(S_.ptr(), "flow");
  const Epetra_MultiVector& phi_c = *S_->GetFieldData("porosity")->ViewComponent("cell");

  S_->GetFieldEvaluator("water_content")->HasFieldChanged(S_.ptr(), "flow");
  const Epetra_MultiVector& wc_c = *S_->GetFieldData("water_content")->ViewComponent("cell");
  const Epetra_MultiVector& wc_prev_c = *S_->GetFieldData("prev_water_content")->ViewComponent("cell");

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
    Functional_AddWaterContentMatrix_(f->Data());
  }

  // calculate normalized residual
  functional_max_norm = 0.0;
  functional_max_cell = 0;

  for (int c = 0; c < ncells_owned; ++c) {
    double factor = mesh_->cell_volume(c) * molar_rho_ * phi_c[0][c] / dtp;
    double tmp = fabs(f_cell[0][c]) / (factor * molar_rho_ * phi_c[0][c]);
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
  const CompositeVector& pres = *S_->GetFieldData("pressure");
  const CompositeVector& temp = *S_->GetFieldData("temperature");

  // Compute conductivities
  Teuchos::RCP<CompositeVector> kvapor_pres = Teuchos::rcp(new CompositeVector(f->Map()));
  Teuchos::RCP<CompositeVector> kvapor_temp = Teuchos::rcp(new CompositeVector(f->Map()));
  CalculateVaporDiffusionTensor_(kvapor_pres, kvapor_temp);

  // Calculate vapor contribution due to temperature.
  // We assume the same DOFs for pressure and temperature. 
  // We assume that field temperature has already essential BCs.
  op_vapor_->Init();
  op_vapor_diff_->Setup(kvapor_temp, Teuchos::null);
  op_vapor_diff_->UpdateMatrices(Teuchos::null, Teuchos::null);
  op_vapor_diff_->ApplyBCs(false, false);

  // -- Calculate residual due to temperature
  CompositeVector g(*f);
  op_vapor_->ComputeNegativeResidual(temp, g);
  f->Update(1.0, g, 1.0);

  // Calculate vapor contribution due to capillary pressure.
  // We elliminate essential BCs to re-use the local Op for PC.
  op_vapor_->Init();
  op_vapor_diff_->Setup(kvapor_pres, Teuchos::null);
  op_vapor_diff_->UpdateMatrices(Teuchos::null, Teuchos::null);
  op_vapor_diff_->ApplyBCs(false, true);

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
  S_->GetFieldEvaluator("molar_density_gas")->HasFieldChanged(S_.ptr(), passwd_);
  const Epetra_MultiVector& n_g = *S_->GetFieldData("molar_density_gas")->ViewComponent("cell");

  S_->GetFieldEvaluator("porosity")->HasFieldChanged(S_.ptr(), passwd_);
  const Epetra_MultiVector& phi = *S_->GetFieldData("porosity")->ViewComponent("cell");

  S_->GetFieldEvaluator("saturation_liquid")->HasFieldChanged(S_.ptr(), passwd_);
  const Epetra_MultiVector& s_l = *S_->GetFieldData("saturation_liquid")->ViewComponent("cell");

  S_->GetFieldEvaluator("molar_density_liquid")->HasFieldChanged(S_.ptr(), passwd_);
  const Epetra_MultiVector& n_l = *S_->GetFieldData("molar_density_liquid")->ViewComponent("cell");

  S_->GetFieldEvaluator("molar_fraction_gas")->HasFieldChanged(S_.ptr(), passwd_);
  const Epetra_MultiVector& mlf_g = *S_->GetFieldData("molar_fraction_gas")->ViewComponent("cell");

  S_->GetFieldEvaluator("molar_fraction_gas")->HasFieldDerivativeChanged(S_.ptr(), passwd_, "temperature");
  const Epetra_MultiVector& dmlf_g_dt = *S_->GetFieldData("dmolar_fraction_gas_dtemperature")->ViewComponent("cell");

  const Epetra_MultiVector& temp = *S_->GetFieldData("temperature")->ViewComponent("cell");
  const Epetra_MultiVector& pres = *S_->GetFieldData("pressure")->ViewComponent("cell");

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

    kp_cell[0][c] = tmp * mlf_g[0][c] / nRT;
    kt_cell[0][c] = tmp * (dmlf_g_dt[0][c] / atm_pressure_ 
                        +  mlf_g[0][c] * pc / (nRT * temp[0][c]));
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
  Teuchos::RCP<const CompositeVector> mu = S_->GetFieldData("viscosity_liquid");

  // update coefficients
  if (update_upwind == FLOW_UPWIND_UPDATE_ITERATION) {
    op_matrix_diff_->UpdateFlux(*solution, *darcy_flux_copy);
    Epetra_MultiVector& flux = *darcy_flux_copy->ViewComponent("face");
    for (int f = 0; f < nfaces_owned; f++) flux[0][f] /= molar_rho_;
  }
  darcy_flux_copy->ScatterMasterToGhosted("face");

  relperm_->Compute(u->Data(), krel_);
  RelPermUpwindFn func1 = &RelPerm::Compute;
  upwind_->Compute(*darcy_flux_upwind, *u->Data(), bc_model, bc_value, *krel_, *krel_, func1);
  Operators::CellToFace_ScaleInverse(mu, krel_);
  krel_->ScaleMasterAndGhosted(molar_rho_);

  // modify relative permeability coefficient for influx faces
  // UpwindInflowBoundary_New(u->Data());

  relperm_->ComputeDerivative(u->Data(), dKdP_);
  RelPermUpwindFn func2 = &RelPerm::ComputeDerivative;
  upwind_->Compute(*darcy_flux_upwind, *u->Data(), bc_model, bc_value, *dKdP_, *dKdP_, func2);
  Operators::CellToFace_ScaleInverse(mu, dKdP_);
  dKdP_->ScaleMasterAndGhosted(molar_rho_);

  double t_old = tp - dtp;
  UpdateSourceBoundaryData(t_old, tp, *u->Data());

  // create diffusion operators
  op_preconditioner_->Init();
  op_preconditioner_diff_->UpdateMatrices(darcy_flux_copy.ptr(), solution.ptr());
  op_preconditioner_diff_->UpdateMatricesNewtonCorrection(darcy_flux_copy.ptr(), solution.ptr());
  op_preconditioner_diff_->ApplyBCs(true, true);

  // add time derivative
  if (dtp > 0.0) {
    S_->GetFieldEvaluator("water_content")->HasFieldDerivativeChanged(S_.ptr(), passwd_, "pressure");
    CompositeVector& dwc_dp = *S_->GetFieldData("dwater_content_dpressure", "water_content");

    op_acc_->AddAccumulationTerm(*u->Data(), dwc_dp, dtp, "cell");
  }

  // Add vapor diffusion. We assume that the corresponding local operator
  // has been already populated during functional evaluation.
 
  // finalize preconditioner
  op_preconditioner_->AssembleMatrix();
  op_preconditioner_->InitPreconditioner(preconditioner_name_, *preconditioner_list_);
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
  int cell_p(0), cell_r(0);

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
  du.Comm().MaxAll(&buf, &error, 1);  // find the global maximum
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
  double max_sat_pert = 0.25;
  double damping_factor = 0.5;
  double reference_pressure = 101325.0;
  
  if (rp_list_->isSublist("clipping parameters")) {
    Teuchos::ParameterList& clip_list = rp_list_->sublist("clipping parameters");
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
    double tmp = duc[0][c];

    if ((fabs(duc[0][c]) > du_pert_max) && (1 - sat > 1e-5)) {
      if (vo_->getVerbLevel() >= Teuchos::VERB_EXTREME) {
        Teuchos::OSTab tab = vo_->getOSTab();
        *vo_->os() << "clip saturation: c=" << c 
                   << " p=" << uc[0][c]
                   << " dp: " << duc[0][c] << " -> " << du_pert_max << std::endl;
      }

      if (duc[0][c] >= 0.0) duc[0][c] = du_pert_max;
      else duc[0][c] = -du_pert_max;
      
      nsat_clipped++;
    }    
  }

  for (int c = 0; c < ncells_owned; c++) {
    double unew = uc[0][c] - duc[0][c];
    double tmp = duc[0][c];

    if ((unew > atm_pressure_) && (uc[0][c] < atm_pressure_)) {
      if (vo_->getVerbLevel() >= Teuchos::VERB_EXTREME) {
	 *vo_->os() << "pressure change: " << uc[0][c] << " -> " << unew << std::endl;
      }
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



