/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  Ugly hackjob to enable direct evaluation of the full model, on a single
  WRM/region.  This is bypassing much of the "niceness" of the framework, but
  seems necessary for solving a cell-wise correction equation.

  Uses intensive, not extensive, forms.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "eos.hh"
#include "wrm_permafrost_model.hh"
#include "vapor_pressure_relation.hh"
#include "pc_ice_water.hh"
#include "pc_liq_atm.hh"
#include "iem.hh"
#include "iem_water_vapor.hh"

#include "permafrost_model.hh"

namespace Amanzi {

bool PermafrostModel::IsSetUp() {
  if (wrm_ == Teuchos::null) return false;
  if (liquid_eos_ == Teuchos::null) return false;
  if (gas_eos_ == Teuchos::null) return false;
  if (ice_eos_ == Teuchos::null) return false;
  if (pc_i_ == Teuchos::null) return false;
  if (pc_l_ == Teuchos::null) return false;
  if (vpr_ == Teuchos::null) return false;
  if (liquid_iem_ == Teuchos::null) return false;
  if (gas_iem_ == Teuchos::null) return false;
  if (ice_iem_ == Teuchos::null) return false;
  if (rock_iem_ == Teuchos::null) return false;
  if (rho_rock_ < 0.) return false;
  if (p_atm_ < -1.e-10) return false;
  return true;
}


int PermafrostModel::EvaluateEnergyAndWaterContent(double T, double p, double poro, AmanziGeometry::Point& result) {
  double eff_p = std::max(p_atm_, p);

  double rho_l = liquid_eos_->MolarDensity(T,eff_p);
  double mass_rho_l = liquid_eos_->MassDensity(T,eff_p);
  double rho_i = ice_eos_->MolarDensity(T,eff_p);
  double rho_g = gas_eos_->MolarDensity(T,eff_p);

  double omega = vpr_->SaturatedVaporPressure(T)/p_atm_;

  double pc_i;
  if (pc_i_->IsMolarBasis()) {
    pc_i = pc_i_->CapillaryPressure(T, rho_l);
  } else {
    double mass_rho_l = liquid_eos_->MassDensity(T,eff_p);
    pc_i = pc_i_->CapillaryPressure(T, mass_rho_l);
  }

  double pc_l = pc_l_->CapillaryPressure(p, p_atm_);

  double sats[3];
  wrm_->saturations(pc_l, pc_i, sats);
  double s_g = sats[0];
  double s_l = sats[1];
  double s_i = sats[2];

  double u_l = liquid_iem_->InternalEnergy(T);
  double u_g = gas_iem_->InternalEnergy(T, omega);
  double u_i = ice_iem_->InternalEnergy(T);

  double u_rock = rock_iem_->InternalEnergy(T);

  // water content
  result[1] = poro * (rho_l * s_l + rho_i * s_i + rho_g * s_g * omega);

  // energy
  result[0] = poro * (u_l * rho_l * s_l + u_i * rho_i * s_i + u_g * rho_g * s_g)
      + (1.0 - poro) * (rho_rock_ * u_rock);

  return 0;
}


int PermafrostModel::EvaluateEnergyAndWaterContentAndJacobian(double T, double p,
        double poro, AmanziGeometry::Point& result, WhetStone::Tensor& jac) {
  return EvaluateEnergyAndWaterContentAndJacobian_FD_(T, p, poro, result, jac);
}


int PermafrostModel::EvaluateEnergyAndWaterContentAndJacobian_FD_(double T, double p,
        double poro, AmanziGeometry::Point& result, WhetStone::Tensor& jac) {
  double eps_T = 1.e-7;
  double eps_p = 1.e-3;

  EvaluateEnergyAndWaterContent(T, p, poro, result);

  AmanziGeometry::Point test(result);
  // d / dT
  EvaluateEnergyAndWaterContent(T+eps_T, p, poro, test);
  jac(0,0) = (test[0] - result[0]) / (eps_T);
  jac(1,0) = (test[1] - result[1]) / (eps_T);

  // d / dp
  EvaluateEnergyAndWaterContent(T, p + eps_p, poro, test);
  jac(0,1) = (test[0] - result[0]) / (eps_p);
  jac(1,1) = (test[1] - result[1]) / (eps_p);

  /*
  // The Jacobian typically goes singular as de/dp and dwc/dp get small,
  // which occurs at ~ 5 > T > 0 and 101 kPa < p < p_atm.
  // If singular, try increasing size and going to a centered diff.
  int count = 0;
  while ((std::abs(jac(0,1)) + std::abs(jac(1,1)) < 1.e-12)
         && count < 4) {
    std::cout << " Zero determinate in calc Jacobian:" << std::endl;
    std::cout << "   [" << jac(0,0) << "," << jac(0,1) << "]" << std::endl;
    std::cout << "   [" << jac(1,0) << "," << jac(1,1) << "]" << std::endl;
    std::cout << "  at T,p = " << T << ", " << p << std::endl;
    std::cout << "  from eps_p = " << eps_p << std::endl;
    eps_p *= 10;
    AmanziGeometry::Point test2(result);
    EvaluateEnergyAndWaterContent(T, p + eps_p, poro, test);
    EvaluateEnergyAndWaterContent(T, p - eps_p, poro, test2);
    jac(0,1) = (test[0] - test2[0]) / (2*eps_p);
    jac(1,1) = (test[1] - test2[1]) / (2*eps_p);
    ++count;
  }

  if (count > 0) {
    std::cout << " Stopping in calc Jacobian:" << std::endl;
    std::cout << "   [" << jac(0,0) << "," << jac(0,1) << "]" << std::endl;
    std::cout << "   [" << jac(1,0) << "," << jac(1,1) << "]" << std::endl;
    std::cout << "  at T,p = " << T << ", " << p << std::endl;
    std::cout << "  from eps_p = " << eps_p << std::endl;
  }
  */

  return 0;
}


/* ----------------------------------------------------------------------
Solves a given energy and water content (at a given, fixed porosity), for
temperature and pressure.

Note this cannot really work for a saturated cell, as d_wc / d{T,p} = 0

Error codes:

  1 = Singular Jacobian at some point in the evaluation.  Often this is
      because the cell is saturated or becomes saturated and below
      freezing.
  2 = Iteration did not converge in max_steps (hard-coded to be 100 for
      now).
---------------------------------------------------------------------- */
int PermafrostModel::InverseEvaluate(double energy, double wc, double poro,
        double& T, double& p) {

  // -- scaling for the norms
  double wc_scale = 10.;
  double e_scale = 10000.;
  double tol = 1.e-6;
  double max_steps = 100;
  double stepnum = 0;

  // get the initial residual
  AmanziGeometry::Point res(2);
  WhetStone::Tensor jac(2,2);
  int ierr = EvaluateEnergyAndWaterContentAndJacobian(T,p,poro,res,jac);
  if (ierr) {
    std::cout << "Error in evaluation: " << ierr << std::endl;
    return ierr + 10;
  }

  AmanziGeometry::Point f(2);
  f[0] = energy;
  f[1] = wc;
  res = res - f;

  // check convergence
  AmanziGeometry::Point scaled_res(res);
  scaled_res[0] = res[0] / e_scale; scaled_res[1] = res[1] / wc_scale;
  double norm = AmanziGeometry::norm(scaled_res);

  bool converged = norm < tol;

  // workspace
  AmanziGeometry::Point x(2);
  AmanziGeometry::Point x_tmp(2);
  x[0] = T; x[1] = p;
  x_tmp[0] = T; x_tmp[1] = p;

  while (!converged) {
    // calculate the update size
    double detJ = jac.determinant();
    AmanziGeometry::Point correction;

    if (std::abs(detJ) < 1.e-12) {
      std::cout << " Zero determinate of Jacobian:" << std::endl;
      std::cout << "   [" << jac(0,0) << "," << jac(0,1) << "]" << std::endl;
      std::cout << "   [" << jac(1,0) << "," << jac(1,1) << "]" << std::endl;
      std::cout << "  at T,p = " << x_tmp[0] << ", " << x_tmp[1] << std::endl;
      std::cout << "  with res(e,wc) = " << res[0] << ", " << res[1] << std::endl;
      return 1;

    } else {
      jac.inverse();
      correction = jac * res;
    }


    // perform the update
    x_tmp = x - correction;
    ierr = EvaluateEnergyAndWaterContentAndJacobian(x_tmp[0],x_tmp[1],poro,res,jac);
    if (ierr) {
      std::cout << "Error in evaluation: " << ierr << std::endl;
      return ierr + 10;
    }

    res = res - f;

    // check convergence and damping
    scaled_res[0] = res[0] / e_scale; scaled_res[1] = res[1] / wc_scale;
    double norm_new = AmanziGeometry::norm(scaled_res);

    double damp = 1.;
    bool backtracking_required = false;
    while (norm_new > norm) {
      backtracking_required = true;

      // backtrack
      damp *= 0.5;
      x_tmp = x - (damp * correction);

      // evaluate the damped value
      ierr = EvaluateEnergyAndWaterContent(x_tmp[0],x_tmp[1],poro,res);
      if (ierr) {
        std::cout << "Error in evaluation: " << ierr << std::endl;
        return ierr + 10;
      }
      res = res - f;

      // check the new residual
      scaled_res[0] = res[0] / e_scale; scaled_res[1] = res[1] / wc_scale;
      norm_new = AmanziGeometry::norm(scaled_res);
    }

    if (backtracking_required) {
      // must recalculate the Jacobian at the new value
      ierr = EvaluateEnergyAndWaterContentAndJacobian(x_tmp[0],x_tmp[1],poro,res,jac);
      if (ierr) {
        std::cout << "Error in evaluation: " << ierr << std::endl;
        return ierr + 10;
      }
      res = res - f;
    }

    // iterate
    x = x_tmp;
    norm = norm_new;

    converged = norm < tol;
    stepnum++;
    if (stepnum > max_steps && !converged) {
      std::cout << " Nonconverged after " << max_steps << " steps with norm (tol) "
                << norm << " (" << tol << ")" << std::endl;
      return 2;
    }
  }

  T = x[0];
  p = x[1];
  return 0;
}

}
