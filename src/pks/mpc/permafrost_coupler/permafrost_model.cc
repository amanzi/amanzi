/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  Ugly hackjob to enable direct evaluation of the full model, on a single
  WRM/region.  This is bypassing much of the "niceness" of the framework, but
  seems necessary for solving a cell-wise correction equation.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "eos.hh"
#include "wrm.hh"
#include "vapor_pressure_relation.hh"
#include "pc_ice_water.hh"
#include "iem.hh"

#include "permafrost_model.hh"

namespace Amanzi {

bool PermafrostModel::IsSetUp() {
  if (wrm_ == Teuchos::null) return false;
  if (liquid_eos_ == Teuchos::null) return false;
  if (gas_eos_ == Teuchos::null) return false;
  if (ice_eos_ == Teuchos::null) return false;
  if (pc_i_ == Teuchos::null) return false;
  if (vpr_ == Teuchos::null) return false;
  if (liquid_iem_ == Teuchos::null) return false;
  if (gas_iem_ == Teuchos::null) return false;
  if (ice_iem_ == Teuchos::null) return false;
  if (rock_iem_ == Teuchos::null) return false;
  if (rho_rock_ < 0.) return false;
  if (p_atm_ < -1.e-10) return false;
  return true;
}


int EvaluateEnergyAndWaterContent(double T, double p, double poro, AmanziMesh::Point& result) {
  double rho_l = liquid_eos_->MolarDensity(T,p);
  double rho_i = ice_eos_->MolarDensity(T,p);
  double rho_g = gas_eos_->MolarDensity(T,p);

  double omega = vpr_->SaturatedVaporPressure(T);

  double pc_i = pc_i_->CapillaryPressure(T, rho_l);
  double one_on_A = wrm_->saturation(pc_i);
  double one_on_B = wrm_->saturation(p_atm_ - p);

  double s_l = 1.0 / (1.0/one_on_A + 1.0/one_on_B - 1.0);
  double s_i = s_l * (1.0/one_on_A - 1.0);
  double s_g = s_l * (1.0/one_on_B - 1.0);

  double u_l = liquid_iem_->InternalEnergy(T);
  double u_g = gas_iem_->InternalEnergy(T);
  double u_i = ice_iem_->InternalEnergy(T);

  double u_rock = rock_iem_->InternalEnergy(T);

  // water content
  result[1] = poro * (rho_l * s_l + rho_i * s_i + rho_g * s_g * omega);

  // energy
  result[0] = poro * (u_l * rho_l * s_l + u_i * rho_i * s_i + u_g * rho_g * s_g * omega)
      + (1.0 - poro) * (rho_rock_ * u_rock);

  return 0;
}



int EvaluateEnergyAndWaterContentAndJacobian(double T, double p, double poro,
        AmanziGeometry::Point& result, WhetStone::Tensor& jac) {
  double eps_T = 1.e-4;
  double eps_p = 1.e-1;

  EvaluateWaterContentAndEnergy(T, p, poro, result);

  AmanziGeometry::Point test(result);
  EvaluateWaterContentAndEnergy(T+eps_T, p, poro, test);

  jac(0,0) = (test[0] - result[0]) / eps_T;
  jac(1,0) = (test[1] - result[1]) / eps_T;

  EvaluateWaterContentAndEnergy(T, p + eps_p, poro, test);

  jac(0,1) = (test[0] - result[0]) / eps_p;
  jac(1,1) = (test[1] - result[1]) / eps_p;

  return 0;
}


}
