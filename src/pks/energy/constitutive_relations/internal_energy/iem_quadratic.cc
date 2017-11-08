/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

u = u0 + a(T - T_ref) + b(T - T_ref)^2 

UNITS: MJ/{mol,kg}
------------------------------------------------------------------------- */

#include "iem_quadratic.hh"

namespace Amanzi {
namespace Energy {

IEMQuadratic::IEMQuadratic(Teuchos::ParameterList& plist) :
    plist_(plist) {
  InitializeFromPlist_();
};

double IEMQuadratic::InternalEnergy(double temp) {
  double dT = temp - T0_;
  return u0_ + (ka_ + kb_*dT) * dT;
};

double IEMQuadratic::DInternalEnergyDT(double temp) {
  double dT = temp - T0_;
  return ka_ + 2.0*kb_*dT;
};

void IEMQuadratic::InitializeFromPlist_() {
  if (plist_.isParameter("quadratic u_0 [J/kg]")) {
    u0_ = 1.e-6 * plist_.get<double>("quadratic u_0 [J/kg]");
    ka_ = 1.e-6 * plist_.get<double>("quadratic a [J/kg-K]");
    kb_ = 1.e-6 * plist_.get<double>("quadratic b [J/kg-K^2]");
    molar_basis_ = false;

  } else {
    u0_ = 1.e-6 * plist_.get<double>("quadratic u_0 [J/mol]");
    ka_ = 1.e-6 * plist_.get<double>("quadratic a [J/mol-K]");
    kb_ = 1.e-6 * plist_.get<double>("quadratic b [J/mol-K^2]");
    molar_basis_ = true;
  }

  T0_ = plist_.get<double>("reference temperature [K]", 273.15);
};

} // namespace
} // namespace
