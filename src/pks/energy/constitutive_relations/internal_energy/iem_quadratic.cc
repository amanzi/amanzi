/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

u = u0 + a(T - T_ref) + b(T - T_ref)^2 

UNITS: J/{mol,kg}
------------------------------------------------------------------------- */

#include "iem_quadratic.hh"

namespace Amanzi {
namespace Energy {
namespace EnergyRelations {

Utils::RegisteredFactory<InternalEnergyModel,IEMQuadratic> IEMQuadratic::factory_("quadratic");

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
    u0_ = plist_.get<double>("quadratic u_0 [J/kg]");
    ka_ = plist_.get<double>("quadratic a [J/kg-K]");
    kb_ = plist_.get<double>("quadratic b [J/kg-K^2]");
    molar_basis_ = false;

  } else {
    u0_ = plist_.get<double>("quadratic u_0 [J/mol]");
    ka_ = plist_.get<double>("quadratic a [J/mol-K]");
    kb_ = plist_.get<double>("quadratic b [J/mol-K^2]");
    molar_basis_ = true;
  }

  T0_ = plist_.get<double>("Reference temperature [K]", 273.15);
};

} // namespace
} // namespace
} // namespace
