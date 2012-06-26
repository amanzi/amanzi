/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Linear internal energy model -- function of Cv and temperature

u = C * (T - T_ref)

UNITS: J/{mol/kg}
------------------------------------------------------------------------- */

#include "iem_linear.hh"

namespace Amanzi {
namespace Energy {
namespace EnergyRelations {

Utils::RegisteredFactory<InternalEnergyModel,IEMLinear> IEMLinear::factory_("linear");

IEMLinear::IEMLinear(Teuchos::ParameterList& plist) :
    plist_(plist) {
  InitializeFromPlist_();
};

double IEMLinear::InternalEnergy(double temp) {
  return Cv_ * (temp - T_ref_);
};

void IEMLinear::InitializeFromPlist_() {
  if (plist_.isParameter("Heat capacity [J/kg-K]")) {
    Cv_ = plist_.get<double>("Heat capacity [J/kg-K]");
    molar_basis_ = false;
  } else {
    Cv_ = plist_.get<double>("Heat capacity [J/mol-K]");
    molar_basis_ = true;
  }

  T_ref_ = plist_.get<double>("Reference temperature [K]", 273.15);
};

} // namespace
} // namespace
} // namespace
