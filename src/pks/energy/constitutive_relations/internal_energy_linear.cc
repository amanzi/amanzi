/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Linear internal energy model -- function of Cv and temperature

See ATS process model documentation's permafrost physical properties
documentation for details.

UNITS: J/{mol/kg}
------------------------------------------------------------------------- */

#include "internal_energy_linear.hh"

namespace Amanzi {
namespace Energy {
namespace EnergyRelations {

InternalEnergyLinear::InternalEnergyLinear(Teuchos::ParameterList& plist) :
    plist_(plist) {
  InitializeFromPlist_();
};

double InternalEnergyLinear::InternalEnergy(double temp) {
  return Cv_ * (temp - T_ref_);
};

void InternalEnergyLinear::InitializeFromPlist_() {
  molar_basis_ = plist_.get<bool>("molar-basis (otherwise, mass-basis)", false);
  Cv_ = plist_.get<double>("heat capacity [J/({kg/mol}-K)]");
  T_ref_ = plist_.get<double>("Reference temperature [K]", 273.15);
};

} // namespace
} // namespace
} // namespace
