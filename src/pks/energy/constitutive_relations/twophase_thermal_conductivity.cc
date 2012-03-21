/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Linear interpolant of thermal conductivity.
------------------------------------------------------------------------- */

#include <cmath>
#include "twophase_thermal_conductivity.hh"

namespace Amanzi {
namespace Energy {
namespace EnergyRelations {

TwophaseThermalConductivity::TwophaseThermalConductivity(Teuchos::ParameterList& plist) :
    plist_(plist) {
  InitializeFromPlist_();
};

double TwophaseThermalConductivity::CalculateConductivity(double poro, double sat_liq, double dens_rock) {
  double k_dry = (d_*(1-poro)*k_rock_ + k_gas_*poro)/(d_*(1-poro) + poro);
  double k_sat = pow(k_rock_,(1-poro)) * pow(k_liquid_,poro);
  double kersten = pow(sat_liq + eps_, alpha_);
  return k_dry + (k_sat - k_dry)*kersten;
};

void TwophaseThermalConductivity::InitializeFromPlist_() {
  d_ = 0.053; // unitless empericial parameter

  eps_ = plist_.get<double>("epsilon");
  alpha_ = plist_.get<double>("unsaturated alpha");
  k_rock_ = plist_.get<double>("thermal conductivity of rock");
  k_liquid_ = plist_.get<double>("thermal conductivity of liquid");
  k_gas_ = plist_.get<double>("thermal conductivity of gas");
};

} // namespace Relations
} // namespace Energy
} // namespace Amanzi
