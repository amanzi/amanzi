/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Linear interpolant of thermal conductivity.
------------------------------------------------------------------------- */

#include <cmath>
#include "thermal_conductivity_threephase_wetdry.hh"

namespace Amanzi {
namespace Energy {
namespace EnergyRelations {

// registry of method
Utils::RegisteredFactory<ThermalConductivityThreePhase,
                         ThermalConductivityThreePhaseWetDry>
        ThermalConductivityThreePhaseWetDry::factory_("three-phase wet/dry");

ThermalConductivityThreePhaseWetDry::ThermalConductivityThreePhaseWetDry(
      Teuchos::ParameterList& plist) : plist_(plist) {
  InitializeFromPlist_();
};

double ThermalConductivityThreePhaseWetDry::ThermalConductivity(double poro,
        double sat_liq, double sat_ice, double temp) {
  double Ki = 831.51/(pow(temp, 1.0552));
  double Kl = 0.5611;
  double k_sat_f = k_sat_u_ * pow(Ki/Kl, poro);

  double kersten_u = pow(sat_liq + eps_, alpha_u_);
  double kersten_f = pow(sat_ice + eps_, alpha_f_);
  return kersten_f * k_sat_f + kersten_u * k_sat_u_
    + (1.0 - kersten_f - kersten_u) * k_dry_;
};

void ThermalConductivityThreePhaseWetDry::InitializeFromPlist_() {

  eps_ = plist_.get<double>("epsilon", 1.e-10);
  alpha_u_ = plist_.get<double>("unsaturated alpha unfrozen");
  alpha_f_ = plist_.get<double>("unsaturated alpha frozen");
  k_dry_ = plist_.get<double>("thermal conductivity, dry");
  k_sat_u_ = plist_.get<double>("thermal conductivity, saturated (unfrozen)");
};

} // namespace Relations
} // namespace Energy
} // namespace Amanzi
