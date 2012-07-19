/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon


Simple model of two-phase thermal conductivity, based upon:

- Interpolation between saturated and dry conductivities via a Kersten number.
- Power-law Kersten number.

See ATS process model documentation's permafrost model for details.
------------------------------------------------------------------------- */

#include <cmath>
#include "thermal_conductivity_twophase_wetdry.hh"

namespace Amanzi {
namespace Energy {
namespace EnergyRelations {

// registry of method
Utils::RegisteredFactory<ThermalConductivityTwoPhase,ThermalConductivityTwoPhaseWetDry>
        ThermalConductivityTwoPhaseWetDry::factory_("two-phase wet/dry");

// constructor
ThermalConductivityTwoPhaseWetDry::ThermalConductivityTwoPhaseWetDry(
      Teuchos::ParameterList& plist) : plist_(plist) {
  InitializeFromPlist_();
};

// do the physics
double ThermalConductivityTwoPhaseWetDry::CalculateConductivity(double poro,
        double sat_liq) {
  double kersten = pow(sat_liq + eps_, alpha_);
  return k_dry_ + (k_wet_ - k_dry_)*kersten;
};

// initialization
void ThermalConductivityTwoPhaseWetDry::InitializeFromPlist_() {
  eps_ = plist_.get<double>("epsilon", 1.e-10);
  alpha_ = plist_.get<double>("unsaturated alpha", 1.0);
  k_dry_ = plist_.get<double>("thermal conductivity, dry");
  k_wet_ = plist_.get<double>("thermal conductivity, wet");
};

} // namespace Relations
} // namespace Energy
} // namespace Amanzi
