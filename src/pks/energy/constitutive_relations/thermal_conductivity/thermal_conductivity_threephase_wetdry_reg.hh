/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Linear interpolant of thermal conductivity.
------------------------------------------------------------------------- */

#include "thermal_conductivity_threephase_wetdry.hh"

namespace Amanzi {
namespace Energy {
namespace EnergyRelations {

// registry of method
Utils::RegisteredFactory<ThermalConductivityThreePhase,
                         ThermalConductivityThreePhaseWetDry>
        ThermalConductivityThreePhaseWetDry::factory_("three-phase wet/dry");

} // namespace Relations
} // namespace Energy
} // namespace Amanzi
