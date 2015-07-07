/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Linear interpolant of thermal conductivity.
------------------------------------------------------------------------- */

#include "thermal_conductivity_threephase_sutra_hacked.hh"

namespace Amanzi {
namespace Energy {
namespace EnergyRelations {

// registry of method
Utils::RegisteredFactory<ThermalConductivityThreePhase,
                         ThermalConductivityThreePhaseSutraHacked>
        ThermalConductivityThreePhaseSutraHacked::factory_("three-phase sutra hacked");


} // namespace Relations
} // namespace Energy
} // namespace Amanzi
