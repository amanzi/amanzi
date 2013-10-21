/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------

   ATS
   Author: Ethan Coon

   Self-registering factory for TC implementations.
   ------------------------------------------------------------------------- */

#include "thermal_conductivity_threephase_factory.hh"

namespace Amanzi {
namespace Energy {
namespace EnergyRelations {

// explicity instantitate the static data of Factory<EOS>
template<> Utils::Factory<ThermalConductivityThreePhase>::map_type* Utils::Factory<ThermalConductivityThreePhase>::map_;

} // namespace
} // namespace
} // namespace

