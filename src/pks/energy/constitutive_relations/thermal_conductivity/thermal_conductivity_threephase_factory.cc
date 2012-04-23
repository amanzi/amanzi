/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------

   ATS
   Author: Ethan Coon

   Self-registering factory for TC implementations.
   ------------------------------------------------------------------------- */

#include <string>
#include "thermal_conductivity_threephase_factory.hh"

namespace Amanzi {
namespace Energy {
namespace EnergyRelations {

// explicity instantitate the static data of Factory<EOS>
template<> Utils::Factory<ThermalConductivityThreePhase>::map_type* Utils::Factory<ThermalConductivityThreePhase>::map_;

// method for instantiating implementations
Teuchos::RCP<ThermalConductivityThreePhase> ThermalConductivityThreePhaseFactory::createThermalConductivityModel(Teuchos::ParameterList& plist) {
  std::string tc_typename = plist.get<std::string>("Thermal Conductivity Type");
  return Teuchos::rcp(CreateInstance(tc_typename, plist));
};

} // namespace
} // namespace
} // namespace

