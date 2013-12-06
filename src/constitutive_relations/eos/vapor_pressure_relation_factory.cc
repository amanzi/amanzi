/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------

   ATS
   Author: Ethan Coon

   Self-registering factory for Vapor Pressure implementations.
   ------------------------------------------------------------------------- */

#include <string>
#include "vapor_pressure_relation_factory.hh"

namespace Amanzi {
namespace Relations {

// explicity instantitate the static data of Factory<VaporPressure>
template<> Utils::Factory<VaporPressureRelation>::map_type* Utils::Factory<VaporPressureRelation>::map_;

// method for instantiating VaporPressure implementations
Teuchos::RCP<VaporPressureRelation> VaporPressureRelationFactory::createVaporPressure(Teuchos::ParameterList& plist) {
  std::string eos_typename = plist.get<std::string>("vapor pressure model type");
  return Teuchos::rcp(CreateInstance(eos_typename, plist));
};

} // namespace
} // namespace

