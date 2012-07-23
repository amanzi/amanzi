/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------

   ATS
   Author: Ethan Coon

   Self-registering factory for Vapor Pressure implementations.
   ------------------------------------------------------------------------- */

#include <string>
#include "vapor_pressure_model_factory.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

// explicity instantitate the static data of Factory<VaporPressureModel>
template<> Utils::Factory<VaporPressureModel>::map_type* Utils::Factory<VaporPressureModel>::map_;

// method for instantiating VaporPressureModel implementations
Teuchos::RCP<VaporPressureModel> VaporPressureModelFactory::createVaporPressureModel(Teuchos::ParameterList& plist) {
  std::string eos_typename = plist.get<std::string>("Vapor Pressure Model Type");
  return Teuchos::rcp(CreateInstance(eos_typename, plist));
};

} // namespace
} // namespace
} // namespace

