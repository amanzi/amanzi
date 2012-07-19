/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------

   ATS
   Author: Ethan Coon

   Self-registering factory for Viscosity implementations.
   ------------------------------------------------------------------------- */

#include <string>
#include "viscosity_model_factory.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

// explicity instantitate the static data of Factory<ViscosityModel>
template<> Utils::FactoryWithState<ViscosityModel>::map_type* Utils::FactoryWithState<ViscosityModel>::map_;

// method for instantiating ViscosityModel implementations
Teuchos::RCP<ViscosityModel> ViscosityModelFactory::createViscosityModel(Teuchos::ParameterList& plist,
        const Teuchos::Ptr<State>& S) {
  std::string visc_typename = plist.get<std::string>("Viscosity Model Type");
  return Teuchos::rcp(CreateInstance(visc_typename, plist, S));
};

} // namespace
} // namespace
} // namespace

