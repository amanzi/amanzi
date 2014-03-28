/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------

   ATS
   Author: Ethan Coon

   Self-registering factory for EOS implementations.
   ------------------------------------------------------------------------- */

#include <string>
#include "surface_relperm_model_factory.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

// explicity instantitate the static data of Factory<SurfaceRelPermModel>
template<> Utils::Factory<SurfaceRelPermModel>::map_type* Utils::Factory<SurfaceRelPermModel>::map_;

// method for instantiating SurfaceRelPermModel implementations
Teuchos::RCP<SurfaceRelPermModel> SurfaceRelPermModelFactory::createModel(Teuchos::ParameterList& plist) {
  std::string type_name = plist.get<std::string>("surface rel perm model type");
  return Teuchos::rcp(CreateInstance(type_name, plist));
};

} // namespace
} // namespace
} // namespace

