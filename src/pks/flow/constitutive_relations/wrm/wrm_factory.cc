/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------

   ATS
   Author: Ethan Coon

   Self-registering factory for WRM implementations.
   ------------------------------------------------------------------------- */

#include <string>
#include "wrm_factory.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

// explicity instantitate the static data of Factory<WRM>
//template<> Factory<WRM>::map_type* Factory<WRM>::map_;
template<> Utils::Factory<WRM>::map_type* Utils::Factory<WRM>::map_;

// method for instantiating WRM implementations
Teuchos::RCP<WRM> WRMFactory::createWRM(Teuchos::ParameterList& plist) {
  std::string wrm_typename = plist.get<std::string>("WRM Type");
  return Teuchos::rcp(CreateInstance(wrm_typename, plist));
};

} // namespace
} // namespace
} // namespace

