/* -*-  mode: c++; indent-tabs-mode: nil -*- */
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

// method for instantiating WRM implementations
Teuchos::RCP<WRM> WRMFactory::createWRM(Teuchos::ParameterList& plist) {
  std::string wrm_typename = plist.get<std::string>("WRM Type");
  return Teuchos::rcp(CreateInstance(wrm_typename, plist));
};

} // namespace
} // namespace
} // namespace

