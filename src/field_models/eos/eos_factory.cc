/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------

   ATS
   Author: Ethan Coon

   Self-registering factory for EOS implementations.
   ------------------------------------------------------------------------- */

#include <string>
#include "eos_factory.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

// explicity instantitate the static data of Factory<EOS>
template<> Utils::FactoryWithState<EOS>::map_type* Utils::FactoryWithState<EOS>::map_;

// method for instantiating EOS implementations
Teuchos::RCP<EOS> EOSFactory::createEOS(Teuchos::ParameterList& plist,
        const Teuchos::Ptr<State>& S) {
  std::string eos_typename = plist.get<std::string>("EOS Type");
  return Teuchos::rcp(CreateInstance(eos_typename, plist, S));
};

} // namespace
} // namespace
} // namespace

