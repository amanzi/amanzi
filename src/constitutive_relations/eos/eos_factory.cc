/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------

   ATS
   Author: Ethan Coon

   Self-registering factory for EOS implementations.
   ------------------------------------------------------------------------- */

#include <string>
#include "eos_factory.hh"

namespace Amanzi {
namespace Relations {

// explicity instantitate the static data of Factory<EOS>
//template<> Factory<EOS>::map_type* Factory<EOS>::map_;
template<> Utils::Factory<EOS>::map_type* Utils::Factory<EOS>::map_;

// method for instantiating EOS implementations
Teuchos::RCP<EOS> EOSFactory::createEOS(Teuchos::ParameterList& plist) {
  std::string eos_typename = plist.get<std::string>("EOS type");
  return Teuchos::rcp(CreateInstance(eos_typename, plist));
};

} // namespace
} // namespace

