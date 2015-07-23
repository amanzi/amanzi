/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------

   ATS
   Author: Ethan Coon

   Self-registering factory for EOS implementations.
   ------------------------------------------------------------------------- */

#include <string>
#include "eos_factory.hh"

// explicity instantitate the static data of Factory<EOS>
template<> 
Amanzi::Utils::Factory<Amanzi::Relations::EOS>::map_type* 
Amanzi::Utils::Factory<Amanzi::Relations::EOS>::map_;

namespace Amanzi {
namespace Relations {

// method for instantiating EOS implementations
Teuchos::RCP<EOS> EOSFactory::createEOS(Teuchos::ParameterList& plist) {
  std::string eos_typename = plist.get<std::string>("EOS type");
  return Teuchos::rcp(CreateInstance(eos_typename, plist));
};

} // namespace
} // namespace

