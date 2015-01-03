/*
  This is the EOS component of the Amanzi code.
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

  Self-registering factory for EOS implementations.
*/

#include <string>
#include "eos_factory.hh"

// explicity instantitate the static data of Factory<EOS>
// template<> 
// Amanzi::Utils::Factory<Amanzi::Relations::EOS>::map_type* 
//    Amanzi::Utils::Factory<Amanzi::Relations::EOS>::map_;


namespace Amanzi {
namespace Relations {

// method for instantiating EOS implementations
Teuchos::RCP<EOS> EOSFactory::createEOS(Teuchos::ParameterList& plist) {
  std::string eos_typename = plist.get<std::string>("EOS type");
  return Teuchos::rcp(CreateInstance(eos_typename, plist));
};

}  // namespace Relations
}  // namespace Amanzi

