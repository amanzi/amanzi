/*
  EOS
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

  Self-registering factory for EOS implementations.
*/

#include <string>
#include "EOSFactory.hh"

// explicity instantitate the static data of Factory<EOS>
// template<> 
// Amanzi::Utils::Factory<Amanzi::EOS>::map_type* 
//    Amanzi::Utils::Factory<Amanzi::EOS>::map_;

namespace Amanzi {
namespace AmanziEOS {

// method for instantiating EOS implementations
Teuchos::RCP<EOS> EOSFactory::CreateEOS(Teuchos::ParameterList& plist) {
  std::string eos_typename = plist.get<std::string>("eos type");
  return Teuchos::rcp(CreateInstance(eos_typename, plist));
};

}  // namespace AmanziEOS
}  // namespace Amanzi

