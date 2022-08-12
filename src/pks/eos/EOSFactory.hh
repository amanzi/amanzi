/*
  EOS
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

  Self-registering factory for EOS implementations.
*/

#ifndef AMANZI_PK_EOS_FACTORY_HH_
#define AMANZI_PK_EOS_FACTORY_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"

namespace Amanzi {
namespace AmanziEOS {

template<typename EOS>
class EOSFactory : public Utils::Factory<EOS> {
 public:
  using Utils::Factory<EOS>::CreateInstance;

  Teuchos::RCP<EOS> Create(Teuchos::ParameterList& plist) {
    std::string name;

    if (plist.isParameter("eos type")) {
      name = plist.get<std::string>("eos type");
    }
    else if (plist.isParameter("com type")) {
      name = plist.get<std::string>("com type");
    }

    return Teuchos::rcp(CreateInstance(name, plist));
  }
};

}  // namespace AmanziEOS
}  // namespace Amanzi

#endif
