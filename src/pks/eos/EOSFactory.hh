/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*
  EOS

  Self-registering factory for EOS implementations.
*/

#ifndef AMANZI_PK_EOS_FACTORY_HH_
#define AMANZI_PK_EOS_FACTORY_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "dbc.hh"
#include "Factory.hh"

namespace Amanzi {
namespace AmanziEOS {

template <typename EOS>
class EOSFactory : public Utils::Factory<EOS> {
 public:
  using Utils::Factory<EOS>::CreateInstance;

  Teuchos::RCP<EOS> Create(Teuchos::ParameterList& plist)
  {
    std::string name;

    if (plist.isParameter("eos type")) {
      name = plist.get<std::string>("eos type");
    } else if (plist.isParameter("com type")) {
      name = plist.get<std::string>("com type");
    }

    auto model = CreateInstance(name, plist);
    if (model == nullptr) AMANZI_ASSERT(false);

    return Teuchos::rcp(model);
  }
};

} // namespace AmanziEOS
} // namespace Amanzi

#endif
