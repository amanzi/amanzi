/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*
  Flow PK

  Self-registering factory for model implementations.
*/

#ifndef AMANZI_FLOW_MODEL_FACTORY_HH_
#define AMANZI_FLOW_MODEL_FACTORY_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"

namespace Amanzi {
namespace Flow {

template <typename Model>
class ModelFactory : public Utils::Factory<Model> {
 public:
  using Utils::Factory<Model>::CreateInstance;

  Teuchos::RCP<Model> Create(Teuchos::ParameterList& plist)
  {
    std::string name;

    if (plist.isParameter("model")) {
      name = plist.get<std::string>("model");
    }

    return Teuchos::rcp(CreateInstance(name, plist));
  }
};

} // namespace Flow
} // namespace Amanzi

#endif
