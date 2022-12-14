/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*
  Flow PK

  Self-registering factory for WRM implementations.
*/

#include <string>
#include "WRMFactory.hh"

namespace Amanzi {
namespace Flow {

// method for instantiating WRM implementations
Teuchos::RCP<WRM>
WRMFactory::Create(Teuchos::ParameterList& plist)
{
  std::string wrm_typename = plist.get<std::string>("water retention model");
  return Teuchos::rcp(CreateInstance(wrm_typename, plist));
};

} // namespace Flow
} // namespace Amanzi
