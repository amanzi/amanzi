/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov
*/

/*
  Transport PK

  Self-registering factory for MDM implementations. Default model
  is "scalar".
*/

#include <string>

#include "MDMFactory.hh"

namespace Amanzi {
namespace Transport {

Teuchos::RCP<MDM>
MDMFactory::Create(Teuchos::ParameterList& plist)
{
  std::string mdm_typename = plist.get<std::string>("model", "scalar");
  Teuchos::ParameterList& tmp_list = plist.sublist("parameters for " + mdm_typename);
  return Teuchos::rcp(CreateInstance(mdm_typename, tmp_list));
};

} // namespace Transport
} // namespace Amanzi
