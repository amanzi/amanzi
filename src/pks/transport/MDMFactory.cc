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
  std::string mdm_typename;
  if (plist.isParameter("model")) {
    mdm_typename = plist.get<std::string>("model");
  } else {
    mdm_typename = plist.get<std::string>("mechanical dispersion type", "isotropic");
  }

  Teuchos::ParameterList& tmp_list = plist.sublist(mdm_typename + " parameters");
  return Teuchos::rcp(CreateInstance(mdm_typename, tmp_list));
};

} // namespace Transport
} // namespace Amanzi
