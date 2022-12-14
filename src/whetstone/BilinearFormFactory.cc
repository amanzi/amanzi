/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Factory of mimetic methods.
*/

#include <string>

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

#include "MeshLight.hh"

#include "BilinearFormFactory.hh"

namespace Amanzi {
namespace WhetStone {

BilinearFormFactory::map_type* BilinearFormFactory::map_; // initialization

Teuchos::RCP<BilinearForm>
BilinearFormFactory::Create(const Teuchos::ParameterList& plist,
                            const Teuchos::RCP<const AmanziMesh::MeshLight>& mesh)
{
  std::string method = plist.get<std::string>("method");
  BFKey key = method;
  map_type::iterator iter = GetMap()->find(key);

  if (iter == GetMap()->end()) {
    std::cout << "Factory: cannot get item of type: " << method << std::endl;
    for (typename map_type::iterator p = GetMap()->begin(); p != GetMap()->end(); ++p) {
      std::cout << "  option: " << p->first << std::endl;
    }
    return Teuchos::null;
  }
  return Teuchos::rcp(iter->second(plist, mesh));
}

} // namespace WhetStone
} // namespace Amanzi
