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

#ifndef AMANZI_BILINEAR_FORM_FACTORY_HH_
#define AMANZI_BILINEAR_FORM_FACTORY_HH_

#include <string>
#include <utility>

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

#include "MeshLight.hh"

#include "BilinearForm.hh"

namespace Amanzi {
namespace WhetStone {

typedef std::string BFKey;

class BilinearFormFactory {
  typedef std::map<BFKey,
                   BilinearForm* (*)(const Teuchos::ParameterList&,
                                     const Teuchos::RCP<const AmanziMesh::MeshLight>&)>
    map_type;

 public:
  static Teuchos::RCP<BilinearForm> Create(const Teuchos::ParameterList& plist,
                                           const Teuchos::RCP<const AmanziMesh::MeshLight>& mesh);

 protected:
  static map_type* GetMap()
  {
    if (!map_) map_ = new map_type;
    return map_;
  }

 private:
  static map_type* map_;
};


template <typename TDerived>
BilinearForm*
CreateT(const Teuchos::ParameterList& plist, const Teuchos::RCP<const AmanziMesh::MeshLight>& mesh)
{
  return new TDerived(plist, mesh);
}

template <typename TDerived>
class RegisteredFactory : public BilinearFormFactory {
 public:
  RegisteredFactory(const std::string& s)
  {
    for (auto it = BilinearFormFactory::GetMap()->begin();
         it != BilinearFormFactory::GetMap()->end();
         ++it) {}
    BilinearFormFactory::GetMap()->insert(
      std::pair<BFKey,
                BilinearForm* (*)(const Teuchos::ParameterList&,
                                  const Teuchos::RCP<const AmanziMesh::MeshLight>&)>(
        s, &CreateT<TDerived>));
    for (auto it = BilinearFormFactory::GetMap()->begin();
         it != BilinearFormFactory::GetMap()->end();
         ++it) {}
  }
};

} // namespace WhetStone
} // namespace Amanzi

#endif
