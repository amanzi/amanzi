/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Factory of various mesh maps
*/

#ifndef AMANZI_MESH_MAPS_FACTORY_HH_
#define AMANZI_MESH_MAPS_FACTORY_HH_

#include "Mesh.hh"

#include "MeshMaps_FEM.hh"
#include "MeshMaps_PEM.hh"
#include "MeshMaps_VEM.hh"

namespace Amanzi {
namespace WhetStone {

class MeshMapsFactory {
 public:
  explicit MeshMapsFactory(){};
  ~MeshMapsFactory(){};

  // select numerical scheme using its name and order
  std::shared_ptr<MeshMaps>
  Create(const Teuchos::ParameterList& plist,
         const Teuchos::RCP<const AmanziMesh::Mesh>& mesh0,
         const Teuchos::RCP<const AmanziMesh::Mesh>& mesh1)
  {
    std::string name = plist.get<std::string>("map name");
    std::shared_ptr<MeshMaps> maps(NULL);

    if (name == "FEM") {
      maps = std::make_shared<MeshMaps_FEM>(mesh0, mesh1);
    } else if (name == "VEM") {
      maps = std::make_shared<MeshMaps_VEM>(mesh0, mesh1, plist);
    } else if (name == "PEM") {
      maps = std::make_shared<MeshMaps_PEM>(mesh0, mesh1);
    }

    return maps;
  }
};

} // namespace WhetStone
} // namespace Amanzi

#endif
