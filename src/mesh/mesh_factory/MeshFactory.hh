/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/
#pragma once

#include "MeshFrameworkFactory.hh"
#include "MeshCache.hh"

namespace Amanzi {
namespace AmanziMesh {

struct MeshFactory : public MeshFrameworkFactory {
  using MeshFrameworkFactory::MeshFrameworkFactory;

  // The same create methods as supported in MeshFrameworkFactory are supported
  // here -- see MeshFrameworkFactory.hh for documentation.
  //
  // This is not virtual -- it hides the framework factory's create on purpose.
  // It calls the hidden create, which creates a MeshFramework, then creates a
  // MeshCache, returning a different pointer type.  It cannot overload because
  // the return type is different!
  //
  template<typename ...Args>
  Teuchos::RCP<Mesh>
  create(Args... args) {
    Teuchos::RCP<MeshFramework> mesh_fw = MeshFrameworkFactory::create(args...);
    auto mesh = Teuchos::rcp(new Mesh(mesh_fw, Teuchos::rcp(new AmanziMesh::MeshFrameworkAlgorithms())));
    return mesh;
  }

  // Need a special one for extracted meshes, since they differ not just in
  // return type but also in parent mesh type (MeshCache vs MeshFramework.
  Teuchos::RCP<Mesh>
  create(const Teuchos::RCP<const Mesh>& parent_mesh,
         const Entity_ID_View& setids,
         const Entity_kind setkind,
         const bool flatten=false);

  // Another for extracted meshes, this time with set names instead of IDs
  Teuchos::RCP<Mesh>
  create(const Teuchos::RCP<const Mesh>& parent_mesh,
         const std::vector<std::string>& setnames,
         const Entity_kind setkind,
         const bool flatten=false);

  // Create a 1D Column Mesh from a columnar structured volume mesh.
  //
  Teuchos::RCP<Mesh>
  createColumn(const Teuchos::RCP<Mesh>& parent, int col_id);

  // Create a MeshSurfaceCell from a MeshFrameworkColumn
  Teuchos::RCP<Mesh>
  createSurfaceCell(const Teuchos::RCP<const Mesh>& parent);

};

} // namespace AmanziMesh
} // namespace Amanzi
