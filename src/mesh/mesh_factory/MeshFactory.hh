/*
  Copyright 2010-201x held jointly by LANL, ORNL, LBNL, and PNNL.
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
  // This is not virtual -- it hides the framework factory's creates on purpose.
  // It calls the hidden create, which creates a MeshFramework, then creates a
  // MeshCache, returning a different pointer type.  It cannot overload because
  // the return type is different!
  //
  Teuchos::RCP<Mesh> create(const std::string& filename)
  {
    Teuchos::RCP<MeshFramework> mesh_fw = MeshFrameworkFactory::create(filename);
    auto mesh = Teuchos::rcp(
      new Mesh(mesh_fw, Teuchos::rcp(new AmanziMesh::MeshFrameworkAlgorithms()), Teuchos::null));
    return mesh;
  }

  Teuchos::RCP<Mesh> create(const double x0,
                            const double y0,
                            const double z0,
                            const double x1,
                            const double y1,
                            const double z1,
                            const int nx,
                            const int ny,
                            const int nz)
  {
    Teuchos::RCP<MeshFramework> mesh_fw =
      MeshFrameworkFactory::create(x0, y0, z0, x1, y1, z1, nx, ny, nz);
    auto mesh = Teuchos::rcp(
      new Mesh(mesh_fw, Teuchos::rcp(new AmanziMesh::MeshFrameworkAlgorithms()), Teuchos::null));
    return mesh;
  }

  Teuchos::RCP<Mesh> create(const double x0,
                            const double y0,
                            const double x1,
                            const double y1,
                            const int nx,
                            const int ny)
  {
    Teuchos::RCP<MeshFramework> mesh_fw = MeshFrameworkFactory::create(x0, y0, x1, y1, nx, ny);
    auto mesh = Teuchos::rcp(
      new Mesh(mesh_fw, Teuchos::rcp(new AmanziMesh::MeshFrameworkAlgorithms()), Teuchos::null));
    return mesh;
  }

  Teuchos::RCP<Mesh> create(const Teuchos::ParameterList& gen_plist)
  {
    Teuchos::RCP<MeshFramework> mesh_fw = MeshFrameworkFactory::create(gen_plist);
    auto mesh = Teuchos::rcp(
      new Mesh(mesh_fw, Teuchos::rcp(new AmanziMesh::MeshFrameworkAlgorithms()), Teuchos::null));
    return mesh;
  }

  // Need a special one for extracted meshes, since they differ not just in
  // return type but also in parent mesh type (MeshCache vs MeshFramework.
  Teuchos::RCP<Mesh> create(const Teuchos::RCP<const Mesh>& parent_mesh,
                            const Entity_ID_View& setids,
                            const Entity_kind setkind,
                            const bool flatten = false);

  // Another for extracted meshes, this time with set names instead of IDs
  Teuchos::RCP<Mesh> create(const Teuchos::RCP<const Mesh>& parent_mesh,
                            const std::vector<std::string>& setnames,
                            const Entity_kind setkind,
                            const bool flatten = false);


  // Create a logical mesh, which is reduced functionality and reduced
  // geometric/topologic connections.
  Teuchos::RCP<Mesh> createLogical(Teuchos::ParameterList& log_plist);


  // Create a 1D Column Mesh from a columnar structured volume mesh.
  //
  Teuchos::RCP<Mesh> createColumn(const Teuchos::RCP<Mesh>& parent, int col_id);

  // Create a MeshSurfaceCell from a MeshFrameworkColumn
  Teuchos::RCP<Mesh> createSurfaceCell(const Teuchos::RCP<const Mesh>& parent);
};

} // namespace AmanziMesh
} // namespace Amanzi
