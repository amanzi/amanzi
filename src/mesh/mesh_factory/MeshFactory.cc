/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

#include "MeshFrameworkFactory.hh"
#include "MeshLogicalFactory.hh"
#include "MeshFrameworkColumn.hh"
#include "MeshSurfaceCell.hh"

#include "MeshFactory.hh"

namespace Amanzi {
namespace AmanziMesh {

Teuchos::RCP<Mesh>
MeshFactory::create(const std::string& filename)
{
  Teuchos::RCP<MeshFramework> mesh_fw = MeshFrameworkFactory::create(filename);
  Teuchos::RCP<Mesh> mesh = Teuchos::rcp(
    new Mesh(mesh_fw, Teuchos::rcp(new AmanziMesh::MeshAlgorithms()), Teuchos::null));

  if (plist_->get<bool>("build columns", false)) {
    mesh->buildColumns();
  } else if (plist_->isParameter("build columns from set")) {
    std::string setname = plist_->get<std::string>("build columns from set");
    mesh->buildColumns(std::vector<std::string>{ setname });
  }
  return mesh;
}


Teuchos::RCP<Mesh>
MeshFactory::create(const double x0,
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
  Teuchos::RCP<Mesh> mesh = Teuchos::rcp(
    new Mesh(mesh_fw, Teuchos::rcp(new AmanziMesh::MeshAlgorithms()), Teuchos::null));

  if (plist_->get<bool>("build columns", false)) {
    mesh->buildColumns();
  } else if (plist_->isParameter("build columns from set")) {
    std::string setname = plist_->get<std::string>("build columns from set");
    mesh->buildColumns(std::vector<std::string>{ setname });
  }
  return mesh;
}


Teuchos::RCP<Mesh>
MeshFactory::create(const double x0,
                    const double y0,
                    const double x1,
                    const double y1,
                    const int nx,
                    const int ny)
{
  Teuchos::RCP<MeshFramework> mesh_fw = MeshFrameworkFactory::create(x0, y0, x1, y1, nx, ny);
  Teuchos::RCP<Mesh> mesh = Teuchos::rcp(
    new Mesh(mesh_fw, Teuchos::rcp(new AmanziMesh::MeshAlgorithms()), Teuchos::null));

  if (plist_->get<bool>("build columns", false)) {
    mesh->buildColumns();
  } else if (plist_->isParameter("build columns from set")) {
    std::string setname = plist_->get<std::string>("build columns from set");
    mesh->buildColumns(std::vector<std::string>{ setname });
  }
  return mesh;
}


Teuchos::RCP<Mesh>
MeshFactory::create(const Teuchos::ParameterList& gen_plist)
{
  Teuchos::RCP<MeshFramework> mesh_fw = MeshFrameworkFactory::create(gen_plist);
  Teuchos::RCP<Mesh> mesh = Teuchos::rcp(
    new Mesh(mesh_fw, Teuchos::rcp(new AmanziMesh::MeshAlgorithms()), Teuchos::null));

  if (plist_->get<bool>("build columns", false)) {
    mesh->buildColumns();
  } else if (plist_->isParameter("build columns from set")) {
    std::string setname = plist_->get<std::string>("build columns from set");
    mesh->buildColumns(std::vector<std::string>{ setname });
  }
  return mesh;
}


Teuchos::RCP<Mesh>
MeshFactory::create(const Teuchos::RCP<const Mesh>& parent_mesh,
                    const MeshFramework::cEntity_ID_View& setids,
                    const Entity_kind setkind,
                    const bool flatten)
{
  if (parent_mesh->getMeshFramework() == Teuchos::null) {
    Errors::Message msg(
      "Cannot create an extracted mesh from a parent whose framework has been deleted.");
  }

  Teuchos::RCP<MeshFramework> mesh_fw =
    MeshFrameworkFactory::create(parent_mesh, setids, setkind, flatten);

  Teuchos::RCP<Mesh> mesh =
    Teuchos::rcp(new Mesh(mesh_fw, Teuchos::rcp(new MeshAlgorithms()), Teuchos::null));
  mesh->setParentMesh(parent_mesh);
  return mesh;
}


Teuchos::RCP<Mesh>
MeshFactory::create(const Teuchos::RCP<const Mesh>& parent_mesh,
                    const std::vector<std::string>& setnames,
                    const Entity_kind setkind,
                    const bool flatten)
{
  if (parent_mesh->getMeshFramework() == Teuchos::null) {
    Errors::Message msg(
      "Cannot create an extracted mesh from a parent whose framework has been deleted.");
  }

  Teuchos::RCP<MeshFramework> mesh_fw =
    MeshFrameworkFactory::create(parent_mesh, setnames, setkind, flatten);
  Teuchos::RCP<Mesh> mesh = Teuchos::null;
  if (mesh_fw != Teuchos::null) {
    mesh = Teuchos::rcp(new Mesh(mesh_fw, Teuchos::rcp(new MeshAlgorithms()), Teuchos::null));
    mesh->setParentMesh(parent_mesh);
  }
  return mesh;
}


Teuchos::RCP<Mesh>
MeshFactory::createLogical(Teuchos::ParameterList& log_plist)
{
  MeshLogicalFactory log_fac(comm_, gm_);
  Teuchos::RCP<MeshFramework> mesh_fw = log_fac.create(log_plist);
  Teuchos::RCP<Mesh> mesh = Teuchos::rcp(
    new Mesh(mesh_fw, Teuchos::rcp(new MeshLogicalAlgorithms()), Teuchos::null));
  return mesh;
}


// Create a 1D Column Mesh from a columnar structured volume mesh.
//
Teuchos::RCP<Mesh>
MeshFactory::createColumn(const Teuchos::RCP<Mesh>& parent_mesh,
                          int col_id,
                          const Teuchos::RCP<Teuchos::ParameterList>& plist)
{
  parent_mesh->buildColumns();

  // create a framework of the extracted 3D column
  Teuchos::RCP<MeshFramework> column_extracted_3D =
    MeshFrameworkFactory::create(parent_mesh,
                                 parent_mesh->columns->getCells<MemSpace_kind::HOST>(col_id),
                                 Entity_kind::CELL,
                                 false);

  // create the MeshFrameworkColumn
  Teuchos::RCP<MeshFramework> column_1D =
    Teuchos::rcp(new MeshFrameworkColumn(column_extracted_3D, plist_));

  // create and return the Mesh
  Teuchos::RCP<Mesh> mesh =
    Teuchos::rcp(new Mesh(column_1D, Teuchos::rcp(new MeshColumnAlgorithms()), plist));
  mesh->setParentMesh(parent_mesh);
  return mesh;
}

// Create a MeshSurfaceCell from a MeshFrameworkColumn
Teuchos::RCP<Mesh>
MeshFactory::createSurfaceCell(const Teuchos::RCP<const Mesh>& parent_mesh)
{
  if (parent_mesh->getMeshFramework() == Teuchos::null) {
    Errors::Message msg(
      "Cannot create a surface cell mesh from a column whose framework has been deleted.");
  }

  Teuchos::RCP<MeshFramework> mesh_surf_cell_fw =
    Teuchos::rcp(new MeshSurfaceCell(parent_mesh->getMeshFramework()));

  // create and return the Mesh
  Teuchos::RCP<Mesh> mesh = Teuchos::rcp(
    new Mesh(mesh_surf_cell_fw, Teuchos::rcp(new MeshAlgorithms()), Teuchos::null));
  mesh->setParentMesh(parent_mesh);
  return mesh;
}


} // namespace AmanziMesh
} // namespace Amanzi
