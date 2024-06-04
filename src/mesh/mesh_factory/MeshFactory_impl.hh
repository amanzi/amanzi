/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

#pragma once

#include "MeshFrameworkFactory.hh"
#include "MeshLogicalFactory.hh"
#include "MeshFrameworkColumn.hh"
#include "MeshSurfaceCell.hh"

#include "MeshFactory_decl.hh"

namespace Amanzi {
namespace AmanziMesh {

template <MemSpace_kind MEM>
auto
MeshFactory::create(const std::string& filename)
{
  Teuchos::RCP<MeshFramework> mesh_fw = MeshFrameworkFactory::create(filename);
  Teuchos::RCP<MeshCacheHost> mesh_on_host = Teuchos::rcp(
    new MeshCacheHost(mesh_fw, Teuchos::rcp(new AmanziMesh::MeshAlgorithms()), Teuchos::null));

  if (plist_->get<bool>("build columns", false)) {
    mesh_on_host->buildColumns();
  } else if (plist_->isParameter("build columns from set")) {
    std::string setname = plist_->get<std::string>("build columns from set");
    mesh_on_host->buildColumns(std::vector<std::string>{ setname });
  }

  if constexpr (MEM == MemSpace_kind::HOST) {
    return mesh_on_host;
  } else {
    auto mesh = Teuchos::rcp(new Mesh(*mesh_on_host));
    return mesh;
  }
}


template <MemSpace_kind MEM>
auto
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
  Teuchos::RCP<MeshCacheHost> mesh_on_host = Teuchos::rcp(
    new MeshCacheHost(mesh_fw, Teuchos::rcp(new AmanziMesh::MeshAlgorithms()), Teuchos::null));

  if (plist_->get<bool>("build columns", false)) {
    mesh_on_host->buildColumns();
  } else if (plist_->isParameter("build columns from set")) {
    std::string setname = plist_->get<std::string>("build columns from set");
    mesh_on_host->buildColumns(std::vector<std::string>{ setname });
  }

  if constexpr (MEM == MemSpace_kind::HOST) {
    return mesh_on_host;
  } else {
    auto mesh = Teuchos::rcp(new Mesh(*mesh_on_host));
    return mesh;
  }
}


template <MemSpace_kind MEM>
auto
MeshFactory::create(const double x0,
                    const double y0,
                    const double x1,
                    const double y1,
                    const int nx,
                    const int ny)
{
  Teuchos::RCP<MeshFramework> mesh_fw = MeshFrameworkFactory::create(x0, y0, x1, y1, nx, ny);
  Teuchos::RCP<MeshCacheHost> mesh_on_host = Teuchos::rcp(
    new MeshCacheHost(mesh_fw, Teuchos::rcp(new AmanziMesh::MeshAlgorithms()), Teuchos::null));

  if (plist_->get<bool>("build columns", false)) {
    mesh_on_host->buildColumns();
  } else if (plist_->isParameter("build columns from set")) {
    std::string setname = plist_->get<std::string>("build columns from set");
    mesh_on_host->buildColumns(std::vector<std::string>{ setname });
  }

  if constexpr (MEM == MemSpace_kind::HOST) {
    return mesh_on_host;
  } else {
    auto mesh = Teuchos::rcp(new Mesh(*mesh_on_host));
    return mesh;
  }
}


template <MemSpace_kind MEM>
auto
MeshFactory::create(const Teuchos::ParameterList& gen_plist)
{
  Teuchos::RCP<MeshFramework> mesh_fw = MeshFrameworkFactory::create(gen_plist);
  Teuchos::RCP<MeshCacheHost> mesh_on_host = Teuchos::rcp(
    new MeshCacheHost(mesh_fw, Teuchos::rcp(new AmanziMesh::MeshAlgorithms()), Teuchos::null));

  if (plist_->get<bool>("build columns", false)) {
    mesh_on_host->buildColumns();
  } else if (plist_->isParameter("build columns from set")) {
    std::string setname = plist_->get<std::string>("build columns from set");
    mesh_on_host->buildColumns(std::vector<std::string>{ setname });
  }

  if constexpr (MEM == MemSpace_kind::HOST) {
    return mesh_on_host;
  } else {
    auto mesh = Teuchos::rcp(new Mesh(*mesh_on_host));
    return mesh;
  }
}


template <typename Mesh_type>
Teuchos::RCP<Mesh_type>
MeshFactory::create(const Teuchos::RCP<const Mesh_type>& parent_mesh,
                    const MeshFramework::cEntity_ID_View& setids,
                    const Entity_kind setkind,
                    const bool flatten)
{
  if (parent_mesh->getMeshFramework() == Teuchos::null) {
    Errors::Message msg(
      "Cannot create an extracted mesh from a parent whose framework has been deleted.");
  }

  Teuchos::RCP<const MeshHost> parent_on_host;
  if constexpr (std::is_same_v<Mesh_type, MeshHost>) {
    parent_on_host = parent_mesh;
  } else {
    parent_on_host = onMemHost(parent_mesh);
  }

  Teuchos::RCP<MeshFramework> mesh_fw =
    MeshFrameworkFactory::create(parent_on_host, setids, setkind, flatten);

  auto mesh_on_host =
    Teuchos::rcp(new MeshCacheHost(mesh_fw, Teuchos::rcp(new MeshAlgorithms()), Teuchos::null));
  mesh_on_host->setParentMesh(parent_on_host);

  if constexpr (std::is_same_v<Mesh_type, MeshHost>) {
    return mesh_on_host;
  } else {
    auto mesh = Teuchos::rcp(new Mesh(*mesh_on_host));
    mesh->setParentMesh(parent_mesh);
    return mesh;
  }
}


template <typename Mesh_type>
Teuchos::RCP<Mesh_type>
MeshFactory::create(const Teuchos::RCP<const Mesh_type>& parent_mesh,
                    const std::vector<std::string>& setnames,
                    const Entity_kind setkind,
                    const bool flatten)
{
  if (parent_mesh->getMeshFramework() == Teuchos::null) {
    Errors::Message msg(
      "Cannot create an extracted mesh from a parent whose framework has been deleted.");
  }

  Teuchos::RCP<const MeshHost> parent_on_host;
  if constexpr (std::is_same_v<Mesh_type, MeshHost>) {
    parent_on_host = parent_mesh;
  } else {
    parent_on_host = onMemHost(parent_mesh);
  }

  Teuchos::RCP<MeshFramework> mesh_fw =
    MeshFrameworkFactory::create(parent_on_host, setnames, setkind, flatten);

  Teuchos::RCP<MeshHost> mesh_on_host = Teuchos::null;
  if (mesh_fw != Teuchos::null) {
    mesh_on_host =
      Teuchos::rcp(new MeshCacheHost(mesh_fw, Teuchos::rcp(new MeshAlgorithms()), Teuchos::null));
    mesh_on_host->setParentMesh(parent_on_host);
  }

  if constexpr (std::is_same_v<Mesh_type, MeshHost>) {
    return mesh_on_host;
  } else {
    Teuchos::RCP<Mesh_type> mesh;
    if (mesh_on_host != Teuchos::null) {
      mesh = Teuchos::rcp(new Mesh(*mesh_on_host));
      mesh->setParentMesh(parent_mesh);
    }
    return mesh;
  }
}


template <MemSpace_kind MEM>
auto
MeshFactory::createLogical(Teuchos::ParameterList& log_plist)
{
  MeshLogicalFactory log_fac(comm_, gm_);
  Teuchos::RCP<MeshFramework> mesh_fw = log_fac.create(log_plist);
  auto mesh_on_host = Teuchos::rcp(
    new MeshCacheHost(mesh_fw, Teuchos::rcp(new MeshLogicalAlgorithms()), Teuchos::null));

  if constexpr (MEM == MemSpace_kind::HOST) {
    return mesh_on_host;
  } else {
    auto mesh = Teuchos::rcp(new MeshCacheDevice(*mesh_on_host));
    return mesh;
  }
}


// Create a 1D Column Mesh from a columnar structured volume mesh.
//
template <typename Mesh_type>
Teuchos::RCP<Mesh_type>
MeshFactory::createColumn(const Teuchos::RCP<Mesh_type>& parent_mesh,
                          int col_id,
                          const Teuchos::RCP<Teuchos::ParameterList>& plist)
{
  // create a framework of the extracted 3D column
  Teuchos::RCP<MeshHost> parent_on_host;
  if constexpr (std::is_same_v<Mesh_type, MeshCacheHost>) {
    parent_on_host = parent_mesh;
  } else {
    parent_on_host = onMemHost(parent_mesh);
  }
  parent_on_host->buildColumns();

  Teuchos::RCP<MeshFramework> column_extracted_3D =
    MeshFrameworkFactory::create(parent_on_host,
                                 parent_on_host->columns.getCells<MemSpace_kind::HOST>(col_id),
                                 Entity_kind::CELL,
                                 false);

  // create the MeshFrameworkColumn
  Teuchos::RCP<MeshFramework> column_1D =
    Teuchos::rcp(new MeshFrameworkColumn(column_extracted_3D, plist_));

  // create and return the Mesh
  auto mesh_on_host =
    Teuchos::rcp(new MeshCacheHost(column_1D, Teuchos::rcp(new MeshColumnAlgorithms()), plist));
  mesh_on_host->setParentMesh(parent_on_host);

  if constexpr (std::is_same_v<Mesh_type, MeshHost>) {
    return mesh_on_host;
  } else {
    auto mesh = Teuchos::rcp(new Mesh(*mesh_on_host));
    mesh->setParentMesh(parent_mesh);
    return mesh;
  }
}

// Create a MeshSurfaceCell from a MeshFrameworkColumn
template <typename Mesh_type>
Teuchos::RCP<Mesh_type>
MeshFactory::createSurfaceCell(const Teuchos::RCP<const Mesh_type>& parent_mesh)
{
  if (parent_mesh->getMeshFramework() == Teuchos::null) {
    Errors::Message msg(
      "Cannot create a surface cell mesh from a column whose framework has been deleted.");
  }

  Teuchos::RCP<const MeshHost> parent_on_host;
  if constexpr (std::is_same_v<Mesh_type, MeshHost>) {
    parent_on_host = parent_mesh;
  } else {
    parent_on_host = onMemHost(parent_mesh);
  }

  Teuchos::RCP<MeshFramework> mesh_surf_cell_fw =
    Teuchos::rcp(new MeshSurfaceCell(parent_on_host->getMeshFramework()));

  // create and return the Mesh
  auto mesh_on_host = Teuchos::rcp(
    new MeshCacheHost(mesh_surf_cell_fw, Teuchos::rcp(new MeshAlgorithms()), Teuchos::null));
  mesh_on_host->setParentMesh(parent_on_host);

  if constexpr (std::is_same_v<Mesh_type, MeshHost>) {
    return mesh_on_host;
  } else {
    auto mesh = Teuchos::rcp(new Mesh(*mesh_on_host));
    mesh->setParentMesh(parent_mesh);
    return mesh;
  }
}


} // namespace AmanziMesh
} // namespace Amanzi
