/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Rao Garimella, others
*/

#include <iostream>
#include "UnitTest++.h"
#include "Teuchos_ParameterList.hpp"

#include "AmanziComm.hh"
#include "MeshFramework.hh"
#include "MeshFrameworkFactory.hh"
#include "MeshFrameworkAudit.hh"
#include "MeshAudit.hh"
#include "MeshException.hh"

#include "framework_meshes.hh"
#include "geometry_harnesses.hh"

using namespace Amanzi;

template <class MeshAudit_type, class Mesh_type>
void
test2D(const Teuchos::RCP<Mesh_type>& mesh)
{
  std::cout << "Pre-deform mesh audit" << std::endl
            << "----------------------------------------" << std::endl;
  testMeshAudit<MeshAudit_type, Mesh_type>(mesh);

  std::cout << "Deform" << std::endl << "----------------------------------------" << std::endl;
  // Deform the mesh
  int nnodes =
    mesh->getNumEntities(AmanziMesh::Entity_kind::NODE, AmanziMesh::Parallel_kind::OWNED);
  typename Mesh_type::Entity_ID_View nodeids("nodeids", nnodes);
  typename Mesh_type::Point_View newpos("newpos", nnodes);

  for (int j = 0; j < nnodes; j++) {
    nodeids[j] = j;
    AmanziGeometry::Point oldcoord(2), newcoord(2);
    oldcoord = mesh->getNodeCoordinate(j);
    newcoord.set(oldcoord[0], 0.5 * oldcoord[1]);
    newpos[j] = newcoord;
  }

  deform(*mesh, nodeids, newpos);

  // If the deformation was successful, the cell volumes should be half
  // of what they were
  int ncells = mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::ALL);

  for (int j = 0; j < ncells; j++) {
    double volume = mesh->getCellVolume(j);
    CHECK_CLOSE(0.5, volume, 1.e-10);
  }
  std::cout << "Post-deform mesh audit" << std::endl
            << "----------------------------------------" << std::endl;
  testMeshAudit<MeshAudit_type, Mesh_type>(mesh);
}

template <class MeshAudit_type, class Mesh_type>
void
test3D(const Teuchos::RCP<Mesh_type>& mesh)
{
  std::cout << "Pre-deform mesh audit" << std::endl
            << "----------------------------------------" << std::endl;
  testMeshAudit<MeshAudit_type, Mesh_type>(mesh);

  std::cout << "Deform" << std::endl << "----------------------------------------" << std::endl;
  // Deform the mesh
  int nnodes =
    mesh->getNumEntities(AmanziMesh::Entity_kind::NODE, AmanziMesh::Parallel_kind::OWNED);
  typename Mesh_type::Entity_ID_View nodeids("nodesids", nnodes);
  typename Mesh_type::Point_View newpos("newpos", nnodes);

  for (int j = 0; j < nnodes; j++) {
    nodeids[j] = j;
    AmanziGeometry::Point oldcoord(3), newcoord(3);
    oldcoord = mesh->getNodeCoordinate(j);
    newcoord.set(oldcoord[0], oldcoord[1], oldcoord[2] / 2.0);
    newpos[j] = newcoord;
  }
  deform(*mesh, nodeids, newpos);

  // check geometry
  std::cout << "Post-deform mesh audit" << std::endl
            << "----------------------------------------" << std::endl;
  testMeshAudit<MeshAudit_type, Mesh_type>(mesh);
  std::cout << "Post-deform mesh geometry" << std::endl
            << "----------------------------------------" << std::endl;
  testGeometryCube<Mesh_type>(mesh, 3, 3, 3);
}


TEST(MESH_CACHED_DEFORM2D)
{
  auto comm = getDefaultComm();

  // We are not including MOAB or SIMPLE since they cannot generate in 2D
  std::vector<AmanziMesh::Framework> frameworks;
  std::vector<std::string> framework_names;

  if (AmanziMesh::framework_enabled(AmanziMesh::Framework::MSTK)) {
    frameworks.push_back(AmanziMesh::Framework::MSTK);
    framework_names.push_back("MSTK");
  }

  for (const auto& frm : frameworks) {
    // Set the framework
    std::cout << std::endl
              << "Testing deformation 2D with " << AmanziMesh::to_string(frm) << std::endl
              << "===========================================" << std::endl;

    // Create the mesh
    auto mesh =
      createStructuredUnitQuad({ frm }, 10, 10, comm, Teuchos::null, Teuchos::null, 10.0, 10.0);
    AmanziMesh::cacheAll(*mesh);

    // deform and test
    test2D<MeshAudit>(mesh);
  }
}

TEST(MESH_CACHED_GENERATED_DEFORM3D)
{
  auto comm = getDefaultComm();
  const int nproc(comm->NumProc());
  if (nproc != 1) {
    std::cout << "Parallel deformation not implemented" << std::endl;
    return;
  }

  std::vector<AmanziMesh::Framework> frameworks;
  std::vector<std::string> framework_names;

  frameworks.push_back(AmanziMesh::Framework::SIMPLE);
  framework_names.push_back("simple");

  if (AmanziMesh::framework_enabled(AmanziMesh::Framework::MSTK)) {
    frameworks.push_back(AmanziMesh::Framework::MSTK);
    framework_names.push_back("MSTK");
  }

  for (const auto& frm : frameworks) {
    std::cout << std::endl
              << "Testing deformation 3D with " << AmanziMesh::to_string(frm) << std::endl
              << "===========================================" << std::endl;

    // start with a mesh that will be deformed into the known mesh coordinates
    auto mesh =
      createStructuredUnitHex({ frm }, 3, 3, 3, comm, Teuchos::null, Teuchos::null, 1.0, 1.0, 2.0);
    AmanziMesh::cacheAll(*mesh);

    test3D<MeshAudit>(mesh);
  }
}
