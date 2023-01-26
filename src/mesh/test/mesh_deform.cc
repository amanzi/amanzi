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
  // Deform the mesh
  int nnodes =
    mesh->getNumEntities(AmanziMesh::Entity_kind::NODE, AmanziMesh::Parallel_type::OWNED);
  AmanziMesh::Entity_ID_List nodeids("nodeids", nnodes);
  AmanziMesh::Point_List newpos("newpos", nnodes);

  for (int j = 0; j < nnodes; j++) {
    nodeids[j] = j;
    AmanziGeometry::Point oldcoord(2), newcoord(2);
    oldcoord = mesh->getNodeCoordinate(j);
    newcoord.set(oldcoord[0], 0.5 * oldcoord[1]);
    newpos[j] = newcoord;
  }

  int ierr = MeshAlgorithms::deform(*mesh, nodeids, newpos);
  CHECK_EQUAL(ierr, 0);

  // If the deformation was successful, the cell volumes should be half
  // of what they were
  int ncells = mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::ALL);

  for (int j = 0; j < ncells; j++) {
    double volume = mesh->getCellVolume(j);
    CHECK_CLOSE(0.5, volume, 1.e-10);
  }
  testMeshAudit<MeshAudit_type, Mesh_type>(mesh);
}


template <class MeshAudit_type, class Mesh_type>
void
test3D(const Teuchos::RCP<Mesh_type>& mesh)
{
  // Deform the mesh
  int nnodes =
    mesh->getNumEntities(AmanziMesh::Entity_kind::NODE, AmanziMesh::Parallel_type::OWNED);
  AmanziMesh::Entity_ID_List nodeids("nodesids", nnodes);
  AmanziMesh::Point_List newpos("newpos", nnodes);

  for (int j = 0; j < nnodes; j++) {
    nodeids[j] = j;
    AmanziGeometry::Point oldcoord(3), newcoord(3);
    oldcoord = mesh->getNodeCoordinate(j);
    newcoord.set(oldcoord[0], oldcoord[1], oldcoord[2] / 2.0);
    newpos[j] = newcoord;
  }
  int ierr = MeshAlgorithms::deform(*mesh, nodeids, newpos);
  CHECK_EQUAL(ierr, 0);

  // check geometry
  testMeshAudit<MeshAudit_type, Mesh_type>(mesh);
  testGeometryCube<Mesh_type>(mesh, 3, 3, 3);
}


TEST(MESH_FRAMEWORK_DEFORM2D)
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
    std::cout << "Testing deformation with " << AmanziMesh::to_string(frm) << std::endl;

    // Create the mesh
    auto mesh = createFrameworkStructuredUnitQuad(
      { frm }, 10, 10, comm, Teuchos::null, Teuchos::null, 10.0, 10.0);

    // deform and test
    test2D<MeshFrameworkAudit>(mesh);
  }
}


TEST(MESH_FRAMEWORK_GENERATED_DEFORM3D)
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
    std::cout << "Testing deformation with " << AmanziMesh::to_string(frm) << std::endl;

    // start with a mesh that will be deformed into the known mesh coordinates
    auto mesh = createFrameworkStructuredUnitHex(
      { frm }, 3, 3, 3, comm, Teuchos::null, Teuchos::null, 1.0, 1.0, 2.0);

    test3D<MeshFrameworkAudit>(mesh);
  }
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
    std::cout << "Testing deformation with " << AmanziMesh::to_string(frm) << std::endl;

    // Create the mesh
    auto mesh =
      createStructuredUnitQuad({ frm }, 10, 10, comm, Teuchos::null, Teuchos::null, 10.0, 10.0);

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
    std::cout << "Testing deformation with " << AmanziMesh::to_string(frm) << std::endl;

    // start with a mesh that will be deformed into the known mesh coordinates
    auto mesh =
      createStructuredUnitHex({ frm }, 3, 3, 3, comm, Teuchos::null, Teuchos::null, 1.0, 1.0, 2.0);

    test3D<MeshAudit>(mesh);
  }
}
