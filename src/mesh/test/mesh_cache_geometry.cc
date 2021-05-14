/*
  Copyright 2010-201x held jointly by LANL, ORNL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

//
// Tests mesh cache on a variety of meshes, making sure that the cache's
// geometry and topological information is set up correctly.
//

#include <UnitTest++.h>

#include "AmanziComm.hh"
#include "MeshFramework.hh"
#include "MeshFrameworkFactory.hh"
#include "MeshAudit_decl.hh"
#include "MeshAudit_impl.hh"

#include "framework_meshes.hh"
#include "cache_meshes.hh"
#include "geometry_harnesses.hh"

using MeshAudit = MeshAudit_<MeshCache, Impl::MeshAudit_Geometry>;

TEST(MESH_CACHE_GEOMETRY_PLANAR)
{
  // only works in MSTK
  std::vector<Framework> frameworks;
  if (framework_enabled(Framework::MSTK)) {
    frameworks.push_back(Framework::MSTK);
  }

  for (const auto& frm : frameworks) {
    std::cout << std::endl
              << "Testing 2D geometry with " << AmanziMesh::framework_names.at(frm) << std::endl
              << "------------------------------------------------" << std::endl;

    auto mesh_fw = createFrameworkStructuredUnitSquare(Preference{frm}, 2);
    auto mesh = createMeshFromFramework(mesh_fw);
    testGeometry2x2<MeshAudit, MeshCache>(mesh);
  }
}


TEST(MESH_CACHE_GEOMETRY_1BOX_GENERATED)
{
  // only works in MSTK
  std::vector<Framework> frameworks;
  if (framework_enabled(Framework::MSTK)) {
    frameworks.push_back(Framework::MSTK);
  }

  if (getDefaultComm()->NumProc() == 1)
    frameworks.push_back(Framework::SIMPLE);

  for (const auto& frm : frameworks) {
    std::cout << std::endl
              << "Testing 3D geometry with " << AmanziMesh::framework_names.at(frm) << std::endl
              << "------------------------------------------------" << std::endl;
    auto mesh_fw = createFrameworkStructuredUnitCube(Preference{frm}, 1);
    auto mesh = createMeshFromFramework(mesh_fw);
    testGeometry1x1x1<MeshAudit, MeshCache>(mesh);
  }
}


TEST(MESH_CACHE_GEOMETRY_1BOX_EXO)
{
  // only works in MSTK or MOAB
  std::vector<Framework> frameworks;
  if (framework_enabled(Framework::MSTK)) {
    frameworks.push_back(Framework::MSTK);
  }
  if (framework_enabled(Framework::MOAB)) {
    frameworks.push_back(Framework::MOAB);
  }

  for (const auto& frm : frameworks) {
    std::cout << std::endl
              << "Testing 3D Box 1x1x1 Exo geometry with " << AmanziMesh::framework_names.at(frm) << std::endl
              << "------------------------------------------------" << std::endl;
    auto mesh_fw = createFrameworkUnstructured(Preference{frm}, "test/hex_1x1x1_sets.exo");
    auto mesh = createMeshFromFramework(mesh_fw);
    testGeometry1x1x1<MeshAudit, MeshCache>(mesh);
  }
}


TEST(MESH_CACHE_GEOMETRY_2BOX)
{
  // only works in MSTK or Simple
  std::vector<Framework> frameworks;
  if (framework_enabled(Framework::MSTK)) {
    frameworks.push_back(Framework::MSTK);
  }

  if (getDefaultComm()->NumProc() == 1)
    frameworks.push_back(Framework::SIMPLE);

  for (const auto& frm : frameworks) {
    std::cout << std::endl
              << "Testing 3D Box 2x2x2 with " << AmanziMesh::framework_names.at(frm) << std::endl
              << "------------------------------------------------" << std::endl;
    auto mesh_fw = createFrameworkStructuredUnitCube(Preference{frm}, 2);
    auto mesh = createMeshFromFramework(mesh_fw);
    testGeometry2x2x2<MeshAudit, MeshCache>(mesh);
  }
}


TEST(MESH_CACHE_GEOMETRY_3BOX_EXO)
{
  // only works in MSTK or MOAB
  std::vector<Framework> frameworks;
  if (framework_enabled(Framework::MSTK)) {
    frameworks.push_back(Framework::MSTK);
  }
  if (framework_enabled(Framework::MOAB)) {
    frameworks.push_back(Framework::MOAB);
  }

  for (const auto& frm : frameworks) {
    std::cout << std::endl
              << "Testing 3D Box 3x3x3 Exo with " << AmanziMesh::framework_names.at(frm) << std::endl
              << "------------------------------------------------" << std::endl;
    auto mesh_fw = createFrameworkUnstructured(Preference{frm}, "test/hex_3x3x3_sets.exo");
    auto mesh = createMeshFromFramework(mesh_fw);
    testMeshAudit<MeshAudit, MeshCache>(mesh);
  }
}


TEST(MESH_CACHE_GEOMETRY_FRACTURE_EXO)
{
  // only works in MSTK or MOAB
  std::vector<Framework> frameworks;
  if (framework_enabled(Framework::MSTK)) {
    frameworks.push_back(Framework::MSTK);
  }
  if (framework_enabled(Framework::MOAB)) {
    frameworks.push_back(Framework::MOAB);
  }

  for (const auto& frm : frameworks) {
    std::cout << std::endl
              << "Testing 3D Fracture Exo with " << AmanziMesh::framework_names.at(frm) << std::endl
              << "------------------------------------------------" << std::endl;
    auto mesh_fw = createFrameworkUnstructured(Preference{frm}, "test/fractures.exo");
    auto mesh = createMeshFromFramework(mesh_fw);
    testMeshAudit<MeshAudit, MeshCache>(mesh);
  }
}


