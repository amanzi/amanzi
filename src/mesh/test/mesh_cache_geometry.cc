/*
  Copyright 2010-201x held jointly by LANL, ORNL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

#include "UnitTest++.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"

#include "AmanziComm.hh"
#include "Mesh.hh"
#include "MeshAudit.hh"

#include "framework_meshes.hh"
#include "geometry_harnesses.hh"

//
// This tests the geometric routines in MeshCache, WITH caching, and having
// removed the framework mesh.  Otherwise it is identical to mesh_geometry.
//
// Effectively this means that all function calls that need a framework will
// use it, and the rest will call the algorithm.
//
TEST(MESH_CACHE_GEOMETRY_PLANAR)
{
  // a 2D, generated, structured quad on the unit square, NX=NY=2
  std::vector<Framework> frameworks;
  // works in MSTK
  if (framework_enabled(Framework::MSTK)) { frameworks.push_back(Framework::MSTK); }

  for (const auto& frm : frameworks) {
    std::cout << std::endl
              << "Testing 2D geometry with " << AmanziMesh::to_string(frm) << std::endl
              << "------------------------------------------------" << std::endl;

    auto mesh = createStructuredUnitQuad(Preference{ frm }, 2, 2);
    // cache, pitch the framework, repeat
    AmanziMesh::cacheAll(*mesh);
    mesh->destroyFramework();

    testMeshAudit<MeshAudit, Mesh>(mesh);
    testGeometryQuad(mesh, 2, 2);
    testExteriorMapsUnitBox(mesh, 2, 2);
  }
}


TEST(MESH_CACHE_GEOMETRY_1CUBE_GENERATED)
{
  // a 3D, generated, structured hex on the unit cube, NX=NY=NZ=1
  // only makes sense in serial
  if (getDefaultComm()->getSize() != 1) return;

  std::vector<Framework> frameworks;
  // works in MSTK & SIMPLE (in serial)
  if (framework_enabled(Framework::MSTK)) { frameworks.push_back(Framework::MSTK); }
  if (getDefaultComm()->getSize() == 1) frameworks.push_back(Framework::SIMPLE);

  for (const auto& frm : frameworks) {
    std::cout << std::endl
              << "Testing 3D geometry with " << AmanziMesh::to_string(frm) << std::endl
              << "------------------------------------------------" << std::endl;
    auto mesh = createStructuredUnitHex(Preference{ frm }, 1, 1, 1);
    // cache, pitch the framework, repeat
    AmanziMesh::cacheAll(*mesh);
    mesh->destroyFramework();

    testMeshAudit<MeshAudit, Mesh>(mesh);
    testGeometryCube(mesh, 1, 1, 1);
    if (frm != Framework::SIMPLE) testExteriorMapsUnitBox(mesh, 1, 1, 1);
  }
}


TEST(MESH_CACHE_GEOMETRY_1CUBE_EXO)
{
  // a 3D exodus file, structured hex on the unit cube, NX=NY=NZ=1
  // only makes sense in serial
  if (getDefaultComm()->getSize() != 1) return;

  // works in MSTK or MOAB
  std::vector<Framework> frameworks;
  if (framework_enabled(Framework::MSTK)) { frameworks.push_back(Framework::MSTK); }
  if (framework_enabled(Framework::MOAB)) { frameworks.push_back(Framework::MOAB); }

  for (const auto& frm : frameworks) {
    std::cout << std::endl
              << "Testing 3D Box 1x1x1 Exo geometry with " << AmanziMesh::to_string(frm)
              << std::endl
              << "------------------------------------------------" << std::endl;
    auto mesh = createUnstructured(Preference{ frm }, "test/hex_1x1x1_sets.exo");
    // cache, pitch the framework, repeat
    AmanziMesh::cacheAll(*mesh);
    mesh->destroyFramework();

    testMeshAudit<MeshAudit, Mesh>(mesh);
    testGeometryCube(mesh, 1, 1, 1);
    if (frm == Framework::MSTK) testExteriorMapsUnitBox(mesh, 1, 1, 1);
  }
}


TEST(MESH_CACHE_GEOMETRY_3CUBE)
{
  // a 3D, generated, structured hex on the unit cube, NX=NY=NZ=3
  // works in MSTK & SIMPLE (in serial)
  std::vector<Framework> frameworks;
  if (framework_enabled(Framework::MSTK)) { frameworks.push_back(Framework::MSTK); }
  if (getDefaultComm()->getSize() == 1) frameworks.push_back(Framework::SIMPLE);

  for (const auto& frm : frameworks) {
    std::cout << std::endl
              << "Testing 3D Box 3x3x3 with " << AmanziMesh::to_string(frm) << std::endl
              << "------------------------------------------------" << std::endl;
    auto mesh = createStructuredUnitHex(Preference{ frm }, 3, 3, 3);
    // cache, pitch the framework, repeat
    AmanziMesh::cacheAll(*mesh);
    mesh->destroyFramework();

    testMeshAudit<MeshAudit, Mesh>(mesh);
    testGeometryCube(mesh, 3, 3, 3);
    if (frm == Framework::MSTK) { testExteriorMapsUnitBox(mesh, 3, 3, 3); }
  }
}


TEST(MESH_CACHE_GEOMETRY_3CUBE_EXO)
{
  // a 3D exodus file, structured hex on the unit cube, NX=NY=NZ=3
  // works in MSTK or MOAB
  std::vector<Framework> frameworks;
  if (framework_enabled(Framework::MSTK)) { frameworks.push_back(Framework::MSTK); }
  if (framework_enabled(Framework::MOAB) && getDefaultComm()->getSize() == 1) {
    // moab only reads exo in serial, otherwise must read par
    frameworks.push_back(Framework::MOAB);
  }

  for (const auto& frm : frameworks) {
    std::cout << std::endl
              << "Testing 3D Box 3x3x3 Exo with " << AmanziMesh::to_string(frm) << std::endl
              << "------------------------------------------------" << std::endl;
    auto mesh = createUnstructured(Preference{ frm }, "test/hex_3x3x3_sets.exo");
    // cache, pitch the framework, repeat
    AmanziMesh::cacheAll(*mesh);
    mesh->destroyFramework();

    testMeshAudit<MeshAudit, Mesh>(mesh);
    testGeometryCube(mesh, 3, 3, 3);
    if (frm == Framework::MSTK) { testExteriorMapsUnitBox(mesh, 3, 3, 3); }
  }
}


//
// NOTE: this should work but is blocked by MSTK ticket #119, which is used for partitioning
//
// TEST(MESH_CACHE_GEOMETRY_3CUBE_PAR)
// {
//   // a 3D exodus file, structured hex on the unit cube, NX=NY=NZ=3, prepartitioned
//   // works in MSTK or MOAB
//   int nprocs = getDefaultComm()->getSize();
//   if (nprocs != 2) return;

//   std::vector<Framework> frameworks;
//   if (framework_enabled(Framework::MSTK)) {
//     frameworks.push_back(Framework::MSTK);
//   }
//   if (framework_enabled(Framework::MOAB)) {
//     frameworks.push_back(Framework::MOAB);
//   }

//   for (const auto& frm : frameworks) {
//     std::cout << std::endl
//               << "Testing 3D Box 3x3x3 Exo with " << AmanziMesh::to_string(frm) << std::endl
//               << "------------------------------------------------" << std::endl;
//     auto mesh = createUnstructured(Preference{frm}, "test/hex_3x3x3.par");
// // cache, pitch the framework, repeat
// AmanziMesh::cacheAll(*mesh);
// mesh->destroyFramework();
//     testMeshAudit<MeshAudit, Mesh>(mesh);
//     testGeometryCube(mesh,3,3,3);
//     if (frm == Framework::MSTK) {
//       testExteriorMapsUnitBox(mesh,3,3,3);
//     }
//   }
// }


TEST(MESH_CACHE_GEOMETRY_2x3CUBE)
{
  // a 3D, generated, structured hex on the unit cube, NX=NY=NZ=3
  // works in MSTK & SIMPLE (in serial)
  std::vector<Framework> frameworks;
  if (framework_enabled(Framework::MSTK)) { frameworks.push_back(Framework::MSTK); }
  if (getDefaultComm()->getSize() == 1) frameworks.push_back(Framework::SIMPLE);

  for (const auto& frm : frameworks) {
    std::cout << std::endl
              << "Testing 3D Box 2x2x3 with " << AmanziMesh::to_string(frm) << std::endl
              << "------------------------------------------------" << std::endl;
    auto mesh = createStructuredUnitHex(Preference{ frm }, 2, 2, 3);
    // cache, pitch the framework, repeat
    AmanziMesh::cacheAll(*mesh);
    mesh->destroyFramework();

    testMeshAudit<MeshAudit, Mesh>(mesh);
    testGeometryCube(mesh, 2, 2, 3);

    if (frm == Framework::MSTK) testExteriorMapsUnitBox(mesh, 2, 2, 3);
  }
}


TEST(MESH_CACHE_GEOMETRY_FRACTURE_EXO)
{
  // only works in MSTK or MOAB
  // Note this only checks the exo mesh, which does not have the fractures in it!
  // Actual fracture capability is checked later, but would depend on this test passing.
  std::vector<Framework> frameworks;
  if (framework_enabled(Framework::MSTK)) { frameworks.push_back(Framework::MSTK); }

  // Not sure what is up with this mesh, but MOAB errors trying to read it.
  // if (framework_enabled(Framework::MOAB)) {
  //   frameworks.push_back(Framework::MOAB);
  // }

  for (const auto& frm : frameworks) {
    std::cout << std::endl
              << "Testing 3D Fracture Exo with " << AmanziMesh::to_string(frm) << std::endl
              << "------------------------------------------------" << std::endl;
    auto mesh = createUnstructured(Preference{ frm }, "test/fractures.exo");
    // cache, pitch the framework, repeat
    AmanziMesh::cacheAll(*mesh);
    mesh->destroyFramework();
    testMeshAudit<MeshAudit, Mesh>(mesh);
  }
}


TEST(MESH_CACHE_GEOMETRY_PINCHOUTS)
{
  // only MSTK can handle this mesh -- it has degeneracies
  std::vector<Framework> frameworks;
  if (framework_enabled(Framework::MSTK)) { frameworks.push_back(Framework::MSTK); }

  for (const auto& frm : frameworks) {
    std::cout << std::endl
              << "Testing 3D Pinchout with " << AmanziMesh::to_string(frm) << std::endl
              << "------------------------------------------------" << std::endl;
    auto mesh = createUnstructured(Preference{ frm }, "test/test_pri_pinchout_mesh.exo");
    // cache, pitch the framework, repeat
    AmanziMesh::cacheAll(*mesh);
    mesh->destroyFramework();
    testMeshAudit<MeshAudit, Mesh>(mesh);
  }
}
