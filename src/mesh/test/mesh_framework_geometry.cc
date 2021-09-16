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
#include "MeshFramework.hh"
#include "MeshFrameworkFactory.hh"
#include "MeshFrameworkAudit.hh"

#include "framework_meshes.hh"
#include "geometry_harnesses.hh"

TEST(MESH_GEOMETRY_PLANAR)
{
  // a 2D, generated, structured quad on the unit square, NX=NY=2
  std::vector<Framework> frameworks;
  // works in MSTK
  if (framework_enabled(Framework::MSTK)) {
    frameworks.push_back(Framework::MSTK);
  }

  for (const auto& frm : frameworks) {
    std::cout << std::endl
              << "Testing 2D geometry with " << AmanziMesh::framework_names.at(frm) << std::endl
              << "------------------------------------------------" << std::endl;

    auto mesh = createFrameworkStructuredUnitQuad(Preference{frm}, 2,2);
    testMeshAudit<MeshFrameworkAudit, MeshFramework>(mesh);
    testGeometryQuad(mesh,2,2);
  }
}


TEST(MESH_GEOMETRY_1CUBE_GENERATED)
{
  // a 3D, generated, structured hex on the unit cube, NX=NY=NZ=1
  // only makes sense in serial
  if (getDefaultComm()->NumProc() != 1) return;

  std::vector<Framework> frameworks;
  // works in MSTK & SIMPLE (in serial)
  if (framework_enabled(Framework::MSTK)) {
    frameworks.push_back(Framework::MSTK);
  }
  if (getDefaultComm()->NumProc() == 1)
    frameworks.push_back(Framework::SIMPLE);

  for (const auto& frm : frameworks) {
    std::cout << std::endl
              << "Testing 3D geometry with " << AmanziMesh::framework_names.at(frm) << std::endl
              << "------------------------------------------------" << std::endl;
    auto mesh = createFrameworkStructuredUnitHex(Preference{frm}, 1,1,1);
    testMeshAudit<MeshFrameworkAudit, MeshFramework>(mesh);
    testGeometryCube(mesh,1,1,1);
  }
}


TEST(MESH_GEOMETRY_1CUBE_EXO)
{
  // a 3D exodus file, structured hex on the unit cube, NX=NY=NZ=1
  // only makes sense in serial
  if (getDefaultComm()->NumProc() != 1) return;

  // works in MSTK or MOAB
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
    auto mesh = createFrameworkUnstructured(Preference{frm}, "test/hex_1x1x1_sets.exo");
    testMeshAudit<MeshFrameworkAudit, MeshFramework>(mesh);
    testGeometryCube(mesh,1,1,1);
  }
}


TEST(MESH_GEOMETRY_3CUBE)
{
  // a 3D, generated, structured hex on the unit cube, NX=NY=NZ=3
  // works in MSTK & SIMPLE (in serial)
  std::vector<Framework> frameworks;
  if (framework_enabled(Framework::MSTK)) {
    frameworks.push_back(Framework::MSTK);
  }
  if (getDefaultComm()->NumProc() == 1)
    frameworks.push_back(Framework::SIMPLE);

  for (const auto& frm : frameworks) {
    std::cout << std::endl
              << "Testing 3D Box 3x3x3 with " << AmanziMesh::framework_names.at(frm) << std::endl
              << "------------------------------------------------" << std::endl;
    auto mesh = createFrameworkStructuredUnitHex(Preference{frm}, 3,3,3);
    testMeshAudit<MeshFrameworkAudit, MeshFramework>(mesh);
    testGeometryCube(mesh,3,3,3);
  }
}


TEST(MESH_GEOMETRY_3CUBE_EXO)
{
  // a 3D exodus file, structured hex on the unit cube, NX=NY=NZ=3
  // works in MSTK or MOAB
  std::vector<Framework> frameworks;
  if (framework_enabled(Framework::MSTK)) {
    frameworks.push_back(Framework::MSTK);
  }
  if (framework_enabled(Framework::MOAB) &&
      getDefaultComm()->NumProc() == 1) {
    // moab only reads exo in serial, otherwise must read par
    frameworks.push_back(Framework::MOAB);
  }

  for (const auto& frm : frameworks) {
    std::cout << std::endl
              << "Testing 3D Box 3x3x3 Exo with " << AmanziMesh::framework_names.at(frm) << std::endl
              << "------------------------------------------------" << std::endl;
    auto mesh = createFrameworkUnstructured(Preference{frm}, "test/hex_3x3x3_sets.exo");
    testMeshAudit<MeshFrameworkAudit, MeshFramework>(mesh);
    testGeometryCube(mesh,3,3,3);
  }
}


//
// NOTE: this should work but is blocked by MSTK ticket #119, which is used for partitioning
//
// TEST(MESH_GEOMETRY_3CUBE_PAR)
// {
//   // a 3D exodus file, structured hex on the unit cube, NX=NY=NZ=3, prepartitioned
//   // works in MSTK or MOAB
//   int nprocs = getDefaultComm()->NumProc();
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
//               << "Testing 3D Box 3x3x3 Exo with " << AmanziMesh::framework_names.at(frm) << std::endl
//               << "------------------------------------------------" << std::endl;
//     auto mesh = createFrameworkUnstructured(Preference{frm}, "test/hex_3x3x3.par");
//     testMeshAudit<MeshFrameworkAudit, MeshFramework>(mesh);
//     testGeometryCube(mesh,3,3,3);
//   }
// }



TEST(MESH_GEOMETRY_2x3CUBE)
{
  // a 3D, generated, structured hex on the unit cube, NX=NY=NZ=3
  // works in MSTK & SIMPLE (in serial)
  std::vector<Framework> frameworks;
  if (framework_enabled(Framework::MSTK)) {
    frameworks.push_back(Framework::MSTK);
  }
  if (getDefaultComm()->NumProc() == 1)
    frameworks.push_back(Framework::SIMPLE);

  for (const auto& frm : frameworks) {
    std::cout << std::endl
              << "Testing 3D Box 2x2x3 with " << AmanziMesh::framework_names.at(frm) << std::endl
              << "------------------------------------------------" << std::endl;
    auto mesh = createFrameworkStructuredUnitHex(Preference{frm}, 2,2,3);
    testMeshAudit<MeshFrameworkAudit, MeshFramework>(mesh);
    testGeometryCube(mesh,2,2,3);
  }
}


TEST(MESH_GEOMETRY_FRACTURE_EXO)
{
  // only works in MSTK or MOAB
  // Note this only checks the exo mesh, which does not have the fractures in it!
  // Actual fracture capability is checked later, but would depend on this test passing.
  std::vector<Framework> frameworks;
  if (framework_enabled(Framework::MSTK)) {
    frameworks.push_back(Framework::MSTK);
  }

  // Not sure what is up with this mesh, but MOAB errors trying to read it.
  // if (framework_enabled(Framework::MOAB)) {
  //   frameworks.push_back(Framework::MOAB);
  // }

  for (const auto& frm : frameworks) {
    std::cout << std::endl
              << "Testing 3D Fracture Exo with " << AmanziMesh::framework_names.at(frm) << std::endl
              << "------------------------------------------------" << std::endl;
    auto mesh = createFrameworkUnstructured(Preference{frm}, "test/fractures.exo");
    testMeshAudit<MeshFrameworkAudit, MeshFramework>(mesh);
  }
}


TEST(MESH_GEOMETRY_PINCHOUTS)
{
  // only MSTK can handle this mesh -- it has degeneracies
  std::vector<Framework> frameworks;
  if (framework_enabled(Framework::MSTK)) {
    frameworks.push_back(Framework::MSTK);
  }

  for (const auto& frm : frameworks) {
    std::cout << std::endl
              << "Testing 3D Pinchout with " << AmanziMesh::framework_names.at(frm) << std::endl
              << "------------------------------------------------" << std::endl;
    auto mesh = createFrameworkUnstructured(Preference{frm}, "test/test_pri_pinchout_mesh.exo");
    testMeshAudit<MeshFrameworkAudit, MeshFramework>(mesh);
  }
}

