/*
  Copyright 2010-202x held jointly by participating institutions.
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
#include "MeshAudit.hh"

#include "framework_meshes.hh"
#include "set_harnesses.hh"

TEST(MESH_SETS_3CUBE)
{
  // create the comm and gm
  auto comm = Amanzi::getDefaultComm();
  std::string infilename = "test/hex_3x3x3.xml";
  Teuchos::ParameterXMLFileReader xmlreader(infilename);
  Teuchos::ParameterList reg_spec(xmlreader.getParameters());
  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, reg_spec, *comm));

  // a 3D, generated, structured hex on the unit cube, NX=NY=NZ=3
  // works in MSTK & SIMPLE (in serial)
  std::vector<Framework> frameworks;
  if (framework_enabled(Framework::MSTK)) { frameworks.push_back(Framework::MSTK); }
  if (getDefaultComm()->NumProc() == 1) frameworks.push_back(Framework::SIMPLE);

  for (const auto& frm : frameworks) {
    std::cout << std::endl
              << "Testing 3D Box 3x3x3 with " << AmanziMesh::to_string(frm) << std::endl
              << "------------------------------------------------" << std::endl;
    auto mesh = createStructuredUnitHex(Preference{ frm }, 3, 3, 3, comm, gm);
    testHexMeshSets3x3x3(mesh, false, frm);
  }
}


TEST(MESH_SETS_3CUBE_EXO)
{
  // create the comm and gm
  auto comm = Amanzi::getDefaultComm();
  std::string infilename = "test/hex_3x3x3.xml";
  Teuchos::ParameterXMLFileReader xmlreader(infilename);
  Teuchos::ParameterList reg_spec(xmlreader.getParameters());
  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, reg_spec, *comm));

  // a 3D exodus file, structured hex on the unit cube, NX=NY=NZ=3
  // works in MSTK or MOAB
  std::vector<Framework> frameworks;
  if (framework_enabled(Framework::MSTK)) { frameworks.push_back(Framework::MSTK); }
  if (framework_enabled(Framework::MOAB) && getDefaultComm()->NumProc() == 1) {
    // moab only reads exo in serial, otherwise must read par
    frameworks.push_back(Framework::MOAB);
  }

  for (const auto& frm : frameworks) {
    std::cout << std::endl
              << "Testing 3D Box 3x3x3 Exo with " << AmanziMesh::to_string(frm) << std::endl
              << "------------------------------------------------" << std::endl;
    auto mesh = createUnstructured(Preference{ frm }, "test/hex_3x3x3_sets.exo", comm, gm);
    testHexMeshSets3x3x3(mesh, true, frm);
  }
}

//
// NOTE: this should work but is blocked by MSTK ticket #119, which is used for partitioning.
// TEST(MESH_SETS_3CUBE_PAR)
// {
//   // create the comm and gm
//   auto comm = Amanzi::getDefaultComm();
//   int nprocs = comm->NumProc();
//   if (nprocs != 2) return;

//   std::string infilename = "test/hex_3x3x3.xml";
//   Teuchos::ParameterXMLFileReader xmlreader(infilename);
//   Teuchos::ParameterList reg_spec(xmlreader.getParameters());
//   auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, reg_spec, *comm));

//   // a 3D exodus file, structured hex on the unit cube, NX=NY=NZ=3, prepartitioned
//   // works in MSTK or MOAB
//   std::vector<Framework> frameworks;
//   if (framework_enabled(Framework::MSTK)) {
//     frameworks.push_back(Framework::MSTK);
//   }
//   if (framework_enabled(Framework::MOAB)) {
//     frameworks.push_back(Framework::MOAB);
//   }

//   for (const auto& frm : frameworks) {
//     std::cout << std::endl
//               << "Testing 3D Box 3x3x3 Par with " << AmanziMesh::to_string(frm) << std::endl
//               << "------------------------------------------------" << std::endl;
//     auto mesh = createUnstructured(Preference{frm}, "test/hex_3x3x3_sets.par", comm, gm);
//     testHexMeshSets3x3x3(mesh, true, frm);
//   }
// }


TEST(MESH_SETS_3QUAD)
{
  // create the comm and gm
  auto comm = Amanzi::getDefaultComm();
  std::string infilename = "test/quad_3x3.xml";
  Teuchos::ParameterXMLFileReader xmlreader(infilename);
  Teuchos::ParameterList reg_spec(xmlreader.getParameters());
  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(2, reg_spec, *comm));

  // a 2D, generated, structured hex on the unit cube, NX=NY=3
  // works in MSTK
  std::vector<Framework> frameworks;
  if (framework_enabled(Framework::MSTK)) { frameworks.push_back(Framework::MSTK); }

  for (const auto& frm : frameworks) {
    std::cout << std::endl
              << "Testing 2D Box 3x3 with " << AmanziMesh::to_string(frm) << std::endl
              << "------------------------------------------------" << std::endl;
    auto mesh = createStructuredUnitQuad(Preference{ frm }, 3, 3, comm, gm);
    testQuadMeshSets3x3(mesh, false, frm, false);
  }
}
