/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

/*
  This test extracts two internal faces sets representing fractures and creates
  a submanifold mesh from these faces.
*/

#include "UnitTest++.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"

#include "AmanziComm.hh"
#include "RegionFactory.hh"
#include "Mesh.hh"
#include "MeshFramework.hh"
#include "MeshFactory.hh"
#include "MeshAudit.hh"

#include "framework_meshes.hh"
#include "geometry_harnesses.hh"

TEST(MESH_FRACTURE_EXTRACTION_GENERATED)
{
  // create the comm and gm
  auto comm = Amanzi::getDefaultComm();
  std::string infilename = "test/fracture.xml";
  Teuchos::ParameterXMLFileReader xmlreader(infilename);
  Teuchos::ParameterList reg_spec(xmlreader.getParameters());
  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, reg_spec, *comm));

  // a 3D, generated, structured hex on the unit cube, NX=NY=NZ=10
  // works in MSTK
  //
  // extract two inner surfaces as fractures
  std::vector<Framework> frameworks;
  if (framework_enabled(Framework::MSTK)) { frameworks.push_back(Framework::MSTK); }

  for (const auto& frm : frameworks) {
    std::cout << std::endl
              << "Testing Fracture Extraction with " << AmanziMesh::to_string(frm) << std::endl
              << "------------------------------------------------" << std::endl;
    auto parent_mesh = createFrameworkStructuredUnitHex(Preference{ frm }, 10, 10, 10, comm, gm);

    // extract the fractures
    std::vector<std::string> setnames{ "fracture 1", "fracture 2" };
    MeshFactory fac(comm, gm);
    fac.set_preference({ frm });
    // Make cache of current mesh
    auto parent_mesh_cache = Teuchos::rcp(
      new Mesh(parent_mesh, Teuchos::rcp(new AmanziMesh::MeshAlgorithms()), Teuchos::null));
    cacheAll(*parent_mesh_cache);

    auto mesh = fac.create(parent_mesh_cache, setnames, AmanziMesh::Entity_kind::FACE, false);
    cacheAll(*mesh);

    // test the surface mesh as a fracture mesh
    // -- mesh audit
    testMeshAudit<MeshAudit, Mesh>(mesh);

    // check we found the right number of things
    int ncells =
      mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
    CHECK_CLOSE_SUMALL(10 * 10 * 2, ncells, *comm);
    int nnodes =
      mesh->getNumEntities(AmanziMesh::Entity_kind::NODE, AmanziMesh::Parallel_kind::OWNED);
    CHECK_CLOSE_SUMALL(11 * 11 * 2 - 11, nnodes, *comm);
    int nfaces =
      mesh->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED);
    CHECK_CLOSE_SUMALL(10 * 11 * 4 - 10, nfaces, *comm);
    auto ents = mesh->getSetEntities(
      "fracture 1", AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
    CHECK_EQUAL(100, ents.size());
  }
}

TEST(MESH_FRACTURE_EXTRACTION_GENERATED_EXTRACTED_MANIFOLD)
{
  // This test is the same as the above, but as it works with the
  // MeshExtractedManifold mesh instead of the framework's extraction
  // constructor, it can be used with any mesh framework that supports edges.

  // create the comm and gm
  auto comm = Amanzi::getDefaultComm();
  std::string infilename = "test/fracture.xml";
  Teuchos::ParameterXMLFileReader xmlreader(infilename);
  Teuchos::ParameterList reg_spec(xmlreader.getParameters());
  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, reg_spec, *comm));

  // a 3D, generated, structured hex on the unit cube, NX=NY=NZ=10
  // works in MSTK
  //
  // extract two inner surfaces as fractures
  std::vector<Framework> frameworks;
  //if (comm->getSize() == 1) {
  //  frameworks.push_back(Framework::SIMPLE);
  //}
  if (framework_enabled(Framework::MSTK)) { frameworks.push_back(Framework::MSTK); }

  for (const auto& frm : frameworks) {
    std::cout << std::endl
              << "Testing Fracture Extraction with " << AmanziMesh::to_string(frm) << std::endl
              << "------------------------------------------------" << std::endl;
    auto fac_list = Teuchos::rcp(new Teuchos::ParameterList());
    fac_list->set("request edges", true);
    auto parent_mesh =
      createFrameworkStructuredUnitHex(Preference{ frm }, 10, 10, 10, comm, gm, fac_list);

    // extract the fractures
    auto fac_plist = Teuchos::rcp(new Teuchos::ParameterList());

    // Extracting multiple sets not supported by MeshExtractedManifold, see #596 part 4
    // std::vector<std::string> setnames{"fracture 1", "fracture 2"};
    std::vector<std::string> setnames{ "fractures" };
    fac_plist->sublist("unstructured")
      .sublist("submesh")
      .set<std::string>("extraction method", "manifold mesh");
    MeshFactory fac(comm, gm, fac_plist);
    fac.set_preference({ frm });
    // Make cache of current mesh
    auto parent_mesh_cache = Teuchos::rcp(
      new Mesh(parent_mesh, Teuchos::rcp(new AmanziMesh::MeshAlgorithms()), Teuchos::null));
    cacheAll(*parent_mesh_cache);

    auto mesh = fac.create(parent_mesh_cache, setnames, AmanziMesh::Entity_kind::FACE, false);
    cacheAll(*mesh);

    // test the surface mesh as a fracture mesh
    // -- mesh audit
    testMeshAudit<MeshAudit, Mesh>(mesh);

    // check we found the right number of things
    int ncells =
      mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
    CHECK_CLOSE_SUMALL(10 * 10 * 2, ncells, *comm);
    int nnodes =
      mesh->getNumEntities(AmanziMesh::Entity_kind::NODE, AmanziMesh::Parallel_kind::OWNED);
    CHECK_CLOSE_SUMALL(11 * 11 * 2 - 11, nnodes, *comm);
    int nfaces =
      mesh->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED);
    CHECK_CLOSE_SUMALL(10 * 11 * 4 - 10, nfaces, *comm);
    auto ents = mesh->getSetEntities(
      "fractures", AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
    CHECK_EQUAL(200, ents.size());
  }
}

TEST(MESH_FRACTURE_EXTRACTION_EXO)
{
  // create the comm and gm
  auto comm = Amanzi::getDefaultComm();
  std::string infilename = "test/fracture2.xml";
  Teuchos::ParameterXMLFileReader xmlreader(infilename);
  Teuchos::ParameterList reg_spec(xmlreader.getParameters());
  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, reg_spec, *comm));

  // a 3D, generated, structured hex on the unit cube, NX=NY=NZ=10
  // works in MSTK
  //
  // extract two inner surfaces as fractures
  std::vector<Framework> frameworks;
  if (framework_enabled(Framework::MSTK)) { frameworks.push_back(Framework::MSTK); }

  for (const auto& frm : frameworks) {
    std::cout << std::endl
              << "Testing Fracture Extraction with " << AmanziMesh::to_string(frm) << std::endl
              << "------------------------------------------------" << std::endl;
    auto parent_mesh =
      createFrameworkUnstructured(Preference{ frm }, "test/mesh_extracted_fracture.exo", comm, gm);

    // extract the fractures
    MeshFactory fac(comm, gm);
    fac.set_preference({ frm });
    auto parent_mesh_cache = Teuchos::rcp(
      new Mesh(parent_mesh, Teuchos::rcp(new AmanziMesh::MeshAlgorithms()), Teuchos::null));
    cacheAll(*parent_mesh_cache);
    auto mesh =
      fac.create(parent_mesh_cache, { "fractures-two" }, AmanziMesh::Entity_kind::FACE, false);
    cacheAll(*mesh);

    // test the surface mesh as a fracture mesh
    // -- mesh audit
    testMeshAudit<MeshAudit, Mesh>(mesh);

    // check we found the right number of things
    int ncells =
      mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
    CHECK_CLOSE_SUMALL(108, ncells, *comm);
    int nnodes =
      mesh->getNumEntities(AmanziMesh::Entity_kind::NODE, AmanziMesh::Parallel_kind::OWNED);
    CHECK_CLOSE_SUMALL(91, nnodes, *comm);
    int nfaces =
      mesh->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED);
    CHECK_CLOSE_SUMALL(198, nfaces, *comm);
    auto ents = mesh->getSetEntities(
      "fractures-two", AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
    CHECK_EQUAL(108, ents.size());
  }
}

TEST(MESH_FRACTURE_EXTRACTION_EXO_MANIFOLD)
{
  // This test is the same as the above, but as it works with the
  // MeshExtractedManifold mesh instead of the framework's extraction
  // constructor, it can be used with any mesh framework that supports edges.

  // create the comm and gm
  auto comm = Amanzi::getDefaultComm();
  std::string infilename = "test/fracture.xml";
  Teuchos::ParameterXMLFileReader xmlreader(infilename);
  Teuchos::ParameterList reg_spec(xmlreader.getParameters());
  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, reg_spec, *comm));

  // a 3D, generated, structured hex on the unit cube, NX=NY=NZ=10
  // works in MSTK
  //
  // extract two inner surfaces as fractures
  std::vector<Framework> frameworks;
  if (framework_enabled(Framework::MSTK)) { frameworks.push_back(Framework::MSTK); }

  // NOTE, moab doesn't support edges, so this can't work here. See #596
  // if (framework_enabled(Framework::MOAB)) {
  //   frameworks.push_back(Framework::MOAB);
  // }

  for (const auto& frm : frameworks) {
    std::cout << std::endl
              << "Testing Fracture Extraction with " << AmanziMesh::to_string(frm) << std::endl
              << "------------------------------------------------" << std::endl;
    auto fac_list = Teuchos::rcp(new Teuchos::ParameterList());
    fac_list->set("request edges", true);
    fac_list->set("request faces", true);
    auto parent_mesh = createFrameworkUnstructured(
      Preference{ frm }, "test/mesh_extracted_fracture.exo", comm, gm, fac_list);

    // extract the fractures
    auto fac_plist = Teuchos::rcp(new Teuchos::ParameterList());
    std::vector<std::string> setnames{ "fractures" };
    fac_plist->sublist("unstructured")
      .sublist("submesh")
      .set<std::string>("extraction method", "manifold mesh");
    MeshFactory fac(comm, gm, fac_plist);
    fac.set_preference({ frm });
    auto parent_mesh_cache = Teuchos::rcp(
      new Mesh(parent_mesh, Teuchos::rcp(new AmanziMesh::MeshAlgorithms()), Teuchos::null));
    cacheAll(*parent_mesh_cache);

    auto mesh = fac.create(parent_mesh_cache, setnames, AmanziMesh::Entity_kind::FACE, false);
    cacheAll(*mesh);

    // test the surface mesh as a fracture mesh
    // -- mesh audit
    testMeshAudit<MeshAudit, Mesh>(mesh);

    // check we found the right number of things
    int ncells =
      mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
    CHECK_CLOSE_SUMALL(108, ncells, *comm);
    int nnodes =
      mesh->getNumEntities(AmanziMesh::Entity_kind::NODE, AmanziMesh::Parallel_kind::OWNED);
    CHECK_CLOSE_SUMALL(91, nnodes, *comm);
    int nfaces =
      mesh->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED);
    CHECK_CLOSE_SUMALL(198, nfaces, *comm);
    auto ents = mesh->getSetEntities(
      "fractures", AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
    CHECK_EQUAL(108, ents.size());
  }
}
