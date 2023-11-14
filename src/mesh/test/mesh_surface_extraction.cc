/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

/*
  This test extracts surface faces -- a set of faces on the boundary of a 3D
  mesh -- and creates a (flattened) 2D mesh from these faces.
*/

#include "UnitTest++.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"

#include "AmanziComm.hh"
#include "RegionFactory.hh"
#include "Mesh.hh"
#include "MeshAudit.hh"
#include "MeshFrameworkFactory.hh"
#include "MeshExtractedManifold.hh"

#include "framework_meshes.hh"
#include "geometry_harnesses.hh"
#include "set_harnesses.hh"

TEST(MESH_SURFACE_EXTRACTION_GENERATED)
{
  // create the comm and gm
  auto comm = Amanzi::getDefaultComm();
  std::string infilename = "test/quad_3x3.xml";
  Teuchos::ParameterXMLFileReader xmlreader(infilename);
  Teuchos::ParameterList reg_spec(xmlreader.getParameters());
  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, reg_spec, *comm));

  // add a region to extract from that is 3D
  Teuchos::ParameterList spec;
  auto& surf_reg_spec = spec.sublist("region: plane");
  Double_List point{ 0.0, 0.0, 1.0 };
  surf_reg_spec.set<Teuchos::Array<double>>("point", point);
  Double_List normal{ 0.0, 0.0, 1.0 };
  surf_reg_spec.set<Teuchos::Array<double>>("normal", normal);
  gm->AddRegion(AmanziGeometry::createRegion("Top Face Plane", gm->size(), spec, *comm));

  // a 3D, generated, structured hex on the unit cube, NX=NY=NZ=3,
  // using the framework extraction, only works in MSTK
  //
  // extract & flatten the top surface to form a 2D mesh, NX=NY=3
  std::vector<Framework> frameworks;
  if (framework_enabled(Framework::MSTK)) { frameworks.push_back(Framework::MSTK); }

  for (const auto& frm : frameworks) {
    std::cout << std::endl
              << "Extracting surface from 3D Generated Box with " << AmanziMesh::to_string(frm)
              << std::endl
              << "------------------------------------------------" << std::endl;

    auto fac_list = Teuchos::rcp(new Teuchos::ParameterList("factory list"));
    fac_list->set<std::string>("partitioner", "zoltan_rcb");
    auto parent_mesh = createStructuredUnitHex(Preference{ frm }, 3, 3, 3, comm, gm, fac_list);

    // extract the surface
    MeshFrameworkFactory fac(comm, gm);
    fac.set_preference({ frm });
    auto top_faces = parent_mesh->getSetEntities(
      "Top Face Plane", AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED);

    auto surface_framework_mesh =
      fac.create(parent_mesh, top_faces, AmanziMesh::Entity_kind::FACE, true);

    // make a MeshCache
    auto mesh = Teuchos::rcp(new Mesh(surface_framework_mesh,
                                      Teuchos::rcp(new AmanziMesh::MeshFrameworkAlgorithms()),
                                      Teuchos::null));
    mesh->setParentMesh(parent_mesh);

    // test the surface mesh as a 3x3 quad mesh
    // -- mesh audit
    testMeshAudit<MeshAudit, Mesh>(mesh);
    // -- geometry
    testGeometryQuad(mesh, 3, 3);
    // -- exterior maps
    testExteriorMapsUnitBox(mesh, 3, 3);
    // -- sets, which should inherit from the parent mesh
    testQuadMeshSets3x3(mesh, false, frm, true);
  }
}


TEST(MESH_SURFACE_EXTRACTION_EXO)
{
  // create the comm and gm
  auto comm = Amanzi::getDefaultComm();
  std::string infilename = "test/quad_3x3.xml";
  Teuchos::ParameterXMLFileReader xmlreader(infilename);
  Teuchos::ParameterList reg_spec(xmlreader.getParameters());
  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, reg_spec, *comm));

  // add a region to extract from that is 3D
  Teuchos::ParameterList spec;
  auto& surf_reg_spec = spec.sublist("region: plane");
  Double_List point{ 0.0, 0.0, 1.0 };
  surf_reg_spec.set<Teuchos::Array<double>>("point", point);
  Double_List normal{ 0.0, 0.0, 1.0 };
  surf_reg_spec.set<Teuchos::Array<double>>("normal", normal);
  gm->AddRegion(AmanziGeometry::createRegion("Top Face Plane", gm->size(), spec, *comm));

  // a 3D, generated, structured hex on the unit cube, NX=NY=NZ=3
  // works in MSTK
  //
  // extract & flatten the top surface to form a 2D mesh, NX=NY=3
  std::vector<Framework> frameworks;
  if (framework_enabled(Framework::MSTK)) { frameworks.push_back(Framework::MSTK); }

  for (const auto& frm : frameworks) {
    std::cout << std::endl
              << "Extracting surface from EXO 3D Box with " << AmanziMesh::to_string(frm)
              << std::endl
              << "------------------------------------------------" << std::endl;
    auto fac_list = Teuchos::rcp(new Teuchos::ParameterList("factory list"));
    fac_list->set<std::string>("partitioner", "zoltan_rcb");
    auto parent_mesh =
      createUnstructured(Preference{ frm }, "test/hex_3x3x3_sets.exo", comm, gm, fac_list);

    // extract the surface
    MeshFrameworkFactory fac(comm, gm);
    fac.set_preference({ frm });
    auto top_faces = parent_mesh->getSetEntities(
      "Top Face Plane", AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED);
    auto surface_framework_mesh =
      fac.create(parent_mesh, top_faces, AmanziMesh::Entity_kind::FACE, true);

    // make a MeshCache
    auto mesh = Teuchos::rcp(new Mesh(surface_framework_mesh,
                                      Teuchos::rcp(new AmanziMesh::MeshFrameworkAlgorithms()),
                                      Teuchos::null));
    mesh->setParentMesh(parent_mesh);

    // test the surface mesh as a 3x3 quad mesh
    // -- mesh audit
    testMeshAudit<MeshAudit, Mesh>(mesh);
    // -- geometry
    testGeometryQuad(mesh, 3, 3);
    // -- exterior maps
    testExteriorMapsUnitBox(mesh, 3, 3);
    // -- sets, which should inherit from the parent mesh
    testQuadMeshSets3x3(mesh, true, frm, true);
  }
}

TEST(MESH_SURFACE_EXTRACTION_GENERATED_EXTRACTED_MANIFOLD)
{
  // Unlike the above tests, which call the framework's extraction constructor,
  // this test creates the internal MeshExtractedManifold object instead.

  // create the comm and gm
  auto comm = Amanzi::getDefaultComm();
  std::string infilename = "test/quad_3x3.xml";
  Teuchos::ParameterXMLFileReader xmlreader(infilename);
  Teuchos::ParameterList reg_spec(xmlreader.getParameters());
  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, reg_spec, *comm));

  // add a region to extract from that is 3D
  Teuchos::ParameterList spec;
  auto& surf_reg_spec = spec.sublist("region: plane");
  Double_List point{ 0.0, 0.0, 1.0 };
  surf_reg_spec.set<Teuchos::Array<double>>("point", point);
  Double_List normal{ 0.0, 0.0, 1.0 };
  surf_reg_spec.set<Teuchos::Array<double>>("normal", normal);
  gm->AddRegion(AmanziGeometry::createRegion("Top Face Plane", gm->size(), spec, *comm));

  // a 3D, generated, structured hex on the unit cube, NX=NY=NZ=3
  // works in MSTK and Simple
  //
  // extract & flatten the top surface to form a 2D mesh, NX=NY=3
  std::vector<Framework> frameworks;
  //if (comm->NumProc() == 1) frameworks.push_back(Framework::SIMPLE);
  if (framework_enabled(Framework::MSTK)) { frameworks.push_back(Framework::MSTK); }

  for (const auto& frm : frameworks) {
    std::cout << std::endl
              << "MeshExtractedManifold from surface of 3D Generated Box with "
              << AmanziMesh::to_string(frm) << std::endl
              << "------------------------------------------------" << std::endl;
    auto fac_list = Teuchos::rcp(new Teuchos::ParameterList("factory list"));
    fac_list->set<std::string>("partitioner", "zoltan_rcb");
    fac_list->set("request edges", true);
    auto parent_mesh =
      createFrameworkStructuredUnitHex(Preference{ frm }, 3, 3, 3, comm, gm, fac_list);

    // extract the surface
    auto fac_plist = Teuchos::rcp(new Teuchos::ParameterList());
    fac_plist->sublist("unstructured")
      .sublist("submesh")
      .set<std::string>("extraction method", "manifold mesh");
    MeshFactory fac(comm, gm, fac_plist);
    fac.set_preference({ frm });

    auto parent_mesh_cache = Teuchos::rcp(new Mesh(
      parent_mesh, Teuchos::rcp(new AmanziMesh::MeshFrameworkAlgorithms()), Teuchos::null));

    std::vector<std::string> setnames({ "Top Face Plane" });
    auto mesh = fac.create(parent_mesh_cache, setnames, AmanziMesh::Entity_kind::FACE, true);

    // test the surface mesh as a 3x3 quad mesh
    // -- mesh audit
    testMeshAudit<MeshAudit, Mesh>(mesh);
    // -- geometry
    testGeometryQuad(mesh, 3, 3);

    // -- exterior maps -- NOT SUPPORTED BY SIMPLE
    if (frm != Framework::SIMPLE) testExteriorMapsUnitBox(mesh, 3, 3);

    // -- sets, which should inherit from the parent mesh
    testQuadMeshSets3x3(mesh, false, frm, true);
  }
}

TEST(MESH_SURFACE_EXTRACTION_EXO_EXTRACTED_MANIFOLD)
{
  // create the comm and gm
  auto comm = Amanzi::getDefaultComm();
  std::string infilename = "test/quad_3x3.xml";
  Teuchos::ParameterXMLFileReader xmlreader(infilename);
  Teuchos::ParameterList reg_spec(xmlreader.getParameters());
  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, reg_spec, *comm));

  // add a region to extract from that is 3D
  Teuchos::ParameterList spec;
  auto& surf_reg_spec = spec.sublist("region: plane");
  Double_List point{ 0.0, 0.0, 1.0 };
  surf_reg_spec.set<Teuchos::Array<double>>("point", point);
  Double_List normal{ 0.0, 0.0, 1.0 };
  surf_reg_spec.set<Teuchos::Array<double>>("normal", normal);
  gm->AddRegion(AmanziGeometry::createRegion("Top Face Plane", gm->size(), spec, *comm));

  // a 3D, generated, structured hex on the unit cube, NX=NY=NZ=3
  // works in MSTK
  //
  // extract & flatten the top surface to form a 2D mesh, NX=NY=3
  std::vector<Framework> frameworks;
  if (framework_enabled(Framework::MSTK)) { frameworks.push_back(Framework::MSTK); }

  // MOAB does not support edges, so cannot work with MeshExtractedManifold.  See #596 part 3
  // if (framework_enabled(Framework::MOAB)) {
  //   frameworks.push_back(Framework::MOAB);
  // }

  for (const auto& frm : frameworks) {
    std::cout << std::endl
              << "MeshExtractedManifold from surface of 3D EXO Box with "
              << AmanziMesh::to_string(frm) << std::endl
              << "------------------------------------------------" << std::endl;
    auto fac_list = Teuchos::rcp(new Teuchos::ParameterList("factory list"));
    fac_list->set<std::string>("partitioner", "zoltan_rcb");
    fac_list->set("request edges", true);
    auto parent_mesh =
      createFrameworkUnstructured(Preference{ frm }, "test/hex_3x3x3_sets.exo", comm, gm, fac_list);

    // extract the surface
    auto fac_plist = Teuchos::rcp(new Teuchos::ParameterList());
    fac_plist->sublist("unstructured")
      .sublist("submesh")
      .set<std::string>("extraction method", "manifold mesh");

    MeshFactory fac(comm, gm, fac_plist);
    fac.set_preference({ frm });
    auto parent_mesh_cache = Teuchos::rcp(new Mesh(
      parent_mesh, Teuchos::rcp(new AmanziMesh::MeshFrameworkAlgorithms()), Teuchos::null));

    std::vector<std::string> setnames({ "Top Face Plane" });
    auto mesh = fac.create(parent_mesh_cache, setnames, AmanziMesh::Entity_kind::FACE, true);

    // test the surface mesh as a 3x3 quad mesh
    // -- mesh audit
    testMeshAudit<MeshAudit, Mesh>(mesh);
    // -- geometry
    testGeometryQuad(mesh, 3, 3);
    // -- exterior maps
    testExteriorMapsUnitBox(mesh, 3, 3);

    // -- sets, which should inherit from the parent mesh
    testQuadMeshSets3x3(mesh, true, frm, true);
  }
}
