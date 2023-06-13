/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

/*
  This test extracts volume cells -- a set of cells in a 3D
  mesh -- and creates a new mesh on this (subset) of cells.
*/

#include <UnitTest++.h>
#include <fstream>

#include "Epetra_Map.h"
#include "AmanziComm.hh"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_Array.hpp"

#include "RegionFactory.hh"
#include "MeshFramework.hh"
#include "MeshFrameworkFactory.hh"

#include "framework_meshes.hh"
#include "geometry_harnesses.hh"
#include "set_harnesses.hh"

using namespace Amanzi;
using namespace Amanzi::AmanziMesh;

// Extract a column from a mesh
TEST(MESH_VOLUME_EXTRACTION_GENERATED)
{
  auto comm = getDefaultComm();
  std::string infilename = "test/hex_3x3x3.xml";
  Teuchos::ParameterXMLFileReader xmlreader(infilename);
  Teuchos::ParameterList reg_spec(xmlreader.getParameters());
  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, reg_spec, *comm));

  // add a region to extract from that is 3D
  Teuchos::ParameterList spec;
  auto& box_reg_spec = spec.sublist("region: box");
  Double_List low{0.0, 0.0, 0.0};
  Double_List high{1.0, 1.0, 1.0};
  box_reg_spec.set<Teuchos::Array<double>>("low coordinate", low);
  box_reg_spec.set<Teuchos::Array<double>>("high coordinate", high);
  gm->AddRegion(AmanziGeometry::createRegion("Unit Hex", gm->size(), spec, *comm));

  // only MSTK currently supports extraction
  std::vector<Framework> frameworks;
  if (framework_enabled(Framework::MSTK)) {
    frameworks.push_back(Framework::MSTK);
  }

  for (const auto& frm : frameworks) {
    std::cout << std::endl
              << "Testing 3D Box 3x3x3 with " << AmanziMesh::to_string(frm) << std::endl
              << "------------------------------------------------" << std::endl;
    auto parent_mesh = createStructuredUnitHex(Preference{frm}, 6,6,3, comm, gm, Teuchos::null, 2,2,1);

    // extract the volume
    MeshFrameworkFactory fac(comm, gm);
    fac.set_preference({frm});

    auto unit_hex_cells = parent_mesh->getSetEntities("Unit Hex",
            AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
    auto vol_framework_mesh = fac.create(parent_mesh,
            unit_hex_cells, AmanziMesh::Entity_kind::CELL);

    // make a MeshCache
    auto mesh = Teuchos::rcp(new Mesh(vol_framework_mesh, Teuchos::rcp(new AmanziMesh::MeshFrameworkAlgorithms()), Teuchos::null));
    mesh->setParentMesh(parent_mesh);

    // test the surface mesh as a 3x3 quad mesh
    // -- mesh audit
    testMeshAudit<MeshAudit, Mesh>(mesh);
    // -- geometry
    testGeometryCube(mesh, 3,3,3);
    // -- exterior maps
    testExteriorMapsUnitBox(mesh,3,3,3);
    // -- sets, which should inherit from the parent mesh
    testHexMeshSets3x3x3(mesh, false, frm);

    // Check that their parents are as expected
    int ncells = mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::ALL);
    for (int i = 0; i < ncells; ++i) {
      int parent_cell = mesh->getEntityParent(AmanziMesh::Entity_kind::CELL,i);
      auto cc = mesh->getCellCentroid(i);
      auto pcc = parent_mesh->getCellCentroid(parent_cell);
      CHECK_CLOSE(cc[0], pcc[0], 1.e-10);
      CHECK_CLOSE(cc[1], pcc[1], 1.e-10);
      CHECK_CLOSE(cc[2], pcc[2], 1.e-10);
    }
  }
}


// Do the same but for an exodus mesh, with sets.
TEST(MESH_VOLUME_EXTRACTION_EXO)
{
  std::string filename("test/hex_3x3x3_sets.exo");
  auto comm = getDefaultComm();
  Teuchos::ParameterList parameterlist;

  // create a sublist name Regions and put a reference to it in
  // reg_spec and other sublists as references.
  Teuchos::ParameterList& reg_spec = parameterlist.sublist("regions");

  Teuchos::ParameterList& top_surface = reg_spec.sublist("Top Surface");
  Teuchos::ParameterList& top_surface_def = top_surface.sublist("region: labeled set");
  top_surface_def.set<std::string>("label","106");
  top_surface_def.set<std::string>("file",filename.c_str());
  top_surface_def.set<std::string>("format","Exodus II");
  top_surface_def.set<std::string>("entity","face");

  Teuchos::ParameterList& side_surface = reg_spec.sublist("Side Surface");
  Teuchos::ParameterList& side_surface_def = side_surface.sublist("region: labeled set");
  side_surface_def.set<std::string>("label","102");
  side_surface_def.set<std::string>("file",filename.c_str());
  side_surface_def.set<std::string>("format","Exodus II");
  side_surface_def.set<std::string>("entity","face");

  Teuchos::ParameterList& r1_surface = reg_spec.sublist("Region 1");
  Teuchos::ParameterList& r1_surface_def = r1_surface.sublist("region: labeled set");
  r1_surface_def.set<std::string>("label","30000");
  r1_surface_def.set<std::string>("file",filename.c_str());
  r1_surface_def.set<std::string>("format","Exodus II");
  r1_surface_def.set<std::string>("entity","cell");

  Teuchos::ParameterList& r2_surface = reg_spec.sublist("Region 2");
  Teuchos::ParameterList& r2_surface_def = r2_surface.sublist("region: labeled set");
  r2_surface_def.set<std::string>("label","20000");
  r2_surface_def.set<std::string>("file",filename.c_str());
  r2_surface_def.set<std::string>("format","Exodus II");
  r2_surface_def.set<std::string>("entity","cell");

  auto gm = Teuchos::rcp(new AmanziGeometry::GeometricModel(3, reg_spec, *comm));

  // only MSTK currently supports extraction
  std::vector<Framework> frameworks;
  if (framework_enabled(Framework::MSTK)) {
    frameworks.push_back(Framework::MSTK);
  }

  for (const auto& frm : frameworks) {
    std::cout << std::endl
              << "Testing 3D Box 3x3x3 with " << AmanziMesh::to_string(frm) << std::endl
              << "------------------------------------------------" << std::endl;
    auto parent_mesh = createUnstructured(Preference{frm}, filename, comm, gm, Teuchos::null);

    // make sure we can get sets on the mesh
    AmanziMesh::Entity_ID_View set_ids = parent_mesh->getSetEntities("Region 1",
            AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::ALL);
    CHECK_EQUAL(9, set_ids.size());
    parent_mesh->buildColumns();

    int ncells = 3;
    AmanziMesh::Entity_ID_View const& cell_list = parent_mesh->columns.cells_.getRowUnmanaged<MemSpace_kind::HOST>(0);
    CHECK_EQUAL(ncells,cell_list.size());

    AmanziMesh::Entity_ID_View const& face_list = parent_mesh->columns.faces_.getRowUnmanaged<MemSpace_kind::HOST>(0);
    CHECK_EQUAL(ncells+1,face_list.size());

    // construct a column mesh by extracting from mesh
    // extract the surface
    MeshFrameworkFactory fac(comm, gm);
    fac.set_preference({frm});
    auto column_mesh_fw = fac.create(parent_mesh, cell_list, AmanziMesh::Entity_kind::CELL);
    auto column_mesh = Teuchos::rcp(new Mesh(column_mesh_fw, Teuchos::rcp(new AmanziMesh::MeshFrameworkAlgorithms()), Teuchos::null));

    // Number of cells in column mesh
    int ncells_col = column_mesh->getNumEntities(AmanziMesh::Entity_kind::CELL,
            AmanziMesh::Parallel_kind::OWNED);
    CHECK_EQUAL(ncells,ncells_col);

    // Number of faces in the column mesh
    int nfaces_col = column_mesh->getNumEntities(AmanziMesh::Entity_kind::FACE,
            AmanziMesh::Parallel_kind::OWNED);
    CHECK_EQUAL(5*ncells + 1, nfaces_col);

    // Check that their parents are as expected
    for (int i = 0; i < ncells_col; ++i) {
      int parent_cell = column_mesh->getEntityParent(AmanziMesh::Entity_kind::CELL,i);
      CHECK_EQUAL(cell_list[i], parent_cell);
    }

    // check we can still get sets
    AmanziMesh::Entity_ID_View set_ids2;
    bool is_valid = column_mesh->isValidSetName("Region 1", AmanziMesh::Entity_kind::CELL);
    CHECK(is_valid);
    set_ids2 = column_mesh->getSetEntities("Region 1", AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::ALL);
    CHECK_EQUAL(1, set_ids2.size());

    is_valid = column_mesh->isValidSetName("Top Surface", AmanziMesh::Entity_kind::FACE);
    CHECK(is_valid);
    set_ids2 = column_mesh->getSetEntities("Top Surface", AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::ALL);
    CHECK_EQUAL(1, set_ids2.size());

    is_valid = column_mesh->isValidSetName("Side Surface", AmanziMesh::Entity_kind::FACE);
    CHECK(is_valid);
    set_ids2 = column_mesh->getSetEntities("Side Surface", AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::ALL);
    CHECK_EQUAL(3, set_ids2.size());
  }
}

