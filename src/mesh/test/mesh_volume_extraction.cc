/*
  Copyright 2010-201x held jointly by LANL, ORNL, LBNL, and PNNL.
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
#include "Mesh.hh"
#include "MeshFactory.hh"

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
  std::vector<double> low{0.0, 0.0, 0.0};
  std::vector<double> high{1.0, 1.0, 1.0};
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
              << "Testing 3D Box 3x3x3 with " << AmanziMesh::framework_names.at(frm) << std::endl
              << "------------------------------------------------" << std::endl;
    auto parent_mesh = createFrameworkStructuredUnitHex(Preference{frm}, 6,6,3, comm, gm, Teuchos::null, false, 2,2,1);

    // extract the volume
    MeshFactory fac(comm, gm);
    fac.set_preference({frm});
    auto mesh = fac.create(parent_mesh, {"Unit Hex"}, AmanziMesh::Entity_kind::CELL);

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
    int ncells = mesh->num_entities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::ALL);
    for (int i = 0; i < ncells; ++i) {
      int parent_cell = mesh->entity_get_parent(AmanziMesh::CELL,i);
      auto cc = mesh->cell_centroid(i);
      auto pcc = parent_mesh->cell_centroid(parent_cell);
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
              << "Testing 3D Box 3x3x3 with " << AmanziMesh::framework_names.at(frm) << std::endl
              << "------------------------------------------------" << std::endl;
    auto parent_mesh = createFrameworkUnstructured(Preference{frm}, filename, comm, gm, Teuchos::null);

    // make sure we can get sets on the mesh
    AmanziMesh::Entity_ID_List set_ids;
    parent_mesh->get_set_entities("Region 1", AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL, &set_ids);
    CHECK_EQUAL(9, set_ids.size());
    parent_mesh->build_columns();

    int ncells = 3;
    AmanziMesh::Entity_ID_List const& cell_list = parent_mesh->cells_of_column(0);
    CHECK_EQUAL(ncells,cell_list.size());

    AmanziMesh::Entity_ID_List const& face_list = parent_mesh->faces_of_column(0);
    CHECK_EQUAL(ncells+1,face_list.size());

    // construct a column mesh by extracting from mesh
    // extract the surface
    MeshFactory fac(comm, gm);
    fac.set_preference({frm});
    auto column_mesh = fac.create(parent_mesh, cell_list, AmanziMesh::Entity_kind::CELL, false);

    // Number of cells in column mesh
    int ncells_col = column_mesh->num_entities(AmanziMesh::CELL,
            AmanziMesh::Parallel_type::OWNED);
    CHECK_EQUAL(ncells,ncells_col);

    // Number of faces in the column mesh
    int nfaces_col = column_mesh->num_entities(AmanziMesh::FACE,
            AmanziMesh::Parallel_type::OWNED);
    CHECK_EQUAL(5*ncells + 1, nfaces_col);

    // Check that their parents are as expected
    for (int i = 0; i < ncells_col; ++i) {
      int parent_cell = column_mesh->entity_get_parent(AmanziMesh::CELL,i);
      CHECK_EQUAL(cell_list[i], parent_cell);
    }

    // check we can still get sets
    AmanziMesh::Entity_ID_List set_ids2;
    bool is_valid = column_mesh->valid_set_name("Region 1", AmanziMesh::CELL);
    CHECK(is_valid);
    column_mesh->get_set_entities("Region 1", AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL, &set_ids2);
    CHECK_EQUAL(1, set_ids2.size());

    set_ids2.clear();
    is_valid = column_mesh->valid_set_name("Top Surface", AmanziMesh::FACE);
    CHECK(is_valid);
    column_mesh->get_set_entities("Top Surface", AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL, &set_ids2);
    CHECK_EQUAL(1, set_ids2.size());

    set_ids2.clear();
    is_valid = column_mesh->valid_set_name("Side Surface", AmanziMesh::FACE);
    CHECK(is_valid);
    column_mesh->get_set_entities("Side Surface", AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL, &set_ids2);
    CHECK_EQUAL(3, set_ids2.size());
  }
}

