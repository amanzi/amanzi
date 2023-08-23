/*
  Copyright 2010-201x held jointly by LANL, ORNL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

// Test: column and surface cell extracted meshes
//   These are 1-off, simplified meshes that make 1D and 0D
//   conceptual meshes for use in subgrid models.  This extracts
//   the column, creates a column mesh, and creates a surface cell
//   mesh, all using the factory, and then checks sets on these.

#include "UnitTest++.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"

#include "errors.hh"
#include "AmanziComm.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"

using namespace Amanzi;

TEST(SURFACE_COLUMN_MESH_3D_UNSTRUCTURED_SETS)
{
  // This test uses the proper nice clean constructors, and looks like what
  // user code will likely do.  It then tests sets and makes sure that set
  // inheritance works from volume mesh to column mesh to surface cell mesh.
  std::string mesh_filename("test/hex_3x3x3_sets.exo");
  auto comm = getDefaultComm();

  std::string xml_filename = "test/hex_3x3x3.xml";
  Teuchos::ParameterXMLFileReader xmlreader(xml_filename);
  Teuchos::ParameterList reg_spec(xmlreader.getParameters());
  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, reg_spec, *comm));

  // make sure we partition with zoltan rcb to get columns
  auto mesh_plist = Teuchos::rcp(new Teuchos::ParameterList("mesh"));
  mesh_plist->set("partitioner", "zoltan_rcb");

  // a 3D exodus file, structured hex on the unit cube, NX=NY=NZ=3
  std::vector<AmanziMesh::Framework> frameworks;
  if (framework_enabled(AmanziMesh::Framework::MSTK)) {
    frameworks.push_back(AmanziMesh::Framework::MSTK);
  }

  for (const auto& frm : frameworks) {
    std::cout << std::endl
              << "Extracting columns and surface cell using " << AmanziMesh::to_string(frm) << std::endl
              << "------------------------------------------------" << std::endl;
    AmanziMesh::MeshFactory fac(comm, gm, mesh_plist);
    auto mesh = fac.create(mesh_filename);
    CHECK_EQUAL(9, mesh->getSetSize("Top Face Plane", AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED));
    CHECK_EQUAL(9, mesh->getSetSize("Top Box", AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED));
    CHECK_EQUAL(9, mesh->getSetSize("Face 106", AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED));
    CHECK_EQUAL(9, mesh->getSetSize("Face 103", AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED));
    CHECK_EQUAL(9, mesh->getSetSize("Top LS", AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED));

    // Create a column mesh from one of the columns
    auto mesh_col = fac.createColumn(mesh, 0, mesh_plist);

    // Create a surface cell from that column
    auto mesh_sc = fac.createSurfaceCell(mesh_col);

    // check geometry of the column
    CHECK_EQUAL(3, mesh_col->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED));
    CHECK_EQUAL(4, mesh_col->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED));
    CHECK_EQUAL(16, mesh_col->getNumEntities(AmanziMesh::Entity_kind::NODE, AmanziMesh::Parallel_kind::OWNED));

    // check regions of the column -- geometric regions
    CHECK_EQUAL(1, mesh_col->getSetSize("Top Face Plane", AmanziMesh::Entity_kind::FACE,
            AmanziMesh::Parallel_kind::ALL));
    CHECK_EQUAL(1, mesh_col->getSetSize("Top Box", AmanziMesh::Entity_kind::CELL,
            AmanziMesh::Parallel_kind::ALL));
    // -- this region is the side, and columns have no side faces!
    CHECK_THROW(mesh_col->getSetSize("West Face Plane", AmanziMesh::Entity_kind::FACE,
            AmanziMesh::Parallel_kind::ALL), Errors::Message);

    // check regions of the column -- labeled sets
    CHECK_EQUAL(1, mesh_col->getSetSize("Face 106", AmanziMesh::Entity_kind::FACE,
            AmanziMesh::Parallel_kind::ALL));
    CHECK_EQUAL(1, mesh_col->getSetSize("Top LS", AmanziMesh::Entity_kind::CELL,
            AmanziMesh::Parallel_kind::ALL));
    // -- this region is the side, and columns have no side faces!
    CHECK_THROW(mesh_col->getSetSize("Face 103", AmanziMesh::Entity_kind::FACE,
            AmanziMesh::Parallel_kind::ALL), Errors::Message);

    // check geometry of the surface cell
    CHECK_EQUAL(1, mesh_sc->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED));
    CHECK_EQUAL(4, mesh_sc->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED));
    CHECK_EQUAL(4, mesh_sc->getNumEntities(AmanziMesh::Entity_kind::NODE, AmanziMesh::Parallel_kind::OWNED));

    // geometric regions
    CHECK_EQUAL(1, mesh_sc->getSetSize("Top Face Plane", AmanziMesh::Entity_kind::CELL,
            AmanziMesh::Parallel_kind::ALL));
    CHECK_EQUAL(1, mesh_sc->getSetSize("Top Box", AmanziMesh::Entity_kind::CELL,
            AmanziMesh::Parallel_kind::ALL));
    // -- this region is the side, and columns have no side faces!
    CHECK_THROW(mesh_sc->getSetSize("West Face Plane", AmanziMesh::Entity_kind::FACE,
            AmanziMesh::Parallel_kind::ALL), Errors::Message);

    // labeled set regions
    CHECK_EQUAL(1, mesh_sc->getSetSize("Face 106", AmanziMesh::Entity_kind::CELL,
            AmanziMesh::Parallel_kind::ALL));
    CHECK_EQUAL(1, mesh_sc->getSetSize("Top LS", AmanziMesh::Entity_kind::CELL,
            AmanziMesh::Parallel_kind::ALL));
    // -- this region is the side, and columns have no side faces!
    CHECK_THROW(mesh_sc->getSetSize("Face 103", AmanziMesh::Entity_kind::FACE,
            AmanziMesh::Parallel_kind::ALL), Errors::Message);

    // check entities
    AmanziMesh::Entity_ID_View cells_in_surf =
      mesh_sc->getSetEntities("Top Face Plane", AmanziMesh::Entity_kind::CELL,
              AmanziMesh::Parallel_kind::OWNED);
    CHECK_EQUAL(1, cells_in_surf.size());
    CHECK_EQUAL(0, cells_in_surf[0]);

    AmanziMesh::Entity_ID_View cells_in_surf2 =
    mesh_sc->getSetEntities("Face 106", AmanziMesh::Entity_kind::CELL,
            AmanziMesh::Parallel_kind::OWNED);
    CHECK_EQUAL(1, cells_in_surf2.size());
    CHECK_EQUAL(0, cells_in_surf2[0]);
  }
}
