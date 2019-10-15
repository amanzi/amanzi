/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

//!

#include <UnitTest++.h>

#include <mpi.h>
#include <iostream>

#include "Teuchos_RCP.hpp"
#include "AmanziComm.hh"

#include "RegionEnumerated.hh"
#include "MeshLogical.hh"
#include "MeshLogicalFactory.hh"
#include "Geometry.hh"
#include "MeshLogicalAudit.hh"

#include "demo_mesh.hh"

void
test_segment_regular(const Teuchos::RCP<const Amanzi::AmanziMesh::Mesh>& m,
                     bool test_region)
{
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;

  MeshLogicalAudit audit(m, std::cout);
  CHECK(!audit.Verify());

  CHECK_EQUAL(4, m->num_entities(CELL, Parallel_type::ALL));
  CHECK_EQUAL(5, m->num_entities(FACE, Parallel_type::ALL));
  CHECK_EQUAL(0.25, m->cell_volume(0, false));
  CHECK_EQUAL(1.0, m->face_area(0));

  CHECK_EQUAL(1.0, m->face_normal(3)[0]);
  CHECK_EQUAL(0.0, m->face_normal(3)[1]);
  CHECK_EQUAL(0.0, m->face_normal(3)[2]);

  Entity_ID_List faces;
  std::vector<int> dirs;
  std::vector<Point> bisectors;
  m->cell_get_faces_and_dirs(2, &faces, dirs);
  CHECK_EQUAL(2, faces.size());
  CHECK_EQUAL(2, faces[0]);
  CHECK_EQUAL(3, faces[1]);
  CHECK_EQUAL(2, dirs.size());
  CHECK_EQUAL(-1, dirs[0]);
  CHECK_EQUAL(1, dirs[1]);

  faces.clear();
  m->cell_get_faces_and_bisectors(2, &faces, &bisectors);
  CHECK_EQUAL(2, faces.size());
  CHECK_EQUAL(2, faces[0]);
  CHECK_EQUAL(3, faces[1]);
  CHECK_EQUAL(2, bisectors.size());
  CHECK_EQUAL(-0.125, bisectors[0][0]);
  CHECK_EQUAL(0., bisectors[0][1]);
  CHECK_EQUAL(0., bisectors[0][2]);
  CHECK_EQUAL(0.125, bisectors[1][0]);
  CHECK_EQUAL(0., bisectors[1][1]);
  CHECK_EQUAL(0., bisectors[1][2]);

  Entity_ID_List cells;
  m->face_get_cells(0, Parallel_type::ALL, &cells);
  CHECK_EQUAL(1, cells.size());
  CHECK_EQUAL(0, cells[0]);

  m->face_get_cells(1, Parallel_type::ALL, &cells);
  CHECK_EQUAL(2, cells.size());
  CHECK_EQUAL(0, cells[0]);
  CHECK_EQUAL(1, cells[1]);

  m->face_get_cells(4, Parallel_type::ALL, &cells);
  CHECK_EQUAL(1, cells.size());
  CHECK_EQUAL(3, cells[0]);


  if (test_region) {
    // check regions
    CHECK_EQUAL(4, m->get_set_size("myregion", CELL, Parallel_type::ALL));
    CHECK_EQUAL(0, m->get_set_size("myregion", FACE, Parallel_type::ALL));

    Entity_ID_List set_ents;
    m->get_set_entities("myregion", CELL, Parallel_type::ALL, &set_ents);
    CHECK_EQUAL(0, set_ents[0]);
    CHECK_EQUAL(2, set_ents[2]);
  }
}


void
test_segment_irregular(const Teuchos::RCP<Amanzi::AmanziMesh::Mesh>& m,
                       bool test_region)
{
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;

  MeshLogicalAudit audit(m, std::cout);
  CHECK(!audit.Verify());

  Teuchos::RCP<const GeometricModel> gm_c = m->geometric_model();
  Teuchos::RCP<GeometricModel> gm =
    Teuchos::rcp_const_cast<GeometricModel>(gm_c);

  Entity_ID_List ents;
  ents.push_back(0);
  ents.push_back(3);

  Teuchos::RCP<Region> enum_rgn =
    Teuchos::rcp(new RegionEnumerated("myregion", 0, "CELL", ents));
  gm->AddRegion(enum_rgn);


  CHECK_EQUAL(2, m->get_set_size("myregion", CELL, Parallel_type::ALL));
  CHECK_EQUAL(0, m->get_set_size("myregion", FACE, Parallel_type::ALL));

  if (test_region) {
    Entity_ID_List set_ents;
    m->get_set_entities("myregion", CELL, Parallel_type::ALL, &set_ents);
    CHECK_EQUAL(0, set_ents[0]);
    CHECK_EQUAL(3, set_ents[1]);
  }
}


void
test_Y(const Teuchos::RCP<Amanzi::AmanziMesh::Mesh>& m, bool test_region)
{
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;

  MeshLogicalAudit audit(m, std::cout);
  CHECK(!audit.Verify());

  CHECK_EQUAL(11, m->num_entities(CELL, Parallel_type::ALL));
  CHECK_EQUAL(15, m->num_entities(FACE, Parallel_type::ALL));

  // surface
  Point zero(0., 0., 0.);
  CHECK_CLOSE(0., norm(zero - m->face_centroid(0)), 1.e-6);

  // branch point
  Point branch(0., 0., -2.5);
  CHECK_CLOSE(0., norm(branch - m->cell_centroid(2)), 1.e-6);
  branch[2] = -3.0;

  Entity_ID_List branch_faces;
  std::vector<int> dirs;
  m->cell_get_faces_and_dirs(2, &branch_faces, dirs);
  CHECK_EQUAL(5, branch_faces.size());


  CHECK_CLOSE(1.e-4, m->cell_volume(0, false), 1.e-8);
  CHECK_CLOSE(1.e-4, m->cell_volume(2, false), 1.e-8);
  CHECK_CLOSE(0.75 * 0.25e-4, m->cell_volume(3, false), 1.e-8);

  CHECK_CLOSE(1.e-4, m->face_area(2), 1.e-8);
  CHECK_CLOSE(.25e-4, m->face_area(3), 1.e-8);

  if (test_region) {
    CHECK_EQUAL(3, m->get_set_size("coarse_root", CELL, Parallel_type::ALL));
    CHECK_EQUAL(8, m->get_set_size("fine_root", CELL, Parallel_type::ALL));
  }
}


// Tests the construction process, ensures it does not crash.
TEST(MESH_LOGICAL_CONSTRUCTION)
{
  using namespace Amanzi::AmanziMesh;
  Teuchos::RCP<Amanzi::AmanziMesh::MeshLogical> mesh =
    Amanzi::Testing::demoMeshLogicalSegmentRegularManual();
}

// Evaulates the manually constructed mesh.
TEST(MESH_LOGICAL_SEGMENT_REGULAR_MANUAL)
{
  test_segment_regular(Amanzi::Testing::demoMeshLogicalSegmentRegularManual(),
                       false);
}


// Evaluates an irregularly space mesh
TEST(MESH_LOGICAL_SEGMENT_IRREGULAR_WITH_SETS)
{
  test_segment_irregular(
    Amanzi::Testing::demoMeshLogicalSegmentIrregularManual(), true);
}


// updates to a Y-mesh
TEST(MESH_LOGICAL_Y_DEFORMED)
{
  using namespace Amanzi::AmanziMesh;
  auto mesh = Amanzi::Testing::demoMeshLogicalYManual();

  MeshLogicalAudit audit(mesh, std::cout);
  CHECK(!audit.Verify());


  std::vector<double> cv, fa;
  std::vector<std::vector<double>> cf_lens;
  std::vector<Amanzi::AmanziGeometry::Point> centroids;
  mesh->get_logical_geometry(&cv, &cf_lens, &fa, &centroids);

  CHECK_CLOSE(1.e-4, fa[0], 1.e-10);
  CHECK_CLOSE(2.0 * 0.25 * 1.e-4, cv[0], 1.e-10);
  CHECK_CLOSE(0.25, cf_lens[0][0], 1.e-10);

  Amanzi::AmanziGeometry::Point zero(0., 0., -0.25);

  Amanzi::AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;
  mesh->cell_get_faces_and_dirs(0, &faces, dirs);
  CHECK_EQUAL(2, faces.size());
  CHECK_EQUAL(2, cf_lens[0].size());


  CHECK_CLOSE(0., Amanzi::AmanziGeometry::norm(zero - centroids[0]), 1.e-6);


  // grow the first cell
  cv[0] *= 2.0;
  fa[0] *= 1.5;

  cf_lens[0][0] *= 2.0;
  cf_lens[0][1] *= 2.0;

  mesh->set_logical_geometry(&cv, &cf_lens, &fa, NULL);

  CHECK_CLOSE(1.5e-4, mesh->face_area(0), 1.e-10);
  CHECK_CLOSE(2.0 * 2.0 * 0.25 * 1.e-4, mesh->cell_volume(0, false), 1.e-10);
}


// Evaluates a Y-mesh
TEST(MESH_LOGICAL_Y_XML_WITH_SETS)
{
  test_Y(Amanzi::Testing::demoMeshLogicalYFromXML(), true);
}


// Evaluates a Y-mesh embedded in another mesh
TEST(MESH_EMBEDDED_Y)
{
  auto m = Amanzi::Testing::demoMeshLogicalYEmbedded();
  Amanzi::AmanziMesh::MeshLogicalAudit audit(m, std::cout);
  CHECK(!audit.Verify());
}

// subgrid model
TEST(MESH_SUBGRID_VARIABLE_TAU)
{
  auto m = Amanzi::Testing::demoMeshLogicalYFromXML("subgrid mesh");
  Amanzi::AmanziMesh::MeshLogicalAudit audit(m, std::cout);
  CHECK(!audit.Verify());
}
