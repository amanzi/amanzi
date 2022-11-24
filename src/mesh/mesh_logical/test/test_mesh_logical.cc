/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

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

using namespace Amanzi::AmanziMesh;
using namespace Amanzi::AmanziGeometry;

bool
POINT_CLOSE(const Point& p1, const Point& p2)
{
  CHECK_EQUAL(p1.dim(), p2.dim());
  CHECK_CLOSE(p1[0], p2[0], 1.e-8);
  CHECK_CLOSE(p1[1], p2[1], 1.e-8);
  if (p1.dim() > 2) CHECK_CLOSE(p1[2], p2[2], 1.e-8);
  return norm(p1 - p2) < 1.e-8;
}

#define CHECK_POINT_CLOSE(p1, p2) CHECK(POINT_CLOSE(p1, p2))

void
test_segment_regular(const Teuchos::RCP<const Amanzi::AmanziMesh::Mesh>& m, bool test_region)
{
  MeshLogicalAudit audit(m, std::cout);
  CHECK(!audit.Verify());

  CHECK_EQUAL(4, m->num_entities(CELL, Parallel_type::ALL));
  CHECK_EQUAL(5, m->num_entities(FACE, Parallel_type::ALL));

  for (int i = 0; i != 4; ++i) {
    CHECK_EQUAL(0.25, m->cell_volume(i));

    Entity_ID_List faces;
    std::vector<int> dirs;
    std::vector<Point> bisectors;
    m->cell_get_faces_and_dirs(i, &faces, &dirs);
    CHECK_EQUAL(2, faces.size());
    CHECK_EQUAL(i, faces[0]);
    CHECK_EQUAL(i + 1, faces[1]);
    CHECK_EQUAL(2, dirs.size());
    CHECK_EQUAL(1, dirs[0]);
    CHECK_EQUAL(-1, dirs[1]);

    faces.clear();
    m->cell_get_faces_and_bisectors(i, &faces, &bisectors);
    CHECK_EQUAL(2, faces.size());
    CHECK_EQUAL(i, faces[0]);
    CHECK_EQUAL(i + 1, faces[1]);
    CHECK_EQUAL(2, bisectors.size());
    CHECK_POINT_CLOSE(Point(0.125, 0., 0.), bisectors[0]);
    CHECK_POINT_CLOSE(Point(-0.125, 0., 0.), bisectors[1]);
  }

  for (int i = 0; i != 5; ++i) {
    CHECK_EQUAL(1.0, m->face_area(0));
    auto normal = m->face_normal(i);
    CHECK_POINT_CLOSE(Point(1., 0., 0.), normal);

    Entity_ID_List cells;
    m->face_get_cells(i, Parallel_type::ALL, &cells);
    if (i == 0) {
      CHECK_EQUAL(1, cells.size());
      CHECK_EQUAL(0, cells[0]);
    } else if (i == 4) {
      CHECK_EQUAL(1, cells.size());
      CHECK_EQUAL(3, cells[0]);
    } else {
      CHECK_EQUAL(2, cells.size());
      CHECK_EQUAL(i - 1, cells[0]);
      CHECK_EQUAL(i, cells[1]);
    }
  }


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
test_segment_irregular(const Teuchos::RCP<Amanzi::AmanziMesh::Mesh>& m, bool test_region)
{
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;

  MeshLogicalAudit audit(m, std::cout);
  CHECK(!audit.Verify());

  Teuchos::RCP<const GeometricModel> gm_c = m->geometric_model();
  Teuchos::RCP<GeometricModel> gm = Teuchos::rcp_const_cast<GeometricModel>(gm_c);

  Entity_ID_List ents;
  ents.push_back(0);
  ents.push_back(3);

  Teuchos::RCP<Region> enum_rgn = Teuchos::rcp(new RegionEnumerated("myregion", 0, "CELL", ents));
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
  std::cout << "branch = " << branch << std::endl;
  std::cout << "centroid = " << m->cell_centroid(2) << std::endl;
  CHECK_CLOSE(0., norm(branch - m->cell_centroid(2)), 1.e-6);
  branch[2] = -3.0;

  Entity_ID_List branch_faces;
  std::vector<int> dirs;
  m->cell_get_faces_and_dirs(2, &branch_faces, &dirs);
  CHECK_EQUAL(5, branch_faces.size());

  for (int i = 0; i != 3; ++i) CHECK_CLOSE(1.e-4, m->cell_volume(i), 1.e-8);
  for (int i = 3; i != 11; ++i) CHECK_CLOSE(.75 * 0.25e-4, m->cell_volume(i), 1.e-8);

  for (int i = 0; i != 3; ++i) CHECK_CLOSE(1.e-4, m->face_area(i), 1.e-8);
  for (int i = 3; i != 15; ++i) CHECK_CLOSE(.25e-4, m->face_area(i), 1.e-8);

  if (test_region) {
    CHECK_EQUAL(3, m->get_set_size("coarse_root", CELL, Parallel_type::ALL));
    CHECK_EQUAL(8, m->get_set_size("fine_root", CELL, Parallel_type::ALL));
  }
}


void
test_2Y(const Teuchos::RCP<Amanzi::AmanziMesh::Mesh>& m, bool test_region)
{
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;

  MeshLogicalAudit audit(m, std::cout);
  CHECK(!audit.Verify());

  CHECK_EQUAL(3, m->num_entities(CELL, Parallel_type::ALL));
  CHECK_EQUAL(5, m->num_entities(FACE, Parallel_type::ALL));

  Entity_ID_List branch_faces;
  std::vector<int> dirs;
  std::vector<Point> bisectors;
  double r22 = sqrt(2.0) / 2.0;

  // check topology/geometry
  std::cout << "  Checking topology/geometry" << std::endl;
  // cell 0
  m->cell_get_faces_and_dirs(0, &branch_faces, &dirs);
  CHECK_EQUAL(3, branch_faces.size());
  CHECK_EQUAL(0, branch_faces[0]);
  CHECK_EQUAL(1, branch_faces[1]);
  CHECK_EQUAL(3, branch_faces[2]);

  CHECK_EQUAL(3, dirs.size());
  CHECK_CLOSE(1, dirs[0], 1.e-8);
  CHECK_CLOSE(-1, dirs[1], 1.e-8);
  CHECK_CLOSE(-1, dirs[2], 1.e-8);

  CHECK_CLOSE(4.0, m->cell_volume(0), 1.e-8);

  m->cell_get_faces_and_bisectors(0, &branch_faces, &bisectors);
  CHECK_EQUAL(3, branch_faces.size());
  CHECK_EQUAL(0, branch_faces[0]);
  CHECK_EQUAL(1, branch_faces[1]);
  CHECK_EQUAL(3, branch_faces[2]);

  CHECK_EQUAL(3, bisectors.size());
  CHECK_POINT_CLOSE(Point(1.0, 0., 0.), bisectors[0]);
  CHECK_POINT_CLOSE(Point(-1.0, 0., 0.), bisectors[1]);
  CHECK_POINT_CLOSE(Point(-1.0, 0., 0.), bisectors[2]);

  // cell 1
  m->cell_get_faces_and_dirs(1, &branch_faces, &dirs);
  CHECK_EQUAL(2, branch_faces.size());
  CHECK_EQUAL(1, branch_faces[0]);
  CHECK_EQUAL(2, branch_faces[1]);

  CHECK_EQUAL(2, dirs.size());
  CHECK_CLOSE(1, dirs[0], 1.e-8);
  CHECK_CLOSE(-1, dirs[1], 1.e-8);

  CHECK_CLOSE(1.0, m->cell_volume(1), 1.e-8);

  m->cell_get_faces_and_bisectors(1, &branch_faces, &bisectors);
  CHECK_EQUAL(2, branch_faces.size());
  CHECK_EQUAL(1, branch_faces[0]);
  CHECK_EQUAL(2, branch_faces[1]);

  CHECK_EQUAL(2, bisectors.size());
  CHECK_POINT_CLOSE(Point(r22 / 2., r22 / 2, 0.), bisectors[0]);
  CHECK_POINT_CLOSE(Point(-r22 / 2., -r22 / 2, 0.), bisectors[1]);

  // cell 2
  m->cell_get_faces_and_dirs(2, &branch_faces, &dirs);
  CHECK_EQUAL(2, branch_faces.size());
  CHECK_EQUAL(3, branch_faces[0]);
  CHECK_EQUAL(4, branch_faces[1]);

  CHECK_EQUAL(2, dirs.size());
  CHECK_CLOSE(1, dirs[0], 1.e-8);
  CHECK_CLOSE(-1, dirs[1], 1.e-8);

  CHECK_CLOSE(1.0, m->cell_volume(2), 1.e-8);

  m->cell_get_faces_and_bisectors(2, &branch_faces, &bisectors);
  CHECK_EQUAL(2, branch_faces.size());
  CHECK_EQUAL(3, branch_faces[0]);
  CHECK_EQUAL(4, branch_faces[1]);

  CHECK_EQUAL(2, bisectors.size());
  CHECK_POINT_CLOSE(Point(r22 / 2., -r22 / 2, 0.), bisectors[0]);
  CHECK_POINT_CLOSE(Point(-r22 / 2., r22 / 2, 0.), bisectors[1]);


  // faces
  Amanzi::AmanziGeometry::Point f0_normal = m->face_normal(0);
  CHECK_POINT_CLOSE(Point(2, 0, 0), f0_normal);
  CHECK_CLOSE(2.0, m->face_area(0), 1.e-8);

  Amanzi::AmanziGeometry::Point f1_normal = m->face_normal(1);
  CHECK_POINT_CLOSE(Point(1, 0, 0), f1_normal);
  CHECK_CLOSE(1.0, m->face_area(1), 1.e-8);

  Amanzi::AmanziGeometry::Point f2_normal = m->face_normal(2);
  CHECK_POINT_CLOSE(Point(sqrt(2) / 2, sqrt(2) / 2, 0), f2_normal);
  CHECK_CLOSE(1.0, m->face_area(2), 1.e-8);

  Amanzi::AmanziGeometry::Point f3_normal = m->face_normal(3);
  CHECK_POINT_CLOSE(Point(1, 0, 0), f3_normal);
  CHECK_CLOSE(1.0, m->face_area(3), 1.e-8);

  Amanzi::AmanziGeometry::Point f4_normal = m->face_normal(4);
  CHECK_POINT_CLOSE(Point(sqrt(2) / 2, -sqrt(2) / 2, 0), f4_normal);
  CHECK_CLOSE(1.0, m->face_area(4), 1.e-8);

  // check centroids
  std::cout << "  Checking centroids" << std::endl;
  CHECK_POINT_CLOSE(Point(2., 0., 0.), m->face_centroid(0));
  CHECK_POINT_CLOSE(Point(1., 0., 0.), m->cell_centroid(0));
  CHECK_POINT_CLOSE(Point(0., 0., 0.), m->face_centroid(1));
  CHECK_POINT_CLOSE(Point(-r22 / 2., -r22 / 2., 0.), m->cell_centroid(1));
  CHECK_POINT_CLOSE(Point(-r22, -r22, 0.), m->face_centroid(2));
  CHECK_POINT_CLOSE(Point(0., 0., 0.), m->face_centroid(3));
  CHECK_POINT_CLOSE(Point(-r22 / 2., r22 / 2., 0.), m->cell_centroid(2));
  CHECK_POINT_CLOSE(Point(-r22, r22, 0.), m->face_centroid(4));

  if (test_region) {
    CHECK_EQUAL(1, m->get_set_size("coarse_root", CELL, Parallel_type::ALL));
    CHECK_EQUAL(2, m->get_set_size("fine_root", CELL, Parallel_type::ALL));
  }
}


// Tests the construction process, ensures it does not crash.
TEST(MESH_LOGICAL_CONSTRUCTION)
{
  using namespace Amanzi::AmanziMesh;
  std::cout << std::endl
            << "TEST: MeshLogical Construction" << std::endl
            << "------------------------------" << std::endl;
  auto m = Amanzi::Testing::demoMeshLogicalSegmentRegularManual();
}

// Evaulates the manually constructed mesh.
TEST(MESH_LOGICAL_SEGMENT_REGULAR_MANUAL)
{
  std::cout << std::endl
            << "TEST: MeshLogical single segment, manual construction" << std::endl
            << "-----------------------------------------------------" << std::endl;
  test_segment_regular(Amanzi::Testing::demoMeshLogicalSegmentRegularManual(), false);
}

// Evaulates the manually constructed mesh.
TEST(MESH_LOGICAL_SEGMENT_REGULAR_XML)
{
  std::cout << std::endl
            << "TEST: MeshLogical single segment, xml construction" << std::endl
            << "-----------------------------------------------------" << std::endl;
  test_segment_regular(Amanzi::Testing::demoMeshLogicalFromXML("regular"), false);
}


// Evaluates an irregularly space mesh
TEST(MESH_LOGICAL_SEGMENT_IRREGULAR_WITH_SETS)
{
  std::cout << std::endl
            << "TEST: MeshLogical single segment, deformed" << std::endl
            << "-----------------------------------------------------" << std::endl;
  test_segment_irregular(Amanzi::Testing::demoMeshLogicalSegmentIrregularManual(), true);
}

// Evaluates a 2Y-mesh
TEST(MESH_LOGICAL_2Y_XML_WITH_SETS)
{
  std::cout << std::endl
            << "TEST: MeshLogical 2Y, from XML" << std::endl
            << "-----------------------------------------------------" << std::endl;
  test_2Y(Amanzi::Testing::demoMeshLogicalFromXML("logical mesh 2Y"), true);
}


// Evaluates a Y-mesh
TEST(MESH_LOGICAL_Y)
{
  std::cout << std::endl
            << "TEST: MeshLogical Y, manual construction" << std::endl
            << "-----------------------------------------------------" << std::endl;
  auto m = Amanzi::Testing::demoMeshLogicalYManual();
  test_Y(m, true);
}

// Evaluates a Y-mesh
TEST(MESH_LOGICAL_Y_XML_WITH_SETS)
{
  std::cout << std::endl
            << "TEST: MeshLogical Y, from XML" << std::endl
            << "-----------------------------------------------------" << std::endl;
  test_Y(Amanzi::Testing::demoMeshLogicalFromXML("logical mesh Y"), true);
}

// Evaluates a Y-mesh embedded in another mesh
TEST(MESH_EMBEDDED_Y)
{
  std::cout << std::endl
            << "TEST: MeshLogical Y, embedded in background mesh" << std::endl
            << "-----------------------------------------------------" << std::endl;
  auto m = Amanzi::Testing::demoMeshLogicalYEmbedded();
  Amanzi::AmanziMesh::MeshLogicalAudit audit(m, std::cout);
  CHECK(!audit.Verify());
}

// subgrid model
TEST(MESH_SUBGRID_VARIABLE_TAU)
{
  std::cout << std::endl
            << "TEST: subgrid mesh in travel time space" << std::endl
            << "-----------------------------------------------------" << std::endl;
  auto m = Amanzi::Testing::demoMeshLogicalFromXML("subgrid mesh");
  Amanzi::AmanziMesh::MeshLogicalAudit audit(m, std::cout);
  CHECK(!audit.Verify());
}
