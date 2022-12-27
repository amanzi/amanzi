/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include <UnitTest++.h>

#include <mpi.h>
#include <iostream>

#include "Teuchos_RCP.hpp"
#include "AmanziComm.hh"

#include "RegionEnumerated.hh"
#include "MeshLogical.hh"
#include "MeshLogicalFactory.hh"
#include "Geometry.hh"

#include "demo_mesh.hh"

void
test_segment_regular(const Teuchos::RCP<Amanzi::AmanziMesh::Mesh>& m,
                     bool test_region) {
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;

  CHECK_EQUAL(4, m->getNumEntities(CELL, Parallel_type::ALL));
  CHECK_EQUAL(5, m->getNumEntities(FACE, Parallel_type::ALL));
  CHECK_EQUAL(0.25, m->getCellVolume(0));
  CHECK_EQUAL(1.0, m->getFaceArea(0));

  CHECK_EQUAL(1.0, m->getFaceNormal(3)[0]);
  CHECK_EQUAL(0.0, m->getFaceNormal(3)[1]);
  CHECK_EQUAL(0.0, m->getFaceNormal(3)[2]);

  Entity_ID_List faces;
  std::vector<int> dirs;
  std::vector<Point> bisectors;
  m->getCellFacesAndDirs(2, faces, &dirs);
  CHECK_EQUAL(2, faces.size());
  CHECK_EQUAL(2, faces[0]);
  CHECK_EQUAL(3, faces[1]);
  CHECK_EQUAL(2, dirs.size());
  CHECK_EQUAL(-1, dirs[0]);
  CHECK_EQUAL(1, dirs[1]);

  faces.clear();
  m->getCellFacesAndBisectors(2, faces, &bisectors);
  CHECK_EQUAL(2, faces.size());
  CHECK_EQUAL(2, faces[0]);
  CHECK_EQUAL(3, faces[1]);
  CHECK_EQUAL(2, bisectors.size());
  CHECK_EQUAL(0.125, bisectors[0][0]);
  CHECK_EQUAL(0., bisectors[0][1]);
  CHECK_EQUAL(0., bisectors[0][2]);
  CHECK_EQUAL(-0.125, bisectors[1][0]);
  CHECK_EQUAL(0., bisectors[1][1]);
  CHECK_EQUAL(0., bisectors[1][2]);

  Entity_ID_List cells;
  m->getFaceCells(0, Parallel_type::ALL, cells);
  CHECK_EQUAL(1, cells.size());
  CHECK_EQUAL(0, cells[0]);

  m->getFaceCells(1, Parallel_type::ALL, cells);
  CHECK_EQUAL(2, cells.size());
  CHECK_EQUAL(0, cells[0]);
  CHECK_EQUAL(1, cells[1]);

  m->getFaceCells(4, Parallel_type::ALL, cells);
  CHECK_EQUAL(1, cells.size());
  CHECK_EQUAL(3, cells[0]);

  
  if (test_region) {
    // check regions
    CHECK_EQUAL(4, m->getSetSize("myregion", CELL, Parallel_type::ALL));
    CHECK_EQUAL(0, m->getSetSize("myregion", FACE, Parallel_type::ALL));

    Entity_ID_List set_ents;
    set_ents = m->getSetEntities("myregion", CELL, Parallel_type::ALL);
    CHECK_EQUAL(0, set_ents[0]);
    CHECK_EQUAL(2, set_ents[2]);
  }
}  


void
test_segment_irregular(const Teuchos::RCP<Amanzi::AmanziMesh::Mesh>& m,
                     bool test_region) {
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;

  Teuchos::RCP<const GeometricModel> gm_c = m->geometric_model();
  Teuchos::RCP<GeometricModel> gm = Teuchos::rcp_const_cast<GeometricModel>(gm_c);
  
  Entity_ID_List ents;
  ents.push_back(0);
  ents.push_back(3);

  Teuchos::RCP<Region> enum_rgn =
    Teuchos::rcp(new RegionEnumerated("myregion", 0, "CELL", ents));
  gm->AddRegion(enum_rgn);


  CHECK_EQUAL(2, m->getSetSize("myregion", CELL, Parallel_type::ALL));
  CHECK_EQUAL(0, m->getSetSize("myregion", FACE, Parallel_type::ALL));

  if (test_region) {
    Entity_ID_List set_ents;
    set_ents = m->getSetEntities("myregion", CELL, Parallel_type::ALL);
    CHECK_EQUAL(0, set_ents[0]);
    CHECK_EQUAL(3, set_ents[1]);
  }
}


void
test_Y(const Teuchos::RCP<Amanzi::AmanziMesh::Mesh>& m,
                     bool test_region) {
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;

  // surface
  Point zero(0.,0.,0.);
  CHECK_CLOSE(0., norm(zero-m->getFaceCentroid(0)), 1.e-6);

  // branch point
  Point branch(0.,0.,-1.25);
  CHECK_CLOSE(0., norm(branch - m->getCellCentroid(2)), 1.e-6);

  Entity_ID_List branch_faces;
  std::vector<int> dirs;
  m->getCellFacesAndDirs(2, branch_faces, &dirs);
  CHECK_EQUAL(5, branch_faces.size());

  // fine root 1 tip
  Point fine1_normal(1.,0.,-1.);
  fine1_normal /= norm(fine1_normal);
  Point tip = branch + 1.25*std::sqrt(2) * fine1_normal;
  CHECK_CLOSE(0., norm(tip-m->getFaceCentroid(4)), 1.e-6);

  if (test_region) {
    CHECK_EQUAL(3, m->getSetSize("coarse_root", CELL, Parallel_type::ALL));
    CHECK_EQUAL(8, m->getSetSize("fine_root", CELL, Parallel_type::ALL));
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
  test_segment_regular(Amanzi::Testing::demoMeshLogicalSegmentRegularManual(), false);
}

// Evaluates the mesh constructed through the factory.
TEST(MESH_LOGICAL_SEGMENT_REGULAR_FACTORY)
{
  test_segment_regular(Amanzi::Testing::demoMeshLogicalSegmentRegular(), true);
  CHECK(*Amanzi::Testing::demoMeshLogicalSegmentRegular() == *Amanzi::Testing::demoMeshLogicalSegmentRegularManual());
}

// Evaluates an irregularly space mesh
TEST(MESH_LOGICAL_SEGMENT_IRREGULAR_WITH_SETS)
{
  test_segment_irregular(Amanzi::Testing::demoMeshLogicalSegmentIrregularManual(), true);
}


// Evaluates a Y-mesh
TEST(MESH_LOGICAL_Y_MANUAL_WITH_SETS)
{
  test_Y(Amanzi::Testing::demoMeshLogicalYManual(), true);
}

// Evaluates a Y-mesh
TEST(MESH_LOGICAL_Y_WITH_SETS)
{
  test_Y(Amanzi::Testing::demoMeshLogicalY(), true);
  CHECK(*Amanzi::Testing::demoMeshLogicalY() == *Amanzi::Testing::demoMeshLogicalYManual());
}


// Evaluates a Y-mesh
TEST(MESH_LOGICAL_Y_XML_WITH_SETS)
{
  test_Y(Amanzi::Testing::demoMeshLogicalYFromXML(), true);
  CHECK(*Amanzi::Testing::demoMeshLogicalYFromXML() == *Amanzi::Testing::demoMeshLogicalYManual());
}


// Evaluates a Y-mesh
TEST(MESH_EMBEDDED_Y)
{
  Amanzi::Testing::demoMeshLogicalYEmbedded();
}
