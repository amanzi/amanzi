/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

#include <UnitTest++.h>

#include <mpi.h>
#include <iostream>

#include "Teuchos_RCP.hpp"
#include "Epetra_MpiComm.h"

#include "EnumeratedSetRegion.hh"
#include "MeshLogical.hh"
#include "MeshLogicalFactory.hh"
#include "Geometry.hh"

Teuchos::RCP<Amanzi::AmanziMesh::MeshLogical> Demo3x3x1_Regular() {
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  Epetra_MpiComm comm(MPI_COMM_WORLD);
  const int nproc(comm.NumProc());
  const int me(comm.MyPID());

  GeometricModelPtr gm = new GeometricModel(3);
  
  // Create the mesh:
  std::vector<double> cell_volumes;
  cell_volumes.push_back(.25);
  cell_volumes.push_back(.25);
  cell_volumes.push_back(.25);
  cell_volumes.push_back(.25);

  std::vector<std::vector<int> > face_cell_list(5);
  face_cell_list[0].push_back(0);
  face_cell_list[1].push_back(0);
  face_cell_list[1].push_back(1);
  face_cell_list[2].push_back(1);
  face_cell_list[2].push_back(2);
  face_cell_list[3].push_back(2);
  face_cell_list[3].push_back(3);
  face_cell_list[4].push_back(3);
  

  std::vector<std::vector<double> > face_cell_lengths;
  face_cell_lengths.resize(5);
  face_cell_lengths[0].push_back(0.125);
  face_cell_lengths[1].resize(2,0.125);
  face_cell_lengths[2].resize(2,0.125);
  face_cell_lengths[3].resize(2,0.125);
  face_cell_lengths[4].push_back(0.125);

  std::vector<Amanzi::AmanziGeometry::Point> face_area_normals;
  Amanzi::AmanziGeometry::Point normal(3);
  normal[0] = 1.0;
  normal[1] = 0.0;                    
  normal[2] = 0.0;
  face_area_normals.resize(5,normal);

  Teuchos::RCP<Amanzi::AmanziMesh::MeshLogical> mesh =
    Teuchos::rcp(new MeshLogical(&comm,cell_volumes,
				 face_cell_list,
				 face_cell_lengths,
				 face_area_normals));
  mesh->set_geometric_model(gm);
  return mesh;    
}

Teuchos::RCP<Amanzi::AmanziMesh::MeshLogical> Demo3x3x1_Regular_Factory() {
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;

  Epetra_MpiComm comm(MPI_COMM_WORLD);

  GeometricModelPtr gm = new GeometricModel(3);
  MeshLogicalFactory fac(&comm, gm);

  Entity_ID_List cells, faces;
  fac.AddSegment(4, 1.0, false, 1.0, true,true, "myregion", &cells, &faces);
  Teuchos::RCP<MeshLogical> mesh = fac.Create();
  return mesh;
}

Teuchos::RCP<Amanzi::AmanziMesh::MeshLogical> Demo3x3x1_Irregular() {
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  const int nproc(comm.NumProc());
  const int me(comm.MyPID());

  GeometricModelPtr gm = new GeometricModel(3);
  
  // Create the mesh:
  std::vector<double> cell_volumes;
  cell_volumes.push_back(2);
  cell_volumes.push_back(3);
  cell_volumes.push_back(4);

  std::vector<std::vector<int> > face_cell_list;
  face_cell_list.resize(4);
  face_cell_list[0].push_back(0);
  face_cell_list[1].push_back(0);
  face_cell_list[1].push_back(1);
  face_cell_list[2].push_back(1);
  face_cell_list[2].push_back(2);
  face_cell_list[3].push_back(2);
  
  std::vector<std::vector<double> > face_cell_lengths;
  face_cell_lengths.resize(4);
  face_cell_lengths[0].push_back(0.5);
  face_cell_lengths[1].push_back(0.5);
  face_cell_lengths[1].push_back(.75);
  face_cell_lengths[2].push_back(.75);
  face_cell_lengths[2].push_back(1.);
  face_cell_lengths[3].push_back(1.);
  
  std::vector<Amanzi::AmanziGeometry::Point> face_area_normals;
  Amanzi::AmanziGeometry::Point normal(3);
  normal[0] = 2.0;
  normal[1] = 0.0;                    
  normal[2] = 0.0;
  face_area_normals.resize(4,normal);

  Teuchos::RCP<Amanzi::AmanziMesh::MeshLogical> mesh =
    Teuchos::rcp(new MeshLogical(&comm,cell_volumes,
				 face_cell_list,
				 face_cell_lengths,
				 face_area_normals));
  mesh->set_geometric_model(gm);
  return mesh;    
}


void
run_test_regular(const Teuchos::RCP<Amanzi::AmanziMesh::Mesh>& m,
                 bool test_region) {
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;

  CHECK_EQUAL(4, m->num_entities(CELL, USED));
  CHECK_EQUAL(5, m->num_entities(FACE, USED));
  CHECK_EQUAL(0.25, m->cell_volume(0));
  CHECK_EQUAL(1.0, m->face_area(0));

  CHECK_EQUAL(1.0, m->face_normal(3)[0]);
  CHECK_EQUAL(0.0, m->face_normal(3)[1]);
  CHECK_EQUAL(0.0, m->face_normal(3)[2]);

  Entity_ID_List faces;
  std::vector<int> dirs;
  std::vector<Point> bisectors;
  m->cell_get_faces_and_dirs(2, &faces, &dirs);
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
  CHECK_EQUAL(0.125, bisectors[0][0]);
  CHECK_EQUAL(0., bisectors[0][1]);
  CHECK_EQUAL(0., bisectors[0][2]);
  CHECK_EQUAL(-0.125, bisectors[1][0]);
  CHECK_EQUAL(0., bisectors[1][1]);
  CHECK_EQUAL(0., bisectors[1][2]);

  Entity_ID_List cells;
  m->face_get_cells(0, USED, &cells);
  CHECK_EQUAL(1, cells.size());
  CHECK_EQUAL(0, cells[0]);

  m->face_get_cells(1, USED, &cells);
  CHECK_EQUAL(2, cells.size());
  CHECK_EQUAL(0, cells[0]);
  CHECK_EQUAL(1, cells[1]);

  m->face_get_cells(4, USED, &cells);
  CHECK_EQUAL(1, cells.size());
  CHECK_EQUAL(3, cells[0]);

  
  if (test_region) {
    // check regions
    CHECK_EQUAL(4, m->get_set_size("myregion", CELL, USED));
    CHECK_EQUAL(0, m->get_set_size("myregion", FACE, USED));

    Entity_ID_List set_ents;
    m->get_set_entities("myregion", CELL, USED, &set_ents);
    CHECK_EQUAL(0, set_ents[0]);
    CHECK_EQUAL(2, set_ents[2]);
  }

}  

TEST(MESH_LOGICAL_CONSTRUCTION)
{
  using namespace Amanzi::AmanziMesh;
  Teuchos::RCP<Amanzi::AmanziMesh::MeshLogical> mesh =
    Demo3x3x1_Irregular();
}

TEST(MESH_LOGICAL_STRUCTURED)
{
  run_test_regular(Demo3x3x1_Regular(), false);
}

TEST(MESH_LOGICAL_STRUCTURED_FACTORY)
{
  run_test_regular(Demo3x3x1_Regular_Factory(), true);
}

TEST(MESH_LOGICAL_WITH_SETS)
{
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  Teuchos::RCP<MeshLogical> mesh =
    Demo3x3x1_Irregular();

  GeometricModelPtr gm = mesh->geometric_model();

  Entity_ID_List ents;
  ents.push_back(0);
  ents.push_back(3);

  Teuchos::RCP<Region> enum_rgn =
    Teuchos::rcp(new EnumeratedSetRegion("myregion",
					   0, "CELL", ents));

  gm->Add_Region(&*enum_rgn);
  mesh->set_geometric_model(&*gm);

  CHECK_EQUAL(2, mesh->get_set_size("myregion", CELL, USED));
  CHECK_EQUAL(0, mesh->get_set_size("myregion", FACE, USED));

  Entity_ID_List set_ents;
  mesh->get_set_entities("myregion", CELL, USED, &set_ents);
  CHECK_EQUAL(0, set_ents[0]);
  CHECK_EQUAL(3, set_ents[1]);
}

