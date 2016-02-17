/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
// -------------------------------------------------------------
/**
 * @file   test_column_mesh.cc
 * @author Rao V. Garimella
 * @date   Thu Mar 18, 2015
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------
// -------------------------------------------------------------

#include <UnitTest++.h>

#include <mpi.h>
#include <iostream>

#include <Epetra_MpiComm.h>

#include "Geometry.hh"
#include "../Mesh_MSTK.hh"
#include "../MeshColumn.hh"
#include "../MeshSurfaceCell.hh"
#include "RegionBox.hh"
#include "RegionLabeledSet.hh"
#include "GeometricModel.hh"


TEST(SURFACE_COLUMN_MESH_3D)
{

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  const int nproc(comm.NumProc());
  const int me(comm.MyPID());


  int nx = 4, ny = 4, nz = 4;
  double lx = 4, ly = 4, lz = 4;
  int dx = 1.0, dy = 1.0, dz = 1.0;

  // create a geometric model with surface region
  Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm =
      Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3));

  Amanzi::AmanziGeometry::Point p0(3), p1(3);
  p0[2] = 3.999;
  p1[0] = 4.; p1[1] = 4.; p1[2] = 5.;
  Teuchos::RCP<Amanzi::AmanziGeometry::RegionBox> r0 =
      Teuchos::rcp(new Amanzi::AmanziGeometry::RegionBox("surface", 0, p0, p1));
  gm->AddRegion(r0);

  // Create the mesh
  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh =
      Teuchos::rcp(new Amanzi::AmanziMesh::Mesh_MSTK(0.0,0.0,0.0,
              lx,ly,lz,nx,ny,nz,
              &comm, gm));

  // Perturb the nodes above the base layer just a bit
  int nnodes = mesh->num_entities(Amanzi::AmanziMesh::NODE,
          Amanzi::AmanziMesh::OWNED);
    
  for (int n = 0; n < nnodes; n++) {
    Amanzi::AmanziGeometry::Point xyz(3);
    mesh->node_get_coordinates(n,&xyz);
    xyz[2] += 0.005*xyz[0]*xyz[1]*xyz[2];
    mesh->node_set_coordinates(n,xyz);
  }

  // Create a column mesh from one of the columns
  Amanzi::AmanziMesh::MeshColumn colmesh(*mesh,10);

  // Extract the surface from this column
  Amanzi::AmanziMesh::MeshSurfaceCell col_surf(colmesh, "surface");

  CHECK_EQUAL(1, col_surf.num_entities(Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::OWNED));
  CHECK_EQUAL(4, col_surf.num_entities(Amanzi::AmanziMesh::FACE, Amanzi::AmanziMesh::OWNED));
  CHECK_EQUAL(4, col_surf.num_entities(Amanzi::AmanziMesh::NODE, Amanzi::AmanziMesh::OWNED));

  Amanzi::AmanziMesh::Entity_ID_List cells_in_surf;
  col_surf.get_set_entities("surface", Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::OWNED, &cells_in_surf);
  CHECK_EQUAL(1, cells_in_surf.size());
  CHECK_EQUAL(0, cells_in_surf[0]);

  CHECK_CLOSE(1.0, col_surf.cell_volume(0), 1.e-9);
  CHECK_CLOSE(1.0, col_surf.face_area(3), 1.e-9);
  
}

TEST(SURFACE_COLUMN_MESH_3D_UNSTRUCTURED)
{

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  const int nproc(comm.NumProc());
  const int me(comm.MyPID());

  // create a geometric model with surface region
  Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm =
      Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3));

  Teuchos::RCP<Amanzi::AmanziGeometry::RegionLabeledSet> r1 =
      Teuchos::rcp(new Amanzi::AmanziGeometry::RegionLabeledSet("surface", 0, "FACE",
              "test/slab-0.05-5x4x25.exo", "Exodus II", "1"));
  gm->AddRegion(r1);

  // Create the mesh
  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh =
      Teuchos::rcp(new Amanzi::AmanziMesh::Mesh_MSTK("test/slab-0.05-5x4x25.exo", &comm, 3, gm));

  CHECK_EQUAL(20, mesh->get_set_size("surface",Amanzi::AmanziMesh::FACE, Amanzi::AmanziMesh::USED));
  
  // Create a column mesh from one of the columns
  Amanzi::AmanziMesh::MeshColumn colmesh(*mesh,10);
  CHECK_EQUAL(1, colmesh.get_set_size("surface",Amanzi::AmanziMesh::FACE, Amanzi::AmanziMesh::USED));

  // Extract the surface from this column
  Amanzi::AmanziMesh::MeshSurfaceCell col_surf(colmesh, "surface");

  CHECK_EQUAL(1, col_surf.num_entities(Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::OWNED));
  CHECK_EQUAL(4, col_surf.num_entities(Amanzi::AmanziMesh::FACE, Amanzi::AmanziMesh::OWNED));
  CHECK_EQUAL(4, col_surf.num_entities(Amanzi::AmanziMesh::NODE, Amanzi::AmanziMesh::OWNED));

  Amanzi::AmanziMesh::Entity_ID_List cells_in_surf;
  col_surf.get_set_entities("surface", Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::OWNED, &cells_in_surf);
  CHECK_EQUAL(1, cells_in_surf.size());
  CHECK_EQUAL(0, cells_in_surf[0]);

  CHECK_CLOSE(6400.0, col_surf.cell_volume(0), 1.e-9);
  CHECK_CLOSE(80.0, col_surf.face_area(3), 1.e-9);
  
}
