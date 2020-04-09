/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

//!

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
#include <fstream>

#include <AmanziComm.hh>

#include "Geometry.hh"
#include "../Mesh_MSTK.hh"
#include "../MeshColumn.hh"
#include "RegionBox.hh"
#include "GeometricModel.hh"


TEST(COLUMN_MESH_3D)
{
  auto comm = Amanzi::getDefaultComm();
  const int nproc(comm->getSize());
  const int me(comm->getRank());


  int nx = 4, ny = 4, nz = 4;
  double lx = 4, ly = 4, lz = 4;
  int dx = 1.0, dy = 1.0, dz = 1.0;

  // create a geometric model with regions
  Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm =
    Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3));

  Amanzi::AmanziGeometry::Point p0(3), p1(3);
  p0[2] = 2.5;
  p1[0] = 4.;
  p1[1] = 4.;
  p1[2] = 5.;
  auto r0 =
    Teuchos::rcp(new Amanzi::AmanziGeometry::RegionBox("myregion", 0, p0, p1));
  gm->AddRegion(r0);


  Amanzi::AmanziGeometry::Point p2(0., 0., 4.), p3(4., 4., 4.);
  auto r1 =
    Teuchos::rcp(new Amanzi::AmanziGeometry::RegionBox("surface", 0, p2, p3));
  gm->AddRegion(r1);


  // Create the mesh
  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh =
    Teuchos::rcp(new Amanzi::AmanziMesh::Mesh_MSTK(
      0.0, 0.0, 0.0, lx, ly, lz, nx, ny, nz, comm, gm));

  CHECK_EQUAL(mesh->build_columns(), 1);

  // Perturb the nodes above the base layer just a bit
  int nnodes = mesh->num_entities(Amanzi::AmanziMesh::NODE,
                                  Amanzi::AmanziMesh::Parallel_type::OWNED);

  for (int n = 0; n < nnodes; n++) {
    Amanzi::AmanziGeometry::Point xyz(3);
    mesh->node_get_coordinates(n, &xyz);
    xyz[2] += 0.005 * xyz[0] * xyz[1] * xyz[2];
    mesh->node_set_coordinates(n, xyz);
  }
 
  // verify in-going topology
  CHECK_EQUAL(16, mesh->num_columns());
  CHECK_EQUAL(4, mesh->cells_of_column(10).size());
  CHECK_EQUAL(5, mesh->faces_of_column(10).size());

  // Create a column mesh from one of the columns
  Amanzi::AmanziMesh::MeshColumn colmesh(mesh, 10);

  // Verify column mesh topology
  int ncells = colmesh.num_entities(Amanzi::AmanziMesh::CELL,
                                    Amanzi::AmanziMesh::Parallel_type::OWNED);
  CHECK_EQUAL(4, ncells);

  int nfaces = colmesh.num_entities(Amanzi::AmanziMesh::FACE,
                                    Amanzi::AmanziMesh::Parallel_type::OWNED);
  CHECK_EQUAL(5, nfaces);

  nnodes = colmesh.num_entities(Amanzi::AmanziMesh::NODE,
                                Amanzi::AmanziMesh::Parallel_type::OWNED);
  CHECK_EQUAL(20, nnodes);

  for (int j = 0; j < ncells; j++) {
    Kokkos::View<Amanzi::AmanziMesh::Entity_ID*> cfaces;
    Kokkos::View<int*> cfdirs;
    colmesh.cell_get_faces_and_dirs(j, cfaces, cfdirs);

    CHECK_EQUAL(2, cfaces.extent(0));
  }
 
  for (int j = 0; j < nfaces; j++) {
    Kokkos::View<Amanzi::AmanziMesh::Entity_ID*> fcells;
    colmesh.face_get_cells(j, Amanzi::AmanziMesh::Parallel_type::OWNED, fcells);

    if (j == 0) {
      CHECK_EQUAL(1, fcells.extent(0));
    } else if (j == nfaces - 1) {
      CHECK_EQUAL(1, fcells.extent(0));
    } else {
      CHECK_EQUAL(2, fcells.extent(0));
    }
  }


  // Verify column mesh geometry
  // centroid of base face

  Amanzi::AmanziGeometry::Point fcenbase(3);
  fcenbase = colmesh.face_centroid(0);

  // area of base face
  double fareabase = colmesh.face_area(0);

  // Make sure centroids of other faces are stacked up
  // exactly above that of the base face

  for (int j = 1; j < nfaces; j++) {
    Amanzi::AmanziGeometry::Point fcen(3);
    fcen = colmesh.face_centroid(j);
    CHECK_EQUAL(fcenbase[0], fcen[0]);
    CHECK_EQUAL(fcenbase[1], fcen[1]);
  }

  // Make sure the normals of the faces are have only a Z component
  for (int j = 0; j < nfaces; j++) {
    Amanzi::AmanziGeometry::Point normal(3);
    normal = colmesh.face_normal(j);
    CHECK_EQUAL(0.0, normal[0]);
    CHECK_EQUAL(0.0, normal[1]);
    CHECK_EQUAL(1.0, fabs(normal[2]));
  }

  // Make sure centroids of cells are stacked up
  // exactly above that of the base face and that
  // their z value is exactly between the z-values
  // of their corresponding bottom and top faces
  for (int i = 0; i < ncells; i++) {
    Amanzi::AmanziGeometry::Point ccen(3), fcen0(3), fcen1(3);

    ccen = colmesh.cell_centroid(i);

    fcen0 = colmesh.face_centroid(i);
    fcen1 = colmesh.face_centroid(i + 1);

    CHECK_EQUAL(fcenbase[0], ccen[0]);
    CHECK_EQUAL(fcenbase[1], ccen[1]);

    CHECK_EQUAL((fcen0[2] + fcen1[2]) / 2.0, ccen[2]);
  }

  // Verify the volume of cells is computed correctly in spite of
  // the topology being special (lateral faces of the cell are
  // omitted). The volume of each cell should be the area of the
  // base face times the distance between the centroids of the lower
  // face of the cell and the upper face of the cell

  for (int j = 0; j < ncells; j++) {
    Kokkos::View<Amanzi::AmanziMesh::Entity_ID*> cfaces;
    colmesh.cell_get_faces(j, cfaces);

    Amanzi::AmanziGeometry::Point locen(3), hicen(3);
    locen = colmesh.face_centroid(cfaces(0));
    hicen = colmesh.face_centroid(cfaces(1));

    double height = norm(hicen - locen);

    double expvolume = fareabase * height;

    double volume = colmesh.cell_volume(j, false);

    CHECK_CLOSE(expvolume, volume, 1.0e-08);
  }


  // verify that the regions have made it through
  Amanzi::AmanziMesh::Entity_ID_List myregion;
  colmesh.get_set_entities("myregion",
                           Amanzi::AmanziMesh::CELL,
                           Amanzi::AmanziMesh::Parallel_type::ALL,
                           myregion);
  CHECK_EQUAL(2, myregion.size());
  CHECK(colmesh.cell_centroid(myregion[0])[2] >= 2.5);
  CHECK(colmesh.cell_centroid(myregion[1])[2] >= 2.5);
  std::cout << "End first" << std::endl;
}


TEST(COLUMN_MESH_3D_FROM_SURFACE)
{
  auto comm = Amanzi::getDefaultComm();
  const int nproc(comm->getSize());
  const int me(comm->getRank());


  int nx = 4, ny = 4, nz = 4;
  double lx = 4, ly = 4, lz = 4;
  int dx = 1.0, dy = 1.0, dz = 1.0;

  // create a geometric model with regions
  Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm =
    Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3));

  Amanzi::AmanziGeometry::Point p0(3), p1(3);
  p0[2] = 2.5;
  p1[0] = 4.;
  p1[1] = 4.;
  p1[2] = 5.;
  auto r0 =
    Teuchos::rcp(new Amanzi::AmanziGeometry::RegionBox("myregion", 0, p0, p1));
  gm->AddRegion(r0);


  Amanzi::AmanziGeometry::Point p2(0., 0., 4.), p3(4., 4., 4.);
  auto r1 =
    Teuchos::rcp(new Amanzi::AmanziGeometry::RegionBox("surface", 0, p2, p3));
  gm->AddRegion(r1);


  // Create the mesh
  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh =
    Teuchos::rcp(new Amanzi::AmanziMesh::Mesh_MSTK(
      0.0, 0.0, 0.0, lx, ly, lz, nx, ny, nz, comm, gm));
  mesh->build_columns("surface");

  int nnodes = mesh->num_entities(Amanzi::AmanziMesh::NODE,
                                  Amanzi::AmanziMesh::Parallel_type::OWNED);

  // verify in-going topology
  CHECK_EQUAL(16, mesh->num_columns());
  CHECK_EQUAL(4, mesh->cells_of_column(10).size());
  CHECK_EQUAL(5, mesh->faces_of_column(10).size());

  // Create a column mesh from one of the columns
  Amanzi::AmanziMesh::MeshColumn colmesh(mesh, 10);


  // Verify column mesh topology
  int ncells = colmesh.num_entities(Amanzi::AmanziMesh::CELL,
                                    Amanzi::AmanziMesh::Parallel_type::OWNED);
  CHECK_EQUAL(4, ncells);

  int nfaces = colmesh.num_entities(Amanzi::AmanziMesh::FACE,
                                    Amanzi::AmanziMesh::Parallel_type::OWNED);
  CHECK_EQUAL(5, nfaces);

  nnodes = colmesh.num_entities(Amanzi::AmanziMesh::NODE,
                                Amanzi::AmanziMesh::Parallel_type::OWNED);
  CHECK_EQUAL(20, nnodes);

  for (int j = 0; j < ncells; j++) {
    Kokkos::View<Amanzi::AmanziMesh::Entity_ID*> cfaces;
    Kokkos::View<int*> cfdirs;
    colmesh.cell_get_faces_and_dirs(j, cfaces, cfdirs);

    CHECK_EQUAL(2, cfaces.extent(0));
  }

  for (int j = 0; j < nfaces; j++) {
    Kokkos::View<Amanzi::AmanziMesh::Entity_ID*> fcells;
    colmesh.face_get_cells(j, Amanzi::AmanziMesh::Parallel_type::OWNED, fcells);

    if (j == 0) {
      CHECK_EQUAL(1, fcells.extent(0));
    } else if (j == nfaces - 1) {
      CHECK_EQUAL(1, fcells.extent(0));
    } else {
      CHECK_EQUAL(2, fcells.extent(0));
    }
  }


  // Verify column mesh geometry
  // centroid of base face

  Amanzi::AmanziGeometry::Point fcenbase(3);
  fcenbase = colmesh.face_centroid(0);

  // area of base face
  double fareabase = colmesh.face_area(0);

  // Make sure centroids of other faces are stacked up
  // exactly above that of the base face

  for (int j = 1; j < nfaces; j++) {
    Amanzi::AmanziGeometry::Point fcen(3);
    fcen = colmesh.face_centroid(j);
    CHECK_EQUAL(fcenbase[0], fcen[0]);
    CHECK_EQUAL(fcenbase[1], fcen[1]);
  }

  // Make sure the normals of the faces are have only a Z component
  for (int j = 0; j < nfaces; j++) {
    Amanzi::AmanziGeometry::Point normal(3);
    normal = colmesh.face_normal(j);
    CHECK_EQUAL(0.0, normal[0]);
    CHECK_EQUAL(0.0, normal[1]);
    CHECK_EQUAL(1.0, fabs(normal[2]));
  }

  // Make sure centroids of cells are stacked up
  // exactly above that of the base face and that
  // their z value is exactly between the z-values
  // of their corresponding bottom and top faces
  for (int i = 0; i < ncells; i++) {
    Amanzi::AmanziGeometry::Point ccen(3), fcen0(3), fcen1(3);

    ccen = colmesh.cell_centroid(i);

    fcen0 = colmesh.face_centroid(i);
    fcen1 = colmesh.face_centroid(i + 1);

    CHECK_EQUAL(fcenbase[0], ccen[0]);
    CHECK_EQUAL(fcenbase[1], ccen[1]);

    CHECK_EQUAL((fcen0[2] + fcen1[2]) / 2.0, ccen[2]);
  }

  // Verify the volume of cells is computed correctly in spite of
  // the topology being special (lateral faces of the cell are
  // omitted). The volume of each cell should be the area of the
  // base face times the distance between the centroids of the lower
  // face of the cell and the upper face of the cell

  for (int j = 0; j < ncells; j++) {
    Kokkos::View<Amanzi::AmanziMesh::Entity_ID*> cfaces;
    colmesh.cell_get_faces(j, cfaces);

    Amanzi::AmanziGeometry::Point locen(3), hicen(3);
    locen = colmesh.face_centroid(cfaces(0));
    hicen = colmesh.face_centroid(cfaces(1));

    double height = norm(hicen - locen);

    double expvolume = fareabase * height;

    double volume = colmesh.cell_volume(j, false);

    CHECK_CLOSE(expvolume, volume, 1.0e-08);
  }


  // verify that the regions have made it through
  Amanzi::AmanziMesh::Entity_ID_List myregion;
  colmesh.get_set_entities("myregion",
                           Amanzi::AmanziMesh::CELL,
                           Amanzi::AmanziMesh::Parallel_type::ALL,
                           myregion);
  CHECK_EQUAL(2, myregion.size());
  CHECK(colmesh.cell_centroid(myregion[0])[2] >= 2.5);
  CHECK(colmesh.cell_centroid(myregion[1])[2] >= 2.5);
}
