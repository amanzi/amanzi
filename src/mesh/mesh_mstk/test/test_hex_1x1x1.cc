/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

//!

#include <UnitTest++.h>

#include <fstream>


#include "../Mesh_MSTK.hh"

#include "AmanziMap.hh"
#include "AmanziComm.hh"

#include "mpi.h"


TEST(MSTK_HEX1)
{
  int i, j, k, err, nc, nv;
  Kokkos::View<Amanzi::AmanziMesh::Entity_ID*> faces("", 6);
  Amanzi::AmanziMesh::Entity_ID_List cellnodes(8), facenodes(4);
  Kokkos::View<int*> facedirs("", 6);
  std::vector<Amanzi::AmanziGeometry::Point> ccoords(8), fcoords(4);

  int NV = 8;
  int NF = 6;
  int NC = 1;
  double xyz[12][3] = { { 0, 0, 0 }, { 1, 0, 0 }, { 0, 1, 0 }, { 1, 1, 0 },
                        { 0, 0, 1 }, { 1, 0, 1 }, { 0, 1, 1 }, { 1, 1, 1 } };

  auto comm = Amanzi::getDefaultComm();

  // Load a single hex from the hex1.exo file

  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh(
    new Amanzi::AmanziMesh::Mesh_MSTK("test/hex_1x1x1_ss.exo", comm));


  // Check number of nodes and their coordinates

  nv = mesh->num_entities(Amanzi::AmanziMesh::NODE,
                          Amanzi::AmanziMesh::Parallel_type::OWNED);
  CHECK_EQUAL(NV, nv);

  for (i = 0; i < nv; i++) {
    Amanzi::AmanziGeometry::Point coords(mesh->space_dimension());

    //    coords.init(mesh->space_dimension());

    mesh->node_get_coordinates(i, &coords);
    CHECK_ARRAY_EQUAL(xyz[i], coords, 3);
  }


  // Check number of cells and their face nodes and their face coordinates

  nc = mesh->num_entities(Amanzi::AmanziMesh::CELL,
                          Amanzi::AmanziMesh::Parallel_type::OWNED);
  CHECK_EQUAL(NC, nc);


  // Check cell coordinates directly

  mesh->cell_get_nodes(0, cellnodes);
  mesh->cell_get_coordinates(0, ccoords);

  for (j = 0; j < 8; j++) {
    CHECK_ARRAY_EQUAL(xyz[cellnodes[j]], ccoords[j], 3);
  }

  mesh->cell_get_faces_and_dirs(0, faces, facedirs);

  for (j = 0; j < 6; j++) {
    mesh->face_get_nodes(faces(j), facenodes);

    int f0 = facenodes[0];
    int f1 = facenodes[1];
    for (k = 2; k < 4; k++) {
      int f2 = facenodes[k];
      Amanzi::AmanziGeometry::Point p0(xyz[f0][0], xyz[f0][1], xyz[f0][2]);
      Amanzi::AmanziGeometry::Point p1(xyz[f1][0], xyz[f1][1], xyz[f1][2]);
      Amanzi::AmanziGeometry::Point p2(xyz[f2][0], xyz[f2][1], xyz[f2][2]);

      Amanzi::AmanziGeometry::Point cross = (p1 - p0) ^ (p2 - p0);
      double area = Amanzi::AmanziGeometry::norm(cross);
      CHECK_CLOSE(area, 1.0, 1e-10);
    }
  }
}
