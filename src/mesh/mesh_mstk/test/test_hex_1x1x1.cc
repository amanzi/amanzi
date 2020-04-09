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
  std::vector<Amanzi::AmanziMesh::Entity_ID> expfacenodes(4);
  Kokkos::View<Amanzi::AmanziMesh::Entity_ID*> faces("", 6);
  Amanzi::AmanziMesh::Entity_ID_List cellnodes(8), facenodes(4);
  Kokkos::View<int*> facedirs("", 6);
  std::vector<Amanzi::AmanziGeometry::Point> ccoords(8), fcoords(4);

  int NV = 8;
  int NF = 6;
  int NC = 1;
  double xyz[12][3] = { { 0, 0, 0 }, { 1, 0, 0 }, { 0, 1, 0 }, { 1, 1, 0 },
                        { 0, 0, 1 }, { 1, 0, 1 }, { 0, 1, 1 }, { 1, 1, 1 } };
  Amanzi::AmanziMesh::Entity_ID local_cellnodes[8] = { 0, 1, 2, 3, 4, 5, 6, 7 };
  Amanzi::AmanziMesh::Entity_ID local_facenodes[6][4] = {
    { 0, 1, 5, 4 }, { 1, 2, 6, 5 }, { 2, 3, 7, 6 },
    { 3, 0, 4, 7 }, { 0, 3, 2, 1 }, { 4, 5, 6, 7 }
  };


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
    mesh->face_get_coordinates(faces(j), fcoords);


    for (k = 0; k < 4; k++) expfacenodes[k] = cellnodes[local_facenodes[j][k]];

    // The order of nodes returned may be different from what we expected
    // So make sure we have a matching node to start with

    int k0 = -1;
    int found = 0;
    for (k = 0; k < 4; k++) {
      if (expfacenodes[k] == facenodes[0]) {
        k0 = k;
        found = 1;
        break;
      }
    }

    CHECK_EQUAL(found, 1);

    if (facedirs(j) == 1) {
      for (k = 0; k < 4; k++) {
        CHECK_EQUAL(expfacenodes[(k0 + k) % 4], facenodes[k]);
        CHECK_ARRAY_EQUAL(xyz[expfacenodes[(k0 + k) % 4]], fcoords[k], 3);
      }
    } else {
      for (k = 0; k < 4; k++) {
        CHECK_EQUAL(expfacenodes[(k0 + 4 - k) % 4], facenodes[k]);
        CHECK_ARRAY_EQUAL(xyz[expfacenodes[(k0 + 4 - k) % 4]], fcoords[k], 3);
      }
    }
  }
}
