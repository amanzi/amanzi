/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

//!

#include <UnitTest++.h>

#include <iostream>

#include "../Mesh_MOAB.hh"


#include "AmanziMap.hh"
#include "AmanziComm.hh"

#include "mpi.h"


TEST(MOAB_HEX_2x2x1)
{
  int i, j, k, err, nc, nf, nv;
  Amanzi::AmanziMesh::Entity_ID_List faces, cnodes, fnodes;
  std::vector<int> facedirs;
  std::vector<Amanzi::AmanziGeometry::Point> ccoords, fcoords;

  int NV = 18;
  int NF = 20;
  int NC = 4;

  double xyz[18][3] = {
    { -0.5, -0.5, 0.25 }, { -0.5, -0.5, -0.25 }, { -0.5, 0, -0.25 },
    { -0.5, 0, 0.25 },    { 0, -0.5, 0.25 },     { 0, -0.5, -0.25 },
    { 0, 0, -0.25 },      { 0, 0, 0.25 },        { -0.5, 0.5, -0.25 },
    { -0.5, 0.5, 0.25 },  { 0, 0.5, -0.25 },     { 0, 0.5, 0.25 },
    { 0.5, -0.5, 0.25 },  { 0.5, -0.5, -0.25 },  { 0.5, 0, -0.25 },
    { 0.5, 0, 0.25 },     { 0.5, 0.5, -0.25 },   { 0.5, 0.5, 0.25 }
  };
  unsigned int cellnodes[4][8] = { { 0, 1, 2, 3, 4, 5, 6, 7 },
                                   { 3, 2, 8, 9, 7, 6, 10, 11 },
                                   { 4, 5, 6, 7, 12, 13, 14, 15 },
                                   { 7, 6, 10, 11, 15, 14, 16, 17 } };
  unsigned int cellfaces[4][6] = { { 8, 4, 16, 0, 10, 17 },
                                   { 16, 5, 12, 1, 11, 18 },
                                   { 9, 6, 19, 2, 17, 14 },
                                   { 19, 7, 13, 3, 18, 15 } };
  int cellfacedirs[4][6] = { { 1, 1, 1, 1, 1, 1 },
                             { -1, 1, 1, 1, 1, 1 },
                             { 1, 1, 1, 1, -1, 1 },
                             { -1, 1, 1, 1, -1, 1 } };
  unsigned int facenodes[20][4] = {
    { 3, 0, 4, 7 },     { 9, 3, 7, 11 },    { 7, 4, 12, 15 },
    { 11, 7, 15, 17 },  { 1, 2, 6, 5 },     { 2, 8, 10, 6 },
    { 5, 6, 14, 13 },   { 6, 10, 16, 14 },  { 0, 1, 5, 4 },
    { 4, 5, 13, 12 },   { 0, 3, 2, 1 },     { 3, 9, 8, 2 },
    { 8, 9, 11, 10 },   { 10, 11, 17, 16 }, { 12, 13, 14, 15 },
    { 15, 14, 16, 17 }, { 2, 3, 7, 6 },     { 4, 5, 6, 7 },
    { 7, 6, 10, 11 },   { 6, 7, 15, 14 }
  };

  unsigned int cfstd[6][4] = {
    { 0, 1, 5, 4 }, // Expected cell-face-node pattern
    { 1, 2, 6, 5 }, { 2, 3, 7, 6 }, { 3, 0, 4, 7 },
    { 0, 3, 2, 1 }, { 4, 5, 6, 7 }
  };

  auto comm = Amanzi::getDefaultComm();


  // Load four hexes from the ExodusII file
  Amanzi::AmanziMesh::Mesh_MOAB mesh("test/hex_2x2x1_ss.exo", comm);


  nv = mesh.num_entities(Amanzi::AmanziMesh::NODE,
                         Amanzi::AmanziMesh::Parallel_type::OWNED);
  CHECK_EQUAL(NV, nv);

  for (i = 0; i < nv; i++) {
    Amanzi::AmanziGeometry::Point coords;

    mesh.node_get_coordinates(i, &coords);
    CHECK_ARRAY_EQUAL(xyz[i], coords, 3);
  }


  nf = mesh.num_entities(Amanzi::AmanziMesh::FACE,
                         Amanzi::AmanziMesh::Parallel_type::OWNED);
  CHECK_EQUAL(NF, nf);

  nc = mesh.num_entities(Amanzi::AmanziMesh::CELL,
                         Amanzi::AmanziMesh::Parallel_type::OWNED);
  CHECK_EQUAL(NC, nc);

  for (i = 0; i < nc; i++) {
    mesh.cell_get_nodes(i, &cnodes);
    mesh.cell_get_faces_and_dirs(i, &faces, facedirs);

    CHECK_ARRAY_EQUAL(cellfaces[i], faces, 6);
    CHECK_ARRAY_EQUAL(cellfacedirs[i], facedirs, 6);

    for (j = 0; j < 6; j++) {
      mesh.face_get_nodes(faces[j], &fnodes);
      mesh.face_get_coordinates(faces[j], &fcoords);

      for (k = 0; k < 4; k++) {
        CHECK_EQUAL(facenodes[faces[j]][k], fnodes[k]);
        CHECK_ARRAY_EQUAL(xyz[facenodes[faces[j]][k]], &(fcoords[k][0]), 3);
      }
    }

    mesh.cell_get_coordinates(i, &ccoords);

    for (j = 0; j < 8; j++) {
      CHECK_EQUAL(cellnodes[i][j], cnodes[j]);
      CHECK_ARRAY_EQUAL(xyz[cellnodes[i][j]], &(ccoords[j][0]), 3);
    }
  }

  // verify global IDs
  Amanzi::AmanziMesh::Entity_ID_List c2f;
  auto cell_map = mesh.cell_map(true);
  auto face_map = mesh.face_map(false);

  for (int c = cell_map->getMinLocalIndex(); c <= cell_map->getMaxLocalIndex();
       c++) {
    CHECK_EQUAL(cell_map->getGlobalElement(c),
                mesh.getGlobalElement(c, Amanzi::AmanziMesh::CELL));

    mesh.cell_get_faces(c, &c2f, true);
    for (int j = 0; j < 6; ++j) {
      int f = face_map->getLocalElement(
        mesh.getGlobalElement(c2f[j], Amanzi::AmanziMesh::FACE));
      CHECK(f == c2f[j]);
    }
  }
}
