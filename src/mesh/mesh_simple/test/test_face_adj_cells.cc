/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

//!

#include <iostream>
#include "stdlib.h"
#include "math.h"
#include "UnitTest++.h"
#include "../Mesh_simple.hh"

#include <AmanziComm.hh>

TEST(FACE_ADJ_CELLS)
{
  using namespace std;
  auto comm = Amanzi::getDefaultComm();
  const unsigned int exp_ncell = 27, exp_nface = 108, exp_nnode = 64;
  const unsigned int exp_nadj[27] = { 3, 4, 3, 4, 5, 4, 3, 4, 3, 4, 5, 4, 5, 6,
                                      5, 4, 5, 4, 3, 4, 3, 4, 5, 4, 3, 4, 3 };
  const int exp_adjcells[27][6] = {
    { 1, 3, 9, -1, -1, -1 },    { 2, 4, 0, 10, -1, -1 },
    { 5, 1, 11, -1, -1, -1 },   { 0, 4, 6, 12, -1, -1 },
    { 1, 5, 7, 3, 13, -1 },     { 2, 8, 4, 14, -1, -1 },
    { 3, 7, 15, -1, -1, -1 },   { 4, 8, 6, 16, -1, -1 },
    { 5, 7, 17, -1, -1, -1 },

    { 10, 12, 0, 18, -1, -1 },  { 11, 13, 9, 1, 19, -1 },
    { 14, 10, 2, 20, -1, -1 },  { 9, 13, 15, 3, 21, -1 },
    { 10, 14, 16, 12, 4, 22 },  { 11, 17, 13, 5, 23, -1 },
    { 12, 16, 6, 24, -1, -1 },  { 13, 17, 15, 7, 25, -1 },
    { 14, 16, 8, 26, -1, -1 },

    { 19, 21, 9, -1, -1, -1 },  { 20, 22, 18, 10, -1, -1 },
    { 23, 19, 11, -1, -1, -1 }, { 18, 22, 24, 12, -1, -1 },
    { 19, 23, 25, 21, 13, -1 }, { 20, 26, 22, 14, -1, -1 },
    { 21, 25, 15, -1, -1, -1 }, { 22, 26, 24, 16, -1, -1 },
    { 23, 25, 17, -1, -1, -1 }
  };


  Amanzi::AmanziMesh::Mesh_simple Mm(
    0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 3, 3, 3, comm);

  CHECK_EQUAL(exp_ncell,
              Mm.num_entities(Amanzi::AmanziMesh::CELL,
                              Amanzi::AmanziMesh::Parallel_type::OWNED));
  CHECK_EQUAL(exp_nface,
              Mm.num_entities(Amanzi::AmanziMesh::FACE,
                              Amanzi::AmanziMesh::Parallel_type::OWNED));
  CHECK_EQUAL(exp_nnode,
              Mm.num_entities(Amanzi::AmanziMesh::NODE,
                              Amanzi::AmanziMesh::Parallel_type::OWNED));


  for (int i = 0; i < exp_ncell; i++) {
    Amanzi::AmanziMesh::Entity_ID_List adjcells;

    Mm.cell_get_face_adj_cells(
      i, Amanzi::AmanziMesh::Parallel_type::OWNED, adjcells);

    unsigned int nadj = adjcells.size();
    CHECK_EQUAL(exp_nadj[i], nadj);

    for (int j = 0; j < nadj; j++) CHECK_EQUAL(exp_adjcells[i][j], adjcells[j]);
  }
}
