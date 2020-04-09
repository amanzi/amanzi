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
#include "Epetra_SerialComm.h"

TEST(NODE_CELL_FACES)
{
  using namespace std;
  auto comm = Amanzi::getDefaultComm();
  const unsigned int exp_nnode = 27;

  Amanzi::AmanziMesh::Mesh_simple Mm(
    0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2, 2, 2, comm);


  for (int i = 0; i < exp_nnode; i++) {
    Amanzi::AmanziMesh::Entity_ID node = i;

    Amanzi::AmanziMesh::Entity_ID_List cells;

    Mm.node_get_cells(node, Amanzi::AmanziMesh::Parallel_type::OWNED, cells);

    unsigned int ncells = cells.size();

    for (int j = 0; j < ncells; j++) {
      Amanzi::AmanziMesh::Entity_ID cell = cells[j];

      Amanzi::AmanziMesh::Entity_ID_List faces;

      Mm.node_get_cell_faces(
        node, cell, Amanzi::AmanziMesh::Parallel_type::OWNED, faces);

      // This is a hex mesh. In any given cell, number of faces
      // connected to a node should be 3

      CHECK_EQUAL(3, faces.size());

      for (int k = 0; k < 3; k++) {
        Amanzi::AmanziMesh::Entity_ID face = faces[k];

        Amanzi::AmanziMesh::Entity_ID_List fnodes;

        Mm.face_get_nodes(face, fnodes);

        unsigned int nfnodes = fnodes.size();

        unsigned int found = 0;

        for (int n = 0; n < nfnodes; n++) {
          if (fnodes[n] == node) {
            found = 1;
            break;
          }
        }

        CHECK_EQUAL(1, found);
      }
    }
  }
}
