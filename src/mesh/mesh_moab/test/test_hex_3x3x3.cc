/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

//!

#include <iostream>

#include "AmanziMap.hh"
#include "AmanziComm.hh"
#include "mpi.h"
#include "UnitTest++.h"

#include "../Mesh_MOAB.hh"

// Unless this example is enhanced, it does lesser testing than
// test_hex_3x3x2.cc

TEST(MOAB_HEX_3x3x3)
{
  using namespace Amanzi;

  int i, j, k, err, nc, nf, nv;
  AmanziMesh::Entity_ID_List faces, nodes;
  AmanziGeometry::Point ccoords, fcoords;

  int NV = 64;
  int NF = 108;
  int NC = 27;

  auto comm = Amanzi::getDefaultComm();

  AmanziMesh::Mesh_MOAB mesh("test/hex_3x3x3_ss.exo", comm);

  nf = mesh.num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  CHECK_EQUAL(NF, nf);

  nc = mesh.num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  CHECK_EQUAL(NC, nc);

  AmanziMesh::Entity_ID_List c2f;
  auto cell_map = mesh.cell_map(false);
  auto face_map = mesh.face_map(false);
  for (int c = cell_map->getMinLocalIndex(); c <= cell_map->getMaxLocalIndex();
       c++) {
    mesh.cell_get_faces(c, &c2f, true);
    for (int j = 0; j < 6; j++) {
      int f = face_map->getLocalElement(c2f[j]);
      CHECK(f == c2f[j]);
    }
  }

  // verify boundary maps
  auto extface_map = mesh.exterior_face_map(false);
  for (int f = extface_map->getMinLocalIndex();
       f <= extface_map->getMaxLocalIndex();
       ++f) {
    const AmanziGeometry::Point& xf = mesh.face_centroid(f);
    CHECK(std::fabs(xf[0]) < 1e-7 || std::fabs(1.0 - xf[0]) < 1e-7 ||
          std::fabs(xf[1]) < 1e-7 || std::fabs(1.0 - xf[1]) < 1e-7 ||
          std::fabs(xf[2]) < 1e-7 || std::fabs(1.0 - xf[2]) < 1e-7);
  }
}
