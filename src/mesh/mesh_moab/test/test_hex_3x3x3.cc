#include <iostream>

#include "Epetra_Map.h"
#include "AmanziComm.hh"
#include "mpi.h"
#include "UnitTest++.h"

#include "../Mesh_MOAB.hh"

// Unless this example is enhanced, it does lesser testing than test_hex_3x3x2.cc

TEST(MOAB_HEX_3x3x3)
{
  using namespace Amanzi;

  int j, nc, nf;
  AmanziMesh::Entity_ID_List faces, nodes;
  AmanziGeometry::Point ccoords, fcoords;

  int NV = 64;
  int NF = 108;
  int NC = 27;

  auto comm = Amanzi::getDefaultComm();

  AmanziMesh::Mesh_MOAB mesh("test/hex_3x3x3_ss.exo",comm);

  nf = mesh.num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  CHECK_EQUAL(NF, nf);
  
  nc = mesh.num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  CHECK_EQUAL(NC, nc);

  AmanziMesh::Entity_ID_List c2f;
  Epetra_Map cell_map(mesh.cell_map(false));
  Epetra_Map face_map(mesh.face_map(false));
  for (int c = cell_map.MinLID(); c <= cell_map.MaxLID(); c++) {
    mesh.cell_get_faces( c, &c2f, true);
    for (j = 0; j < 6; j++) {
       int f = face_map.LID(c2f[j]);
       CHECK(f == c2f[j]);
    }
  }
  
  // verify boundary maps
  Epetra_Map extface_map(mesh.exterior_face_map(false));
  for (int f = extface_map.MinLID(); f <= extface_map.MaxLID(); ++f) {
    const AmanziGeometry::Point& xf = mesh.face_centroid(f);
    CHECK(std::fabs(xf[0]) < 1e-7 || std::fabs(1.0 - xf[0]) < 1e-7 ||
          std::fabs(xf[1]) < 1e-7 || std::fabs(1.0 - xf[1]) < 1e-7 ||
          std::fabs(xf[2]) < 1e-7 || std::fabs(1.0 - xf[2]) < 1e-7);
  }
}

