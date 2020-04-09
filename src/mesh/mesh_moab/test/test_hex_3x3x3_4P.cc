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


TEST(MOAB_HEX_3x3x3_4P)
{
  using namespace Amanzi;

  int i, j, k, err, nc, nf, nv;
  unsigned int faces[6], nodes[8];
  int facedirs[6];
  double ccoords[24], fcoords[12];

  int NVowned[4] = { 16, 16, 16, 16 };
  int NFowned[4] = { 16, 26, 26, 40 };
  int NCowned[4] = { 3, 6, 6, 12 };
  int NVused[4] = { 36, 48, 48, 64 };
  int NFused[4] = { 52, 75, 75, 108 };
  int NCused[4] = { 12, 18, 18, 27 };
  int NVghost[4] = { 20, 32, 32, 48 };
  int NFghost[4] = { 36, 49, 49, 68 };
  int NCghost[4] = { 9, 12, 12, 15 };

  auto comm = Amanzi::getDefaultComm();

  int rank, size;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  CHECK_EQUAL(4, size);

  if (rank == 0) {
    int DebugWait = 0;
    while (DebugWait)
      ;
  }

  // Load a single hex from the hex1.exo file

  AmanziMesh::Mesh_MOAB mesh("test/hex_3x3x3_ss_4P.h5m", comm);

  nv = mesh.num_entities(AmanziMesh::NODE, AmanziMesh::Parallel_type::OWNED);
  CHECK_EQUAL(NVowned[rank], nv);

  nf = mesh.num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  CHECK_EQUAL(NFowned[rank], nf);

  nc = mesh.num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  CHECK_EQUAL(NCowned[rank], nc);

  nv = mesh.num_entities(AmanziMesh::NODE, AmanziMesh::Parallel_type::ALL);
  CHECK_EQUAL(NVused[rank], nv);

  nf = mesh.num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);
  CHECK_EQUAL(NFused[rank], nf);

  nc = mesh.num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);
  CHECK_EQUAL(NCused[rank], nc);

  nv = mesh.num_entities(AmanziMesh::NODE, AmanziMesh::Parallel_type::GHOST);
  CHECK_EQUAL(NVghost[rank], nv);

  nf = mesh.num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::GHOST);
  CHECK_EQUAL(NFghost[rank], nf);

  nc = mesh.num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::GHOST);
  CHECK_EQUAL(NCghost[rank], nc);


  AmanziMesh::Entity_ID_List c2f;
  std::vector<int> c2fdirs;
  auto cell_map = mesh.cell_map(false);
  auto face_map = mesh.face_map(true);

  for (int c = cell_map->getMinLocalIndex(); c <= cell_map->getMaxLocalIndex();
       c++) {
    CHECK_EQUAL(cell_map->getGlobalElement(c),
                mesh.getGlobalElement(c, AmanziMesh::CELL));
    mesh.cell_get_faces_and_dirs(c, &c2f, c2fdirs);

    for (int j = 0; j < 6; j++) {
      int f = face_map->getLocalElement(
        mesh.getGlobalElement(c2f[j], AmanziMesh::FACE));
      CHECK_EQUAL(f, c2f[j]);
      CHECK_EQUAL(1, abs(c2fdirs[j]));
    }
  }

  // verify boundary maps: owned
  int gid, g;
  {
    auto extface_map = mesh.exterior_face_map(false);
    int nfaces(extface_map->getMaxLocalIndex() + 1), nall;
    Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, 1, &nfaces, &nall);
    CHECK_EQUAL(nall, 54);

    for (int f = extface_map->getMinLocalIndex();
         f <= extface_map->getMaxLocalIndex();
         ++f) {
      gid = extface_map->getGlobalElement(f);
      g = face_map->getLocalElement(gid);
      const AmanziGeometry::Point& xf = mesh.face_centroid(g);

      CHECK(std::fabs(xf[0]) < 1e-7 || std::fabs(1.0 - xf[0]) < 1e-7 ||
            std::fabs(xf[1]) < 1e-7 || std::fabs(1.0 - xf[1]) < 1e-7 ||
            std::fabs(xf[2]) < 1e-7 || std::fabs(1.0 - xf[2]) < 1e-7);
    }
  }

  // verify boundary maps: owned + ghost
  {
    auto extface_map = mesh.exterior_face_map(true);
    for (int f = extface_map->getMinLocalIndex();
         f <= extface_map->getMaxLocalIndex();
         ++f) {
      gid = extface_map->getGlobalElement(f);
      g = face_map->getLocalElement(gid);
      const AmanziGeometry::Point& xf = mesh.face_centroid(g);

      CHECK(std::fabs(xf[0]) < 1e-7 || std::fabs(1.0 - xf[0]) < 1e-7 ||
            std::fabs(xf[1]) < 1e-7 || std::fabs(1.0 - xf[1]) < 1e-7 ||
            std::fabs(xf[2]) < 1e-7 || std::fabs(1.0 - xf[2]) < 1e-7);
    }
  }
}
