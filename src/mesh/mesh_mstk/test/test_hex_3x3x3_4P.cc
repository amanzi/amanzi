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
#include <fstream>

#include "../Mesh_MSTK.hh"
#include "MeshAudit.hh"

#include "AmanziMap.hh"
#include "AmanziComm.hh"


TEST(MSTK_HEX_3x3x3_4P)
{
  int i, j, k, err, nc, nf, nv;
  std::vector<Amanzi::AmanziMesh::Entity_ID> faces(6), nodes(8);
  std::vector<Amanzi::AmanziGeometry::Point> ccoords(8), fcoords(4);

  auto comm = Amanzi::getDefaultComm();
  int rank = comm->getRank();
  int size = comm->getSize();
  CHECK_EQUAL(4, size);

  if (size != 4) {
    std::cerr << "Test must be run with 4 processors" << std::endl;
    //    return;
  }

  //  if (rank == 0) {
  int DebugWait = 0;
  while (DebugWait)
    ;
  //  }

  // Load a single hex from the hex1.exo file

  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh(
    new Amanzi::AmanziMesh::Mesh_MSTK("test/hex_3x3x3_sets.exo", comm));


  Kokkos::View<Amanzi::AmanziMesh::Entity_ID*,Kokkos::HostSpace> c2f("", 6);
  auto cell_map = mesh->cell_map(false);
  auto face_map = mesh->face_map(true);

  for (int c = cell_map->getMinLocalIndex(); c <= cell_map->getMaxLocalIndex();
       c++) {
    CHECK_EQUAL(cell_map->getGlobalElement(c),
                mesh->getGlobalElement(c, Amanzi::AmanziMesh::CELL));
    mesh->cell_get_faces<Kokkos::HostSpace>(c, c2f);

    for (int j = 0; j < 6; j++) {
      int f = face_map->getLocalElement(
        mesh->getGlobalElement(c2f(j), Amanzi::AmanziMesh::FACE));
      CHECK_EQUAL(f, c2f(j));
      if (f != c2f(j)) {
        std::cout << std::endl;
        std::cout << "Processor ID " << rank << std::endl;
        std::cout << "Cell ID " << cell_map->getGlobalElement(c) << std::endl;
        std::cout << "Problem face c2f[j] = " << c2f(j) << " GID = "
                  << mesh->getGlobalElement(c2f(j), Amanzi::AmanziMesh::FACE)
                  << " f = " << f << std::endl;
        std::cout << std::endl;
      }
    }
  }


  std::stringstream fname;
  fname << "test/mstk_hex_3x3x3_4P." << rank << ".out";
  std::ofstream fout(fname.str().c_str());
  Amanzi::MeshAudit auditor(mesh, fout);
  auditor.Verify();
}
