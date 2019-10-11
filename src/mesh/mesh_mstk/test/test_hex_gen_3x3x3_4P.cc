#include <UnitTest++.h>

#include <fstream>

#include "../Mesh_MSTK.hh"

#include "AmanziMap.hh"
#include "AmanziComm.hh"

#include "MeshAudit.hh"


// Test for generation of hex mesh distributed over 4 processors

TEST(MSTK_HEX_GEN_3x3x3_4P)
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

  // if (rank == 0) {
  int DebugWait = 0;
  while (DebugWait)
    ;
  // }

  // Create a 3x3x3 cell hex mesh

  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh(new Amanzi::AmanziMesh::Mesh_MSTK(
    0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 3, 3, 3, comm));


  Kokkos::View<Amanzi::AmanziMesh::Entity_ID*> c2f("", 6);
  auto cell_map = mesh->cell_map(false);
  auto face_map = mesh->face_map(true);

  for (int c = cell_map->getMinLocalIndex(); c <= cell_map->getMaxLocalIndex();
       c++) {
    CHECK_EQUAL(cell_map->getGlobalElement(c),
                mesh->getGlobalElement(c, Amanzi::AmanziMesh::CELL));
    mesh->cell_get_faces(c, c2f);

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
  fname << "test/mstk_hex_gen_3x3x3_4P." << rank << ".out";
  std::ofstream fout(fname.str().c_str());
  Amanzi::MeshAudit auditor(mesh, fname);
  auditor.Verify();
}
