#include <UnitTest++.h>
#include <fstream>
#include <fstream>

#include "../Mesh_MSTK.hh"
#include "MeshAudit.hh"

#include "Epetra_Map.h"
#include "AmanziComm.hh"



TEST(MSTK_HEX_3x3x3_4P)
{
  std::vector<Amanzi::AmanziMesh::Entity_ID> faces(6), nodes(8);
  std::vector<Amanzi::AmanziGeometry::Point> ccoords(8), fcoords(4);

  auto comm = Amanzi::getDefaultComm();
  int rank = comm->MyPID();
  int size = comm->NumProc();
  CHECK_EQUAL(4,size);

  if (size != 4) {
    std::cerr << "Test must be run with 4 processors" << std::endl;
  }

  // Load a single hex from the hex1.exo file
  Teuchos::RCP<Amanzi::AmanziMesh::MeshFramework> mesh(new Amanzi::AmanziMesh::Mesh_MSTK("test/hex_3x3x3_sets.exo",comm));

  std::stringstream fname;
  fname << "test/mstk_hex_3x3x3_4P." << rank << ".out";
  std::ofstream fout(fname.str().c_str());
  Amanzi::MeshAudit auditor(mesh,fout);
  auditor.Verify();
}

