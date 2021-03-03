#include <UnitTest++.h>

#include <fstream>

#include "../Mesh_MSTK.hh"

#include "Epetra_Map.h"
#include "AmanziComm.hh"

#include "MeshAudit.hh"


// Test for generation of hex mesh distributed over 4 processors

TEST(MSTK_HEX_GEN_3x3x3_4P)
{
  std::vector<Amanzi::AmanziMesh::Entity_ID> faces(6), nodes(8);
  std::vector<Amanzi::AmanziGeometry::Point> ccoords(8), fcoords(4);

  auto comm = Amanzi::getDefaultComm();
  int rank = comm->MyPID();
  int size = comm->NumProc();
  CHECK_EQUAL(4,size);
  if (size != 4) {
    std::cerr << "Test must be run with 4 processors" << std::endl;
    //    return;
  }

  // if (rank == 0) {
  int DebugWait = 0;
  while (DebugWait);
  // }

  Teuchos::RCP<Amanzi::AmanziMesh::MeshFramework> mesh(new Amanzi::AmanziMesh::Mesh_MSTK(0.0,0.0,0.0,1.0,1.0,1.0,3,3,3,comm));

  std::stringstream fname;
  fname << "test/mstk_hex_gen_3x3x3_4P." << rank << ".out";
  std::ofstream fout(fname.str().c_str());
  Amanzi::MeshAudit auditor(mesh,fname);
  auditor.Verify();  
}

