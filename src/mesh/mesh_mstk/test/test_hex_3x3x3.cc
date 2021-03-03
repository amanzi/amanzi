#include <UnitTest++.h>

#include <fstream>

#include "../Mesh_MSTK.hh"

#include "MeshAudit.hh"

#include "Epetra_Map.h"
#include "AmanziComm.hh"


TEST(MSTK_HEX_3x3x3)
{
  int j, nc, nf;
  // int NV = 64;
  int NF = 108;
  int NC = 27;

  auto comm = Amanzi::getDefaultComm();

  // Load a mesh consisting of 3x3x3 elements 

  Teuchos::RCP<Amanzi::AmanziMesh::MeshFramework> mesh(new Amanzi::AmanziMesh::Mesh_MSTK("test/hex_3x3x3_sets.exo",comm));

  nf = mesh->getNumEntities(Amanzi::AmanziMesh::Entity_kind::FACE,Amanzi::AmanziMesh::Parallel_type::OWNED);
  CHECK_EQUAL(NF,nf);
  nc = mesh->getNumEntities(Amanzi::AmanziMesh::Entity_kind::CELL,Amanzi::AmanziMesh::Parallel_type::OWNED);
  CHECK_EQUAL(NC,nc);

  std::stringstream fname;
  fname << "test/mstk_hex_3x3x3.out";
  std::ofstream fout(fname.str().c_str());
  Amanzi::MeshAudit auditor(mesh,fout);
  auditor.Verify();
}

