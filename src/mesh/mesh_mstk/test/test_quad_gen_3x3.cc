#include <UnitTest++.h>

#include <fstream>

#include "../Mesh_MSTK.hh"


#include "Epetra_Map.h"
#include "AmanziComm.hh"

#include "MeshAudit.hh"

// Test generation of quad mesh in serial

TEST(MSTK_QUAD_GEN_3x3)
{
  int j, nc, nf, nv;
  int NV = 16;
  int NF = 24;
  int NC = 9;

  auto comm = Amanzi::getDefaultComm();


  // Load a mesh consisting of 3x3 elements

  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh(new Amanzi::AmanziMesh::Mesh_MSTK(0.0,0.0,1.0,1.0,3,3,comm));

  nv = mesh->num_entities(Amanzi::AmanziMesh::NODE,Amanzi::AmanziMesh::Parallel_type::OWNED);
  CHECK_EQUAL(NV,nv);
  
  nf = mesh->num_entities(Amanzi::AmanziMesh::FACE,Amanzi::AmanziMesh::Parallel_type::OWNED);
  CHECK_EQUAL(NF,nf);
  
  nc = mesh->num_entities(Amanzi::AmanziMesh::CELL,Amanzi::AmanziMesh::Parallel_type::OWNED);
  CHECK_EQUAL(NC,nc);


  std::vector<Amanzi::AmanziMesh::Entity_ID>  c2f(6);
  Epetra_Map cell_map(mesh->cell_map(false));
  Epetra_Map face_map(mesh->face_map(false));
  for (int c=cell_map.MinLID(); c<=cell_map.MaxLID(); c++)
    {
      CHECK_EQUAL(cell_map.GID(c),mesh->GID(c,Amanzi::AmanziMesh::CELL));
      mesh->cell_get_faces(c, &c2f, true);
      for (j=0; j<4; j++)
	{
	  int f = face_map.LID(c2f[j]);
	  CHECK( f == c2f[j] );
	}

    }

  std::ofstream fout("test/mstk_quad_gen_4x4.out");
  Amanzi::MeshAudit auditor(mesh,fout);
  auditor.Verify();
}

