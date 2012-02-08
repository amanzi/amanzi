#include <UnitTest++.h>

#include <iostream>

#include "../Mesh_MSTK.hh"


#include "Epetra_Map.h"
#include "Epetra_MpiComm.h"

#include "MeshAudit.hh"

#include "mpi.h"

// Test generation of quad mesh in serial

TEST(MSTK_QUAD_GEN_4x4)
{

  int i, j, k, err, nc, nf, nv;
  unsigned int faces[6], nodes[8];
  int facedirs[6];

  int NV = 16;
  int NF = 24;
  int NC = 9;

  Teuchos::RCP<Epetra_MpiComm> comm(new Epetra_MpiComm(MPI_COMM_WORLD));


  // Load a mesh consisting of 3x3 elements (4x4 nodes)

  Amanzi::AmanziMesh::Mesh_MSTK mesh(0.0,0.0,1.0,1.0,3,3,comm.get());

  nv = mesh.num_entities(Amanzi::AmanziMesh::NODE,Amanzi::AmanziMesh::OWNED);
  CHECK_EQUAL(NV,nv);
  
  nf = mesh.num_entities(Amanzi::AmanziMesh::FACE,Amanzi::AmanziMesh::OWNED);
  CHECK_EQUAL(NF,nf);
  
  nc = mesh.num_entities(Amanzi::AmanziMesh::CELL,Amanzi::AmanziMesh::OWNED);
  CHECK_EQUAL(NC,nc);


  std::vector<unsigned int>  c2f(6);
  std::vector<int> c2fdirs(6);
  Epetra_Map cell_map(mesh.cell_epetra_map(false));
  Epetra_Map face_map(mesh.face_epetra_map(false));
  for (int c=cell_map.MinLID(); c<=cell_map.MaxLID(); c++)
    {
      CHECK_EQUAL(cell_map.GID(c),mesh.GID(c,Amanzi::AmanziMesh::CELL));
      mesh.cell_get_faces(c, &c2f);
      mesh.cell_get_face_dirs(c, &c2fdirs);
      for (int j=0; j<4; j++)
	{
	  int f = face_map.LID(c2f[j]);
	  CHECK( f == c2f[j] );
	}

    }

  //  std::ofstream fout("mstk_quad_gen_4x4.out");
  //  Amanzi::MeshAudit auditor(mesh,fout);
  //  auditor.Verify();

}

