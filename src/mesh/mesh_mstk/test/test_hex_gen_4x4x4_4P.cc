#include <UnitTest++.h>

#include <iostream>

#include "../Mesh_MSTK.hh"

#include "Epetra_Map.h"
#include "Epetra_MpiComm.h"

// #include "MeshAudit.hh"

#include "mpi.h"

// Test for generation of hex mesh distributed over 4 processors

TEST(MSTK_HEX_GEN_4x4x4_4P)
{

  int i, j, k, err, nc, nf, nv;
  std::vector<Amanzi::AmanziMesh::Entity_ID> faces(6), nodes(8);
  std::vector<int> facedirs(6);
  std::vector<Amanzi::AmanziGeometry::Point> ccoords(8), fcoords(4);

  int NVowned[4] = {24,20,16,4}; // updated
  int NFowned[4] = {29,31,29,19}; // updated
  int NCowned[4] = {6,7,7,7}; // updated
  int NVused[4] = {48,48,60,64}; // updated
  int NFused[4] = {75,75,98,108}; // updated
  int NCused[4] = {18,18,24,27}; // updated
  int NVghost[4] = {24,28,44,60}; // updated
  int NFghost[4] = {46,44,69,89}; // updated
  int NCghost[4] = {12,11,17,20}; // updated

			      
  Teuchos::RCP<Epetra_MpiComm> comm(new Epetra_MpiComm(MPI_COMM_WORLD));
			      

  int rank, size;

  int initialized;

  MPI_Initialized(&initialized);
  
  if (!initialized)
    MPI_Init(NULL,NULL);

  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  CHECK_EQUAL(4,size);

  if (size != 4) {
    cerr << "Test must be run with 4 processors" << std::endl;
    //    return;
  }

  // if (rank == 0) {
  int DebugWait = 0;
  while (DebugWait);
  // }

  // Create a 3x3x3 cell hex mesh

  Amanzi::AmanziMesh::Mesh_MSTK mesh(0.0,0.0,0.0,1.0,1.0,1.0,3,3,3,comm.get());

  nv = mesh.num_entities(Amanzi::AmanziMesh::NODE,Amanzi::AmanziMesh::OWNED);  
  CHECK_EQUAL(NVowned[rank],nv);
  
  nf = mesh.num_entities(Amanzi::AmanziMesh::FACE,Amanzi::AmanziMesh::OWNED);  
  CHECK_EQUAL(NFowned[rank],nf);
  
  nc = mesh.num_entities(Amanzi::AmanziMesh::CELL,Amanzi::AmanziMesh::OWNED);
  CHECK_EQUAL(NCowned[rank],nc);

  nv = mesh.num_entities(Amanzi::AmanziMesh::NODE,Amanzi::AmanziMesh::USED);  
  CHECK_EQUAL(NVused[rank],nv);
  
  nf = mesh.num_entities(Amanzi::AmanziMesh::FACE,Amanzi::AmanziMesh::USED);  
  CHECK_EQUAL(NFused[rank],nf);
  
  nc = mesh.num_entities(Amanzi::AmanziMesh::CELL,Amanzi::AmanziMesh::USED);
  CHECK_EQUAL(NCused[rank],nc);

  nv = mesh.num_entities(Amanzi::AmanziMesh::NODE,Amanzi::AmanziMesh::GHOST);  
  CHECK_EQUAL(NVghost[rank],nv);
  
  nf = mesh.num_entities(Amanzi::AmanziMesh::FACE,Amanzi::AmanziMesh::GHOST);  
  CHECK_EQUAL(NFghost[rank],nf);
  
  nc = mesh.num_entities(Amanzi::AmanziMesh::CELL,Amanzi::AmanziMesh::GHOST);
  CHECK_EQUAL(NCghost[rank],nc);


  std::vector<Amanzi::AmanziMesh::Entity_ID>  c2f(6);
  std::vector<int> c2fdirs(6);
  Epetra_Map cell_map(mesh.cell_epetra_map(false));
  Epetra_Map face_map(mesh.face_epetra_map(true));

  for (int c=cell_map.MinLID(); c<=cell_map.MaxLID(); c++)
    {
      CHECK_EQUAL(cell_map.GID(c),mesh.GID(c,Amanzi::AmanziMesh::CELL));
      mesh.cell_get_faces(c, &c2f);
      mesh.cell_get_face_dirs(c, &c2fdirs);

      for (int j=0; j<6; j++)
	{
	  int f = face_map.LID(mesh.GID(c2f[j],Amanzi::AmanziMesh::FACE));
	  CHECK_EQUAL( f,c2f[j] );
	  if (f != c2f[j]) {
	    cout << std::endl;
	    cout << "Processor ID " << rank << std::endl;
	    cout << "Cell ID " << cell_map.GID(c) << std::endl;
	    cout << "Problem face c2f[j] = " << c2f[j] << " GID = " << mesh.GID(c2f[j],Amanzi::AmanziMesh::FACE) << " f = " << f << std::endl;
	    cout << std::endl;
	  }
	}

    }

  //  std::stringstream fname;
  //  fname << "mstk_hex_gen_4x4x4_4P." << rank << ".out";
  //  std::ofstream fout(fname.str().c_str());
  //  Amanzi::MeshAudit auditor(mesh,fname);
  //  auditor.verify();  

}

