#include <UnitTest++.h>

#include <iostream>

#include "../Mesh_MSTK.hh"

#include "Epetra_Map.h"
#include "Epetra_MpiComm.h"

#include "MeshAudit.hh"

#include "mpi.h"


TEST(MSTK_QUAD_GEN_4x4x4_4P)
{

  int i, j, k, err, nc, nf, nv;
  std::vector<Amanzi::AmanziMesh::Entity_ID> faces(6), nodes(8);
  std::vector<int> facedirs(4);
  std::vector<Amanzi::AmanziGeometry::Point> ccoords(8), fcoords(4);

  int NVowned[4] = {6,4,2,4}; 
  int NFowned[4] = {7,6,4,7}; 
  int NCowned[4] = {2,2,2,3}; 
  int NVused[4] = {12,12,16,12}; 
  int NFused[4] = {17,17,24,17}; 
  int NCused[4] = {6,6,9,6}; 
  int NVghost[4] = {6,8,14,8}; 
  int NFghost[4] = {10,11,20,10}; 
  int NCghost[4] = {4,4,7,3}; 

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

  //  if (rank == 0) {
  int DebugWait = 0;
  while (DebugWait);
    //  }


  // Load a single hex from the hex1.exo file

  Amanzi::AmanziMesh::Mesh_MSTK mesh(0.0,0.0,1.0,1.0,3,3,comm.get());


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


  std::vector<Amanzi::AmanziMesh::Entity_ID>  c2f(4);
  std::vector<int> c2fdirs(4);
  Epetra_Map cell_map(mesh.cell_epetra_map(false));
  Epetra_Map face_map(mesh.face_epetra_map(true));

  for (int c=cell_map.MinLID(); c<=cell_map.MaxLID(); c++)
    {
      CHECK_EQUAL(cell_map.GID(c),mesh.GID(c,Amanzi::AmanziMesh::CELL));
      mesh.cell_get_faces_and_dirs(c, true, &c2f, &c2fdirs);

      for (int j=0; j<4; j++)
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
  //  fname << "mstk_quad_gen_4x4_4P." << rank << ".out";
  //  std::ofstream fout(fname.str().c_str());
  //  Amanzi::MeshAudit auditor(mesh,fout);
  //  auditor.Verify();


}

