#include <UnitTest++.h>

#include <iostream>

#include "../Mesh_MOAB.hh"


#include "Epetra_Map.h"
#include "Epetra_MpiComm.h"

#include "mpi.h"

// Unless this example is enhanced, it does lesser testing than test_hex_3x3x2.cc

TEST(MOAB_HEX_4x4x4)
{

  int i, j, k, err, nc, nf, nv;
  unsigned int faces[6], nodes[8];
  int facedirs[6];
  double ccoords[24], fcoords[12];

  int NV = 64;
  int NF = 108;
  int NC = 27;


  std::auto_ptr<Epetra_MpiComm> comm(new Epetra_MpiComm(MPI_COMM_WORLD));

  // Load a single hex from the hex1.exo file

  Amanzi::AmanziMesh::Mesh_MOAB mesh("test/hex_4x4x4_ss.exo",comm.get());

  nf = mesh.count_entities(Amanzi::AmanziMesh::FACE,Amanzi::AmanziMesh::OWNED);
  CHECK_EQUAL(NF,nf);
  
  nc = mesh.count_entities(Amanzi::AmanziMesh::CELL,Amanzi::AmanziMesh::OWNED);
  CHECK_EQUAL(NC,nc);


  std::vector<unsigned int>  c2f(6);
  std::vector<int> c2fdirs(6);
  Epetra_Map cell_map(mesh.cell_epetra_map(false));
  Epetra_Map face_map(mesh.face_epetra_map(false));
  for (int c=cell_map.MinLID(); c<=cell_map.MaxLID(); c++)
    {
      mesh.cell_get_faces_and_dirs( c, &c2f, &c2fdirs, true);
      for (int j=0; j<6; j++)
	{
	  int f = face_map.LID(c2f[j]);
	  CHECK( f == c2f[j] );
	}

    }
  
  int ns;

  // Verify cell sets

  ns = mesh.num_sets(Amanzi::AmanziMesh::CELL);
  CHECK_EQUAL(3,ns);

  std::vector<unsigned int> csetids(3);
  unsigned int expcsetids[3] = {10000,20000,30000};

  mesh.get_set_ids(Amanzi::AmanziMesh::CELL,&csetids);

  CHECK_ARRAY_EQUAL(expcsetids,csetids,3);

  
  unsigned int csetsize, expcsetsizes[3] = {9,9,9};

  unsigned int expcsetcells[3][9] = {{0,1,2,3,4,5,6,7,8},
				     {9,10,11,12,13,14,15,16,17},
				     {18,19,20,21,22,23,24,25,26}};

  for (i = 0; i < ns; i++) {
    unsigned int setcells[9];

    csetsize = mesh.get_set_size(csetids[i],Amanzi::AmanziMesh::CELL,Amanzi::AmanziMesh::OWNED);
    CHECK_EQUAL(expcsetsizes[i],csetsize);


    mesh.get_set(csetids[i],Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::OWNED, setcells, setcells+csetsize);
    
    CHECK_ARRAY_EQUAL(expcsetcells[i],setcells,csetsize);
  }


  // Verify the sidesets

  ns = mesh.num_sets(Amanzi::AmanziMesh::FACE);
  CHECK_EQUAL(21,ns);
  
  std::vector<unsigned int> fsetids(21);
  unsigned int  expfsetids[21]={1,101,102,103,104,105,106,
				10001,10002,20001,30001,10003,
				20002,30002,10004,20003,30003,
				10005,20004,30004,30005};
  mesh.get_set_ids(Amanzi::AmanziMesh::FACE,&fsetids);
  
  CHECK_ARRAY_EQUAL(expfsetids,fsetids,21);

  // We won't check the super-sideset, sideset number 1 with 48 faces
  unsigned int fsetsize, expfsetsizes[20] = {9,9,9,9,9,9,
				    9,3,3,3,3,
                                    3,3,3,3,
                                    3,3,3,3,9};

  unsigned int expfsetfaces[20][9] = {{0,1,2,3,4,5,6,7,8},
				      {9,10,11,12,13,14,15,16,17},
				      {18,19,20,21,22,23,24,25,26},
				      {27,28,29,30,31,32,33,34,35},
				      {36,37,38,39,40,41,42,43,44},
				      {45,46,47,48,49,50,51,52,53},
				      {0,1,2,3,4,5,6,7,8},
				      {9,10,12,0,0,0,0,0,0},
				      {11,13,15,0,0,0,0,0,0},
				      {14,16,17,0,0,0,0,0,0},
				      {18,20,23,0,0,0,0,0,0},
				      {19,22,25,0,0,0,0,0,0},
				      {21,24,26,0,0,0,0,0,0},
				      {27,28,30,0,0,0,0,0,0},
				      {29,31,33,0,0,0,0,0,0},
				      {32,34,35,0,0,0,0,0,0},
				      {36,38,41,0,0,0,0,0,0},
				      {37,40,43,0,0,0,0,0,0},
				      {39,42,44,0,0,0,0,0,0},
				      {45,46,47,48,49,50,51,52,53}};
			   


  for (i = 0; i < ns-1; i++) {
    unsigned int setfaces[9];

    CHECK_EQUAL(true,mesh.valid_set_id(fsetids[i+1],Amanzi::AmanziMesh::FACE));
   
    fsetsize = mesh.get_set_size(fsetids[i+1],Amanzi::AmanziMesh::FACE,Amanzi::AmanziMesh::OWNED);
    CHECK_EQUAL(expfsetsizes[i],fsetsize);


    mesh.get_set(fsetids[i+1],Amanzi::AmanziMesh::FACE, Amanzi::AmanziMesh::OWNED, setfaces, setfaces+fsetsize);
    
    CHECK_ARRAY_EQUAL(expfsetfaces[i],setfaces,fsetsize);
  }


}

