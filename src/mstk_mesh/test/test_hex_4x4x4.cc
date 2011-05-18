#include <UnitTest++.h>

#include <iostream>

#include "../Mesh_MSTK.hh"


#include "Epetra_Map.h"
#include "Epetra_MpiComm.h"

#include "mpi.h"

// Unless this example is enhanced, it does lesser testing than test_hex_3x3x2.cc

TEST(MSTK_HEX_4x4x4)
{

  int i, j, k, err, nc, nf, nv;
  unsigned int faces[6], nodes[8];
  int facedirs[6];

  int NV = 64;
  int NF = 108;
  int NC = 27;


  // Load a mesh consisting of 3x3x3 elements (4x4x4 nodes)

  Amanzi::AmanziMesh::Mesh_MSTK mesh("test/hex_4x4x4_ss.exo",MPI_COMM_WORLD,3);

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

  unsigned int expcsetids[3] = {10000,20000,30000};
  std::vector<Amanzi::AmanziMesh::Set_ID> csetids(3);

  mesh.get_set_ids(Amanzi::AmanziMesh::CELL,&csetids);

  CHECK_ARRAY_EQUAL(expcsetids,csetids,3);

  
  unsigned int csetsize, expcsetsizes[3] = {9,9,9};

  unsigned int expcsetcells[3][9] = {{0,1,2,3,4,5,6,7,8},
				     {9,10,11,12,13,14,15,16,17},
				     {18,19,20,21,22,23,24,25,26}};

  for (i = 0; i < ns; i++) {
    std::vector<Amanzi::AmanziMesh::Entity_ID> setcells(9);

    csetsize = mesh.get_set_size(csetids[i],Amanzi::AmanziMesh::CELL,Amanzi::AmanziMesh::OWNED);
    CHECK_EQUAL(expcsetsizes[i],csetsize);


    mesh.get_set_entities(csetids[i],Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::OWNED, &setcells);
    
    CHECK_ARRAY_EQUAL(expcsetcells[i],setcells,csetsize);
  }


  // Verify the sidesets

  ns = mesh.num_sets(Amanzi::AmanziMesh::FACE);
  CHECK_EQUAL(21,ns);
  
  unsigned int expfsetids[21]={1,101,102,103,104,105,106,
                               10001,10002,20001,30001,10003,
			       20002,30002,10004,20003,30003,
			       10005,20004,30004,30005};
  std::vector<Amanzi::AmanziMesh::Set_ID> fsetids(21);

  mesh.get_set_ids(Amanzi::AmanziMesh::FACE,&fsetids);
  
  CHECK_ARRAY_EQUAL(expfsetids,fsetids,21);

  // We won't check the super-sideset, sideset number 1 with 48 faces
  unsigned int fsetsize, expfsetsizes[20] = {9,9,9,9,9,9,
				    9,3,3,3,3,
                                    3,3,3,3,
                                    3,3,3,3,9};

  unsigned int expfsetfaces[20][9] = {{4,19,9,32,23,14,36,27,40},
				      {0,6,42,11,47,75,51,80,84},
				      {3,45,18,78,57,31,90,67,100},
				      {12,25,52,38,62,85,72,95,105},
				      {30,66,35,99,70,39,103,73,106},
				      {79,83,91,87,94,101,97,104,107},
				      {4,19,9,32,23,14,36,27,40},
				      {0,6,11,0,0,0,0,0,0},
				      {42,47,51,0,0,0,0,0,0},
				      {75,80,84,0,0,0,0,0,0},
				      {3,18,31,0,0,0,0,0,0},
				      {45,57,67,0,0,0,0,0,0},
				      {78,90,100,0,0,0,0,0,0},
				      {12,25,38,0,0,0,0,0,0},
				      {52,62,72,0,0,0,0,0,0},
				      {85,95,105,0,0,0,0,0,0},
				      {30,35,39,0,0,0,0,0,0},
				      {66,70,73,0,0,0,0,0,0},
				      {99,103,106,0,0,0,0,0,0},
				      {79,83,91,87,94,101,97,104,107}};
			   


  for (i = 0; i < ns-1; i++) {
    std::vector<Amanzi::AmanziMesh::Entity_ID> setfaces(9);

    CHECK_EQUAL(true,mesh.valid_set_id(fsetids[i+1],Amanzi::AmanziMesh::FACE));
   
    fsetsize = mesh.get_set_size(fsetids[i+1],Amanzi::AmanziMesh::FACE,Amanzi::AmanziMesh::OWNED);
    CHECK_EQUAL(expfsetsizes[i],fsetsize);


    mesh.get_set_entities(fsetids[i+1],Amanzi::AmanziMesh::FACE, Amanzi::AmanziMesh::OWNED, &setfaces);
    
    CHECK_ARRAY_EQUAL(expfsetfaces[i],setfaces,fsetsize);
  }


}

