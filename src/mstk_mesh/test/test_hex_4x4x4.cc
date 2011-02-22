#include <UnitTest++.h>

#include <iostream>

#include "../Mesh_maps_mstk.hh"
#include "../../mesh_data/Entity_kind.hh"


#include "Epetra_Map.h"
#include "Epetra_MpiComm.h"

#include "mpi.h"

// Unless this example is enhanced, it does lesser testing than test_hex_3x3x2.cc

TEST(MSTK_HEX_4x4x4)
{

  int i, j, k, err, nc, nf, nv;
  unsigned int faces[6], nodes[8];
  int facedirs[6];
  double ccoords[24], fcoords[12];

  int NV = 64;
  int NF = 108;
  int NC = 27;


  // Load a mesh consisting of 3x3x3 elements (4x4x4 nodes)

  Mesh_maps_mstk mesh("test/hex_4x4x4_ss.exo",MPI_COMM_WORLD);

  nf = mesh.count_entities(Mesh_data::FACE,OWNED);
  CHECK_EQUAL(NF,nf);
  
  nc = mesh.count_entities(Mesh_data::CELL,OWNED);
  CHECK_EQUAL(NC,nc);


  std::vector<unsigned int>  c2f(6);
  std::vector<int> c2fdirs(6);
  Epetra_Map cell_map(mesh.cell_map(false));
  Epetra_Map face_map(mesh.face_map(false));
  for (int c=cell_map.MinLID(); c<=cell_map.MaxLID(); c++)
    {
      mesh.cell_to_faces( c, c2f.begin(), c2f.end() );
      mesh.cell_to_face_dirs(c, c2fdirs.begin(), c2fdirs.end() );
      for (int j=0; j<6; j++)
	{
	  int f = face_map.LID(c2f[j]);
	  CHECK( f == c2f[j] );
	}

    }
  
  int ns;

  // Verify cell sets

  ns = mesh.num_sets(Mesh_data::CELL);
  CHECK_EQUAL(3,ns);

  unsigned int csetids[3], expcsetids[3] = {10000,20000,30000};

  mesh.get_set_ids(Mesh_data::CELL,csetids,csetids+3);

  CHECK_ARRAY_EQUAL(expcsetids,csetids,3);

  
  unsigned int csetsize, expcsetsizes[3] = {9,9,9};

  unsigned int expcsetcells[3][9] = {{0,1,2,3,4,5,6,7,8},
				     {9,10,11,12,13,14,15,16,17},
				     {18,19,20,21,22,23,24,25,26}};

  for (i = 0; i < ns; i++) {
    unsigned int setcells[9];

    csetsize = mesh.get_set_size(csetids[i],Mesh_data::CELL,OWNED);
    CHECK_EQUAL(expcsetsizes[i],csetsize);


    mesh.get_set(csetids[i],Mesh_data::CELL, OWNED, setcells, setcells+csetsize);
    
    CHECK_ARRAY_EQUAL(expcsetcells[i],setcells,csetsize);
  }


  // Verify the sidesets

  ns = mesh.num_sets(Mesh_data::FACE);
  CHECK_EQUAL(21,ns);
  
  unsigned int fsetids[21], expfsetids[21]={1,101,102,103,104,105,106,
					    10001,10002,20001,30001,10003,
					    20002,30002,10004,20003,30003,
					    10005,20004,30004,30005};
  mesh.get_set_ids(Mesh_data::FACE,fsetids,fsetids+21);
  
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
    unsigned int setfaces[9];

    CHECK_EQUAL(true,mesh.valid_set_id(fsetids[i+1],Mesh_data::FACE));
   
    fsetsize = mesh.get_set_size(fsetids[i+1],Mesh_data::FACE,OWNED);
    CHECK_EQUAL(expfsetsizes[i],fsetsize);


    mesh.get_set(fsetids[i+1],Mesh_data::FACE, OWNED, setfaces, setfaces+fsetsize);
    
    CHECK_ARRAY_EQUAL(expfsetfaces[i],setfaces,fsetsize);
  }


}

