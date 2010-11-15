#include <UnitTest++.h>

#include <iostream>

#include "../Mesh_maps_moab.hh"
#include "../../mesh_data/Entity_kind.hh"


#include "Epetra_Map.h"
#include "Epetra_MpiComm.h"

#include "mpi.h"


TEST(MOAB_HEX_4x4x4_4P)
{

  int i, j, k, err, nc, nf, nv;
  unsigned int faces[6], nodes[8];
  int facedirs[6];
  double ccoords[24], fcoords[12];

  int NVowned[4] = {16,16,16,16};
  int NFowned[4] = {16,26,26,40};
  int NCowned[4] = {3,6,6,12};
  int NVused[4] = {16,24,24,36};
  int NFused[4] = {16,29,29,52};
  int NCused[4] = {3,6,6,12};
  int NVghost[4] = {20,32,32,48};
  int NFghost[4] = {36,49,49,68};
  int NCghost[4] = {9,12,12,15};

  int cfdirs[4][12][6] = {{{1,1,1,1,1,1},
			   {1,1,1,1,-1,1},
			   {1,1,1,1,-1,1},
			   {0,0,0,0,0,0},
			   {0,0,0,0,0,0},
			   {0,0,0,0,0,0},
			   {0,0,0,0,0,0},
			   {0,0,0,0,0,0},
			   {0,0,0,0,0,0},
			   {0,0,0,0,0,0},
			   {0,0,0,0,0,0},
			   {0,0,0,0,0,0},},
			  
			  {{1,1,1,1,1,1},
			   {-1,1,1,1,1,1},
			   {1,1,1,1,-1,1},
			   {-1,1,1,1,-1,1},
			   {1,1,1,1,-1,1},
			   {-1,1,1,1,-1,1},
			   {0,0,0,0,0,0},
			   {0,0,0,0,0,0},
			   {0,0,0,0,0,0},
			   {0,0,0,0,0,0},
			   {0,0,0,0,0,0},
			   {0,0,0,0,0,0}},
			  
			  {{1,1,1, 1, 1,1},
			   {1,1,1,-1, 1,1},
			   {1,1,1, 1,-1,1},
			   {1,1,1,-1,-1,1},
			   {1,1,1, 1,-1,1},
			   {1,1,1,-1,-1,1},
			   {0,0,0,0,0,0},
			   {0,0,0,0,0,0},
			   {0,0,0,0,0,0},
			   {0,0,0,0,0,0},
			   {0,0,0,0,0,0},
			   {0,0,0,0,0,0}},
			  
			  
			  {{ 1,1,1, 1, 1,1},
			   { 1,1,1,-1, 1,1},
			   {-1,1,1, 1, 1,1},
			   {-1,1,1,-1, 1,1},
			   { 1,1,1, 1,-1,1},
			   { 1,1,1,-1,-1,1},
			   {-1,1,1, 1,-1,1},
			   {-1,1,1,-1,-1,1},
			   { 1,1,1, 1,-1,1},
			   { 1,1,1,-1,-1,1},
			   {-1,1,1, 1,-1,1},
			   {-1,1,1,-1,-1,1}}
  };

  int cfaces[4][27][6] = {{{8,1,2,11,0,7},
			   {9,3,4,12,7,15},
			   {10,5,6,13,15,14},
			   {0,0,0,0,0,0},
			   {0,0,0,0,0,0},
			   {0,0,0,0,0,0},
			   {0,0,0,0,0,0},
			   {0,0,0,0,0,0},
			   {0,0,0,0,0,0},
			   {0,0,0,0,0,0},
			   {0,0,0,0,0,0},
			   {0,0,0,0,0,0}},
 
			  {{26,1,14,8,0,15},
			   {14,3,16,10,2,19},
			   {27,4,22,9,15,23},
			   {22,5,17,12,19,24},
			   {28,6,25,11,23,20},
			   {25,7,18,13,24,21},
			   {0,0,0,0,0,0},
			   {0,0,0,0,0,0},
			   {0,0,0,0,0,0},
			   {0,0,0,0,0,0},
			   {0,0,0,0,0,0},
			   {0,0,0,0,0,0}},
			  
			  {{7,14,0,26,1,15},
			   {9,17,2,14,3,16},
			   {10,22,4,27,15,23},
			   {11,18,5,22,16,24},
			   {12,25,6,28,23,20},
			   {13,19,8,25,24,21},
			   {0,0,0,0,0,0},
			   {0,0,0,0,0,0},
			   {0,0,0,0,0,0},
			   {0,0,0,0,0,0},
			   {0,0,0,0,0,0},
			   {0,0,0,0,0,0}},
			  
			  
			  {{43,11,17,42,0,18},
			   {41,4,19,11,2,20},
			   {17,25,10,40,1,26},
			   {19,5,13,25,3,27},
			   {44,28,29,45,18,30},
			   {46,6,31,28,20,32},
			   {29,33,12,47,26,34},
			   {31,7,15,33,27,35},
			   {48,36,37,49,30,21},
			   {50,8,38,36,32,22},
			   {37,39,14,51,34,23},
			   {38,9,16,39,35,24}}};

			      
			      

  int rank, size;

  MPI_Init(NULL,NULL);

  // Load a single hex from the hex1.exo file

  Mesh_maps_moab mesh("test/hex_4x4x4_ss_4P.h5m",MPI_COMM_WORLD);

  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  CHECK_EQUAL(4,size);

  nv = mesh.count_entities(Mesh_data::NODE,OWNED);  
  CHECK_EQUAL(NVowned[rank],nv);
  
  nf = mesh.count_entities(Mesh_data::FACE,OWNED);  
  CHECK_EQUAL(NFowned[rank],nf);
  
  nc = mesh.count_entities(Mesh_data::CELL,OWNED);
  CHECK_EQUAL(NCowned[rank],nc);

  nv = mesh.count_entities(Mesh_data::NODE,USED);  
  CHECK_EQUAL(NVused[rank],nv);
  
  nf = mesh.count_entities(Mesh_data::FACE,USED);  
  CHECK_EQUAL(NFused[rank],nf);
  
  nc = mesh.count_entities(Mesh_data::CELL,USED);
  CHECK_EQUAL(NCused[rank],nc);

  nv = mesh.count_entities(Mesh_data::NODE,GHOST);  
  CHECK_EQUAL(NVghost[rank],nv);
  
  nf = mesh.count_entities(Mesh_data::FACE,GHOST);  
  CHECK_EQUAL(NFghost[rank],nf);
  
  nc = mesh.count_entities(Mesh_data::CELL,GHOST);
  CHECK_EQUAL(NCghost[rank],nc);


  std::vector<unsigned int>  c2f(6);
  std::vector<int> c2fdirs(6);
  Epetra_Map cell_map(mesh.cell_map(false));
  Epetra_Map face_map(mesh.face_map(true));

  for (int c=cell_map.MinLID(); c<=cell_map.MaxLID(); c++)
    {
      CHECK_EQUAL(cell_map.GID(c),mesh.GID(c,Mesh_data::CELL));
      mesh.cell_to_faces( c, c2f.begin(), c2f.end() );
      mesh.cell_to_face_dirs(c, c2fdirs.begin(), c2fdirs.end());

      for (int j=0; j<6; j++)
	{
	  int f = face_map.LID(mesh.GID(c2f[j],Mesh_data::FACE));
	  CHECK_EQUAL( f,c2f[j] );
	  CHECK_EQUAL(cfaces[rank][c][j],f);
	  CHECK_EQUAL(cfdirs[rank][c][j],c2fdirs[j]);
	}

    }
  

  // Verify cell sets

  int ns;
  ns = mesh.num_sets(Mesh_data::CELL);
  CHECK_EQUAL(3,ns);

  unsigned int csetids[3], expcsetids[3] = {10000,20000,30000};

  mesh.get_set_ids(Mesh_data::CELL,csetids,csetids+3);

  CHECK_ARRAY_EQUAL(expcsetids,csetids,3);

  MPI_Finalize();

}

