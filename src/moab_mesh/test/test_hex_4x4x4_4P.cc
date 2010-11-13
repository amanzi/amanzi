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
  Epetra_Map cell_map(mesh.cell_map(false));
  Epetra_Map face_map(mesh.face_map(true));

  for (int c=cell_map.MinLID(); c<=cell_map.MaxLID(); c++)
    {
      CHECK_EQUAL(cell_map.GID(c),mesh.GID(c,Mesh_data::CELL));
      mesh.cell_to_faces( c, c2f.begin(), c2f.end() );
      for (int j=0; j<6; j++)
	{
	  int f = face_map.LID(mesh.GID(c2f[j],Mesh_data::FACE));
	  CHECK( f == c2f[j] );
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

