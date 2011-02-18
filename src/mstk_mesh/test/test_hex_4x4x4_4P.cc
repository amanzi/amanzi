#include <UnitTest++.h>

#include <iostream>

#include "../Mesh_maps_mstk.hh"
#include "../../mesh_data/Entity_kind.hh"


#include "Epetra_Map.h"
#include "Epetra_MpiComm.h"

#include "mpi.h"


TEST(MSTK_HEX_4x4x4_4P)
{

  int i, j, k, err, nc, nf, nv;
  unsigned int faces[6], nodes[8];
  int facedirs[6];
  double ccoords[24], fcoords[12];

  int NVowned[4] = {24,20,16,4}; // updated
  int NFowned[4] = {29,31,29,19}; // updated
  int NCowned[4] = {6,7,7,7}; // updated
  int NVused[4] = {48,48,60,64}; // updated
  int NFused[4] = {75,75,98,108}; // updated
  int NCused[4] = {18,18,24,27}; // updated
  int NVghost[4] = {24,28,44,60}; // updated
  int NFghost[4] = {46,44,69,89}; // updated
  int NCghost[4] = {12,11,17,20}; // updated

			      
			      

  int rank, size;

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

  Mesh_maps_mstk mesh("test/hex_4x4x4_ss.exo",MPI_COMM_WORLD);


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
	  if (f != c2f[j]) {
	    cout << std::endl;
	    cout << "Processor ID " << rank << std::endl;
	    cout << "Cell ID " << cell_map.GID(c) << std::endl;
	    cout << "Problem face c2f[j] = " << c2f[j] << " GID = " << mesh.GID(c2f[j],Mesh_data::FACE) << " f = " << f << std::endl;
	    cout << std::endl;
	  }
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

