#include <UnitTest++.h>

#include <iostream>


#include "../Mesh_maps_mstk.hh"
#include "../../mesh_data/Entity_kind.hh"

#include "Epetra_Map.h"
#include "Epetra_MpiComm.h"

#include "mpi.h"


TEST(MSTK_HEX1)
{

  int i, j, k, err, nc, nv;
  unsigned int faces[6], facenodes[4], cellnodes[8], expfacenodes[4];
  int facedirs[6];
  double ccoords[24], fcoords[12];

  int NV = 8;
  int NF = 6;
  int NC = 1;
  double xyz[12][3] = {{0, 0, 0},
		       {1, 0, 0},
		       {0, 1, 0},
		       {1, 1, 0},
		       {0, 0, 1},
		       {1, 0, 1}, 
		       {0, 1, 1},
		       {1, 1, 1}};
  unsigned int local_cellnodes[8] = {0,1,2,3,4,5,6,7};
  unsigned int local_facenodes[6][4] = {{0,1,5,4},
					{1,2,6,5},
					{2,3,7,6},
					{3,0,4,7},
					{0,3,2,1},
					{4,5,6,7}};


  // Load a single hex from the hex1.exo file

  Mesh_maps_mstk mesh("test/hex_2x2x2_ss.exo",MPI_COMM_WORLD);


  // Check number of nodes and their coordinates

  nv = mesh.count_entities(Mesh_data::NODE, OWNED);
  CHECK_EQUAL(NV,nv);

  for (i = 0; i < nv; i++) {
    double coords[3];
    
    mesh.node_to_coordinates(i,coords,coords+6);
    CHECK_ARRAY_EQUAL(xyz[i],coords,3);
  }


  // Check number of cells and their face nodes and their face coordinates
  
  nc = mesh.count_entities(Mesh_data::CELL,OWNED);
  CHECK_EQUAL(NC,nc);


  // Check cell coordinates directly

  mesh.cell_to_nodes(0,cellnodes, cellnodes+8);
  mesh.cell_to_coordinates(0,ccoords,ccoords+24);
    
  for (j = 0; j < 8; j++) {
    CHECK_ARRAY_EQUAL(xyz[cellnodes[j]],&(ccoords[3*j]),3);
  }


    
  mesh.cell_to_faces(0,faces,faces+6);
  mesh.cell_to_face_dirs(0,facedirs,facedirs+6);

  for (j = 0; j < 6; j++) {
    mesh.face_to_nodes(faces[j],facenodes,facenodes+4);
    mesh.face_to_coordinates(faces[j],fcoords,fcoords+12);


    for (k = 0; k < 4; k++)
      expfacenodes[k] = cellnodes[local_facenodes[j][k]];

    // The order of nodes returned may be different from what we expected
    // So make sure we have a matching node to start with

    int k0 = -1;
    int found = 0;
    for (k = 0; k < 4; k++) {
      if (expfacenodes[k] == facenodes[0]) {
	k0 = k;
	found = 1;
	break;
      }
    }

    CHECK_EQUAL(found,1); 

    if (facedirs[j] == 1) {
      for (k = 0; k < 4; k++) {
	CHECK_EQUAL(expfacenodes[(k0+k)%4],facenodes[k]);
	CHECK_ARRAY_EQUAL(xyz[expfacenodes[(k0+k)%4]],&(fcoords[3*k]),3);
      }
    }
    else {
      for (k = 0; k < 4; k++) {
	CHECK_EQUAL(expfacenodes[(k0+4-k)%4],facenodes[k]);
	CHECK_ARRAY_EQUAL(xyz[expfacenodes[(k0+4-k)%4]],&(fcoords[3*k]),3);
      }
    }

  }
      

  // Verify the sidesets

  int ns;
  ns = mesh.num_sets(Mesh_data::FACE);
  CHECK_EQUAL(7,ns);

  unsigned int setids[7], expsetids[7]={1,101,102,103,104,105,106};
  mesh.get_set_ids(Mesh_data::FACE,setids,setids+7);
  
  CHECK_ARRAY_EQUAL(expsetids,setids,7);

  unsigned int setsize, expsetsizes[7] = {6,1,1,1,1,1,1};
  unsigned int expsetfaces[7][6] = {{0,1,2,3,4,5},
				    {4,0,0,0,0,0},
				    {5,0,0,0,0,0},
				    {0,0,0,0,0,0},
				    {1,0,0,0,0,0},
				    {2,0,0,0,0,0},
				    {3,0,0,0,0,0}};


  for (i = 0; i < ns; i++) {
    unsigned int setfaces[6];

    setsize = mesh.get_set_size(setids[i],Mesh_data::FACE,OWNED);
    CHECK_EQUAL(expsetsizes[i],setsize);


    mesh.get_set(setids[i],Mesh_data::FACE, OWNED, setfaces, setfaces+setsize);
    
    int allfound = 1;
    for (j = 0; j < setsize; j++) {
      int found = 0;
      for (k = 0; k < setsize; k++) {
	if (expsetfaces[i][j] == setfaces[k]) {
	  found = 1;
	  break;
	}
      }
      if (!found) {
	allfound = 0;
	cerr << "Could not find FACE " << expsetfaces[j] << " in returned set\n";
	CHECK_EQUAL(found,1);
	break;
      }
    }
    
  }

  
}

