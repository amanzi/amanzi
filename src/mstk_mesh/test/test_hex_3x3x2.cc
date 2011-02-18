#include <UnitTest++.h>

#include <iostream>

#include "../Mesh_maps_mstk.hh"
#include "../../mesh_data/Entity_kind.hh"


#include "Epetra_Map.h"
#include "Epetra_MpiComm.h"

#include "mpi.h"


TEST(MSTK_HEX_3x3x2)
{

  int i, j, k, err, nc, nf, nv;
  unsigned int faces[6], cnodes[8], fnodes[6], expfacenodes[4];
  int facedirs[6];
  double ccoords[24], fcoords[12];

  int NV = 18;
  int NF = 20;
  int NC = 4;
  double xyz[18][3] = {{-0.5,-0.5, 0.25},
		       {-0.5,-0.5,-0.25},
		       {-0.5, 0,  -0.25},
		       {-0.5, 0,   0.25},
		       { 0,  -0.5, 0.25},
		       { 0,  -0.5,-0.25},
		       { 0,   0,  -0.25},
		       { 0,   0,   0.25},
		       {-0.5, 0.5,-0.25},
		       {-0.5, 0.5, 0.25},
		       { 0,   0.5,-0.25},
		       { 0,   0.5, 0.25},
		       { 0.5,-0.5, 0.25},
		       { 0.5,-0.5,-0.25},
		       { 0.5, 0,  -0.25},
		       { 0.5, 0,   0.25},
		       { 0.5, 0.5,-0.25},
		       { 0.5, 0.5, 0.25}};

  unsigned int cnstd[8] = {0,1,2,3,4,5,6,7};
  unsigned int cfstd[6][4] = {{0,1,5,4},    // Expected cell-face-node pattern
			      {1,2,6,5},
			      {2,3,7,6},
			      {3,0,4,7},
			      {0,3,2,1},
			      {4,5,6,7}};


  // Load a simple 2 element hex mesh 

  Mesh_maps_mstk mesh("test/hex_3x3x2_ss.exo",MPI_COMM_WORLD);


  nv = mesh.count_entities(Mesh_data::NODE,OWNED);
  CHECK_EQUAL(NV,nv);

  for (i = 0; i < nv; i++) {
    double coords[3];

    mesh.node_to_coordinates(i,coords,coords+6);
    CHECK_ARRAY_EQUAL(xyz[i],coords,3);
  }


  nf = mesh.count_entities(Mesh_data::FACE,OWNED);
  CHECK_EQUAL(NF,nf);
  
  nc = mesh.count_entities(Mesh_data::CELL,OWNED);
  CHECK_EQUAL(NC,nc);
    
  for (i = 0; i < nc; i++) {
    mesh.cell_to_nodes(i,cnodes,cnodes+8);
    mesh.cell_to_faces(i,faces,faces+6);
    mesh.cell_to_face_dirs(i,facedirs,facedirs+6);

    for (j = 0; j < 6; j++) {

      mesh.face_to_nodes(faces[j],fnodes,fnodes+4);
      mesh.face_to_coordinates(faces[j],fcoords,fcoords+12);
      

      for (k = 0; k < 4; k++)
	expfacenodes[k] = cnodes[cfstd[j][k]];

      // The order of nodes returned may be different from what we expected
      // So make sure we have a matching node to start with
      
      int k0 = -1;
      int found = 0;
      for (k = 0; k < 4; k++) {
	if (expfacenodes[k] == fnodes[0]) {
	  k0 = k;
	  found = 1;
	  break;
	}
      }
      
      CHECK_EQUAL(found,1); 
      
      if (facedirs[j] == 1) {
	for (k = 0; k < 4; k++) {
	  CHECK_EQUAL(expfacenodes[(k0+k)%4],fnodes[k]);
	  CHECK_ARRAY_EQUAL(xyz[expfacenodes[(k0+k)%4]],&(fcoords[3*k]),3);
	}
      }
      else {
	for (k = 0; k < 4; k++) {
	  CHECK_EQUAL(expfacenodes[(k0+4-k)%4],fnodes[k]);
	  CHECK_ARRAY_EQUAL(xyz[expfacenodes[(k0+4-k)%4]],&(fcoords[3*k]),3);
	}
      }
      
    }

  
    mesh.cell_to_coordinates(i,ccoords,ccoords+24);
    
    for (j = 0; j < 8; j++) {
      CHECK_ARRAY_EQUAL(xyz[cnodes[cnstd[j]]],&(ccoords[3*j]),3);
    }
  }


  // Verify the sidesets

  int ns;
  ns = mesh.num_sets(Mesh_data::FACE);
  CHECK_EQUAL(7,ns);

  unsigned int setids[7], expsetids[7]={1,101,102,103,104,105,106};
  mesh.get_set_ids(Mesh_data::FACE,setids,setids+7);
  
  CHECK_ARRAY_EQUAL(expsetids,setids,7);

  unsigned int setsize, expsetsizes[7] = {16,4,4,2,2,2,2};
  unsigned int expsetfaces[7][16] = {{3,8,14,18,1,6,12,16,0,11,4,9,7,17,15,19},
				     {3,8,14,18,0,0,0,0,0,0,0,0,0,0,0,0},
				     {1,6,12,16,0,0,0,0,0,0,0,0,0,0,0,0},
				     {0,11,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
				     {4,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
				     {7,17,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
				     {15,19,0,0,0,0,0,0,0,0,0,0,0,0,0,0}};


  for (i = 0; i < ns; i++) {
    unsigned int setfaces[16];

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
	cerr << "Could not find FACE " << expsetfaces[i][j] << " in returned set\n";
	CHECK_EQUAL(found,1);
	break;
      }
    }
    
  }



  std::vector<unsigned int>  c2f(6);
  Epetra_Map cell_map(mesh.cell_map(true));
  Epetra_Map face_map(mesh.face_map(false));

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

}

