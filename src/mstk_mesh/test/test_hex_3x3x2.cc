#include <UnitTest++.h>

#include <iostream>

#include "../Mesh_MSTK.hh"


#include "Epetra_Map.h"
#include "Epetra_MpiComm.h"

#include "mpi.h"

using namespace Amanzi;
using namespace AmanziMesh;
using namespace AmanziGeometry;


TEST(MSTK_HEX_3x3x2)
{

  int i, j, k, err, nc, nf, nv;
  std::vector<Entity_ID> faces(6), cnodes(8), fnodes(6), expfacenodes(4);
  std::vector<int> facedirs(6);
  std::vector<Point> ccoords(8), fcoords(4);

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

  Mesh_MSTK mesh("test/hex_3x3x2_ss.exo",MPI_COMM_WORLD,3);


  nv = mesh.num_entities(AmanziMesh::NODE,AmanziMesh::OWNED);
  CHECK_EQUAL(NV,nv);

  for (i = 0; i < nv; i++) {
    Point coords;

    coords.init(mesh.space_dimension()); 

    mesh.node_get_coordinates(i,&coords);
    CHECK_EQUAL(xyz[i][0],coords[0]);
    CHECK_EQUAL(xyz[i][1],coords[1]);
    CHECK_EQUAL(xyz[i][2],coords[2]);
  }


  nf = mesh.num_entities(AmanziMesh::FACE,AmanziMesh::OWNED);
  CHECK_EQUAL(NF,nf);
  
  nc = mesh.num_entities(AmanziMesh::CELL,AmanziMesh::OWNED);
  CHECK_EQUAL(NC,nc);
    
  for (i = 0; i < nc; i++) {
    mesh.cell_get_nodes(i,&cnodes);
    mesh.cell_get_faces(i,&faces);
    mesh.cell_get_face_dirs(i,&facedirs);

    for (j = 0; j < 6; j++) {

      mesh.face_get_nodes(faces[j],&fnodes);
      mesh.face_get_coordinates(faces[j],&fcoords);
      

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
	  CHECK_ARRAY_EQUAL(xyz[expfacenodes[(k0+k)%4]],fcoords[k],3);
	}
      }
      else {
	for (k = 0; k < 4; k++) {
	  CHECK_EQUAL(expfacenodes[(k0+4-k)%4],fnodes[k]);
	  CHECK_ARRAY_EQUAL(xyz[expfacenodes[(k0+4-k)%4]],fcoords[k],3);
	}
      }
      
    }

  
    mesh.cell_get_coordinates(i,&ccoords);
    
    for (j = 0; j < 8; j++)
      CHECK_ARRAY_EQUAL(xyz[cnodes[cnstd[j]]],ccoords[j],3);
  }


  // Verify the sidesets

  int ns;
  ns = mesh.num_sets(AmanziMesh::FACE);
  CHECK_EQUAL(7,ns);

  unsigned int expsetids[7]={1,101,102,103,104,105,106};
  std::vector<Set_ID> setids(7);
  mesh.get_set_ids(AmanziMesh::FACE,&setids);
  
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
    std::vector<Entity_ID> setfaces(16);

    setsize = mesh.get_set_size(setids[i],AmanziMesh::FACE,AmanziMesh::OWNED);
    CHECK_EQUAL(expsetsizes[i],setsize);


    mesh.get_set_entities(setids[i],AmanziMesh::FACE, AmanziMesh::OWNED, &setfaces);
    
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



  std::vector<Entity_ID>  c2f(6);
  Epetra_Map cell_map(mesh.cell_epetra_map(true));
  Epetra_Map face_map(mesh.face_epetra_map(false));

  for (int c=cell_map.MinLID(); c<=cell_map.MaxLID(); c++)
    {
      CHECK_EQUAL(cell_map.GID(c),mesh.GID(c,AmanziMesh::CELL));
      mesh.cell_get_faces(c, &c2f);
      for (int j=0; j<6; j++)
  	{
  	  int f = face_map.LID(mesh.GID(c2f[j],AmanziMesh::FACE));
  	  CHECK( f == c2f[j] );
  	}

    }

}

