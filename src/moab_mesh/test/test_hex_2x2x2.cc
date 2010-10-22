#include <UnitTest++.h>

#include <iostream>

#include "../Mesh_maps.hh"
#include "../../mesh_data/Entity_kind.hh"
#include "../Element_category.hh"


#include "Epetra_Map.h"
#include "Epetra_MpiComm.h"

#include "mpi.h"


TEST(MOAB_HEX1)
{

  using namespace std;
  using namespace MOAB_mesh;


  int i, j, k, err, nc, nv;
  int faces[6], nodes[8], facedirs[6];
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
  int cellnodes[8] = {0,1,3,2,4,5,7,6};
  int facenodes[6][4] = {{0,1,5,4},
			 {1,3,7,5},
			 {3,2,6,7},
			 {2,0,4,6},
			 {0,2,3,1},
			 {4,5,7,6}};


  // Load a single hex from the hex1.exo file

  Mesh_maps mesh("hex_2x2x2_ss.exo",MPI_COMM_WORLD);


  nv = mesh.count_entities(Mesh_data::NODE,OWNED);
  CHECK_EQUAL(NV,nv);

  for (i = 0; i < nv; i++) {
    double coords[3];
    
    mesh.node_to_coordinates(i,coords,coords+6);
    CHECK_ARRAY_EQUAL(xyz[i],coords,3);
  }

  
  nc = mesh.count_entities(Mesh_data::CELL,OWNED);
  CHECK_EQUAL(NC,nc);

    
  mesh.cell_to_faces(0,faces,faces+6);
  mesh.cell_to_face_dirs(0,facedirs,facedirs+6);

  for (j = 0; j < 6; j++) {
    mesh.face_to_nodes(faces[j],nodes,nodes+4);
    mesh.face_to_coordinates(faces[j],fcoords,fcoords+12);
      
    for (k = 0; k < 4; k++) {
      CHECK_EQUAL(facenodes[j][k],nodes[k]);
      CHECK_ARRAY_EQUAL(xyz[facenodes[j][k]],&(fcoords[3*k]),3);
    }
  }
      
  mesh.cell_to_nodes(0,nodes,nodes+8);
  mesh.cell_to_coordinates(0,ccoords,ccoords+24);
    
  for (j = 0; j < 8; j++) {
    CHECK_EQUAL(cellnodes[j],nodes[j]);
    CHECK_ARRAY_EQUAL(xyz[cellnodes[j]],&(ccoords[3*j]),3);
  }

}

