#include <UnitTest++.h>

#include <iostream>

#include "../Mesh_maps.hh"
#include "../../mesh_data/Entity_kind.hh"
#include "../Element_category.hh"


#include "Epetra_Map.h"
#include "Epetra_MpiComm.h"

#include "mpi.h"


TEST(MOAB_HEX_3x3x2)
{

  using namespace std;
  using namespace MOAB_mesh;


  int i, j, k, err, nc, nf, nv;
  int faces[6], nodes[8], facedirs[6];
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
  int cellnodes[4][8] = {{0,1,2,3,4,5,6,7},
			 {3,2,8,9,7,6,10,11},
			 {4,5,6,7,12,13,14,15},
			 {7,6,10,11,15,14,16,17}};
  int facenodes[20][4] = {{0,1,5,4},
			  {1,2,6,5},
			  {2,3,7,6},
			  {3,0,4,7},
			  {0,3,2,1},
			  {4,5,6,7},
			  {2,8,10,6},
			  {8,9,11,10},
			  {9,3,7,11},
			  {3,9,8,2},
			  {7,6,10,11},
			  {4,5,13,12},
			  {5,6,14,13},
			  {6,7,15,14},
			  {7,4,12,15},
			  {12,13,14,15},
			  {6,10,16,14},
			  {10,11,17,16},
			  {11,7,15,17},
			  {15,14,16,17}};


  // Load a single hex from the hex1.exo file

  Mesh_maps mesh("hex_3x3x2_ss.exo",MPI_COMM_WORLD);


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
    mesh.cell_to_faces(i,faces,faces+6);
    mesh.cell_to_face_dirs(i,facedirs,facedirs+6);

    for (j = 0; j < 6; j++) {
      mesh.face_to_nodes(faces[j],nodes,nodes+4);
      mesh.face_to_coordinates(faces[j],fcoords,fcoords+12);
      
      for (k = 0; k < 4; k++) {
	CHECK_EQUAL(facenodes[faces[j]][k],nodes[k]);
	CHECK_ARRAY_EQUAL(xyz[facenodes[faces[j]][k]],&(fcoords[3*k]),3);
      }
    }

  
    mesh.cell_to_nodes(i,nodes,nodes+8);
    mesh.cell_to_coordinates(i,ccoords,ccoords+24);
    
    for (j = 0; j < 8; j++) {
      CHECK_EQUAL(cellnodes[i][j],nodes[j]);
      CHECK_ARRAY_EQUAL(xyz[cellnodes[i][j]],&(ccoords[3*j]),3);
    }
  }

}

