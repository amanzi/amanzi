#include <UnitTest++.h>

#include <iostream>


#include "../Mesh_MOAB.hh"

#include "Epetra_Map.h"
#include "Epetra_MpiComm.h"

#include "mpi.h"


TEST(MOAB_HEX1)
{

  int i, j, k, err, nc, nv;
  unsigned int faces[6], nodes[8];
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
  unsigned int cellnodes[8] = {0,1,3,2,4,5,7,6};
  unsigned int facenodes[6][4] = {{0,1,5,4},
				  {1,3,7,5},
				  {3,2,6,7},
				  {2,0,4,6},
                                  {0,2,3,1},
				  {4,5,7,6}};


  std::auto_ptr<Epetra_MpiComm> comm(new Epetra_MpiComm(MPI_COMM_WORLD));

  // Load a single hex from the hex1.exo file

  Amanzi::AmanziMesh::Mesh_MOAB mesh("test/hex_2x2x2_ss.exo",comm.get());


  // Check number of nodes and their coordinates

  nv = mesh.count_entities(Amanzi::AmanziMesh::NODE, Amanzi::AmanziMesh::OWNED);
  CHECK_EQUAL(NV,nv);

  for (i = 0; i < nv; i++) {
    double coords[3];
    
    mesh.node_to_coordinates(i,coords,coords+3);
    CHECK_ARRAY_EQUAL(xyz[i],coords,3);
  }


  // Check number of cells and their face nodes and their face coordinates
  
  nc = mesh.count_entities(Amanzi::AmanziMesh::CELL,Amanzi::AmanziMesh::OWNED);
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
      
  // Check cell nodes and cell coordinates directly

  mesh.cell_to_nodes(0,nodes,nodes+8);
  mesh.cell_to_coordinates(0,ccoords,ccoords+24);
    
  for (j = 0; j < 8; j++) {
    CHECK_EQUAL(cellnodes[j],nodes[j]);
    CHECK_ARRAY_EQUAL(xyz[cellnodes[j]],&(ccoords[3*j]),3);
  }


  // Verify the sidesets

  int ns;
  ns = mesh.num_sets(Amanzi::AmanziMesh::FACE);
  CHECK_EQUAL(7,ns);

  std::vector<unsigned int> setids(7);
  unsigned int expsetids[7]={1,101,102,103,104,105,106};
  mesh.get_set_ids(Amanzi::AmanziMesh::FACE,&setids);
  
  CHECK_ARRAY_EQUAL(expsetids,setids,7);

  unsigned int setsize, expsetsizes[7] = {6,1,1,1,1,1,1};
  unsigned int expsetfaces[7][6] = {{0,1,2,3,4,5},
				    {0,0,0,0,0,0},
				    {1,0,0,0,0,0},
				    {2,0,0,0,0,0},
				    {3,0,0,0,0,0},
				    {4,0,0,0,0,0},
				    {5,0,0,0,0,0}};


  for (i = 0; i < ns; i++) {
    unsigned int setfaces[6];

    setsize = mesh.get_set_size(setids[i],Amanzi::AmanziMesh::FACE,Amanzi::AmanziMesh::OWNED);
    CHECK_EQUAL(expsetsizes[i],setsize);


    mesh.get_set(setids[i],Amanzi::AmanziMesh::FACE, Amanzi::AmanziMesh::OWNED, setfaces, setfaces+setsize);
    
    CHECK_ARRAY_EQUAL(expsetfaces[i],setfaces,setsize);
  }

  
}

