#include <iostream>

#include "Epetra_Map.h"
#include "Epetra_MpiComm.h"
#include "UnitTest++.h"

#include "../Mesh_MOAB.hh"

#include "mpi.h"


TEST(MOAB_HEX1)
{
  int i, j, k, err, nc, nv;
  Amanzi::AmanziMesh::Entity_ID_List nodes;
  Amanzi::AmanziMesh::Entity_ID_List faces;
  std::vector<Amanzi::AmanziGeometry::Point> ccoords, fcoords;

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


  std::shared_ptr<Epetra_MpiComm> comm_(new Epetra_MpiComm(MPI_COMM_WORLD));

  // Load a single hex-mesh
  Amanzi::AmanziMesh::Mesh_MOAB mesh("test/hex_1x1x1_ss.exo",comm_.get());

  // Check number of nodes and their coordinates
  nv = mesh.num_entities(Amanzi::AmanziMesh::NODE, Amanzi::AmanziMesh::Parallel_type::OWNED);
  CHECK_EQUAL(NV,nv);

  for (i = 0; i < nv; i++) {
    Amanzi::AmanziGeometry::Point coords;
    
    mesh.node_get_coordinates(i,&coords);
    CHECK_ARRAY_EQUAL(xyz[i],coords,3);
  }


  // Check number of cells and their face nodes and their face coordinates
  nc = mesh.num_entities(Amanzi::AmanziMesh::CELL,Amanzi::AmanziMesh::Parallel_type::OWNED);
  CHECK_EQUAL(NC,nc);

  mesh.cell_get_faces(0,&faces,true);

  for (j = 0; j < 6; j++) {
    mesh.face_get_nodes(faces[j],&nodes);
    mesh.face_get_coordinates(faces[j],&fcoords);

    for (k = 0; k < 4; k++) {
      CHECK_EQUAL(facenodes[j][k],nodes[k]);
      CHECK_ARRAY_EQUAL(xyz[facenodes[j][k]],&(fcoords[k][0]),3);
    }
  }
      
  // Check cell nodes and cell coordinates directly
  mesh.cell_get_nodes(0,&nodes);
  mesh.cell_get_coordinates(0,&ccoords);
    
  for (j = 0; j < 8; j++) {
    CHECK_EQUAL(cellnodes[j],nodes[j]);
    CHECK_ARRAY_EQUAL(xyz[cellnodes[j]],&(ccoords[j][0]),3);
  }

  // Verify the sidesets
  // std::vector<unsigned int> setids(7);
  // unsigned int expsetids[7]={1,101,102,103,104,105,106};
  // mesh.get_set_ids(Amanzi::AmanziMesh::FACE,&setids);
  
  // CHECK_ARRAY_EQUAL(expsetids,setids,7);

  // unsigned int setsize, expsetsizes[7] = {6,1,1,1,1,1,1};
  // unsigned int expsetfaces[7][6] = {{0,1,2,3,4,5},
  //       			    {0,0,0,0,0,0},
  //       			    {1,0,0,0,0,0},
  //       			    {2,0,0,0,0,0},
  //       			    {3,0,0,0,0,0},
  //       			    {4,0,0,0,0,0},
  //       			    {5,0,0,0,0,0}};

  // for (i = 0; i < ns; i++) {
  //   unsigned int setfaces[6];

  //   setsize = mesh.get_set_size(setids[i],Amanzi::AmanziMesh::FACE,Amanzi::AmanziMesh::Parallel_type::OWNED);
  //   CHECK_EQUAL(expsetsizes[i],setsize);

  //   mesh.get_set(setids[i],Amanzi::AmanziMesh::FACE, Amanzi::AmanziMesh::Parallel_type::OWNED, setfaces, setfaces+setsize);
    
  //   CHECK_ARRAY_EQUAL(expsetfaces[i],setfaces,setsize);
  // }
}

