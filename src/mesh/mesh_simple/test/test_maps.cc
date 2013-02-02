#include <iostream>
#include "stdlib.h"
#include "math.h"
#include "UnitTest++.h"
#include "../Mesh_simple.hh"
#include <Epetra_Comm.h>
#include <Epetra_MpiComm.h>
#include "Epetra_SerialComm.h"
#include "GenerationSpec.hh"

SUITE (MeshSimple) {
TEST(MAPS) {
  
  using namespace std;

#ifdef HAVE_MPI
  Epetra_MpiComm *comm = new Epetra_MpiComm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm *comm = new Epetra_SerialComm();
#endif


  Amanzi::AmanziMesh::Mesh_simple Mm(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1, 1, 1, comm); 

  double xc[] = { 2.0, 2.0, 2.0 };
  Mm.node_set_coordinates(7,xc);

  int expcellnodes[8] = {0,1,3,2,4,5,7,6};
  int expnodecoords[8][3] = {{0.0,0.0,0.0},
                             {1.0,0.0,0.0},
                             {0.0,1.0,0.0},
                             {1.0,1.0,0.0},
                             {0.0,0.0,1.0},
                             {1.0,0.0,1.0},
                             {0.0,1.0,1.0},
                             {2.0,2.0,2.0}};
  int expfacenodes[6][4] = {{0,1,3,2},
                            {4,5,7,6},
                            {0,1,5,4},
                            {2,3,7,6},
                            {0,2,6,4},
                            {1,3,7,5}};

  int expfacedirs[6] = {-1,1,1,-1,-1,1};
                              

  CHECK_EQUAL(1,Mm.num_entities(Amanzi::AmanziMesh::CELL,Amanzi::AmanziMesh::OWNED));
  CHECK_EQUAL(6,Mm.num_entities(Amanzi::AmanziMesh::FACE,Amanzi::AmanziMesh::OWNED));
  CHECK_EQUAL(8,Mm.num_entities(Amanzi::AmanziMesh::NODE,Amanzi::AmanziMesh::OWNED));

  vector<Amanzi::AmanziGeometry::Point> x(8);
  vector<Amanzi::AmanziMesh::Entity_ID> nodes(8);
  vector<Amanzi::AmanziMesh::Entity_ID> faces(6);
  vector<int> face_dirs(6);
  
  for (Amanzi::AmanziMesh::Entity_ID i=0; i<Mm.num_entities(Amanzi::AmanziMesh::CELL,Amanzi::AmanziMesh::OWNED); i++)
    {

      Mm.cell_get_nodes(i, &nodes);

      CHECK_EQUAL(8,nodes.size());
      CHECK_ARRAY_EQUAL(expcellnodes,nodes,8);
      
      for (int j=0; j<8; j++) {
	Mm.node_get_coordinates(nodes[j],&(x[j]));
        CHECK_ARRAY_EQUAL(expnodecoords[expcellnodes[j]],x[j],3);
      }

      Mm.cell_get_faces_and_dirs(i, &faces, &face_dirs, true);
      double xx[4][3];
      for (int j=0; j<6; j++) {
        Amanzi::AmanziMesh::Entity_ID_List fnodes;

        Mm.face_get_nodes(faces[j],&fnodes);
        CHECK_ARRAY_EQUAL(expfacenodes[faces[j]],fnodes,4);

	Mm.face_get_coordinates(faces[j],&x);
	        
	for (int k=0; k<4; k++) {
          CHECK_ARRAY_EQUAL(expnodecoords[expfacenodes[faces[j]][k]],x[k],3);
        }

      }


      Mm.cell_get_coordinates(i, &x);
      CHECK_EQUAL(8,x.size());
      for (int k = 0; k < 8; k++)
        CHECK_ARRAY_EQUAL(expnodecoords[expcellnodes[k]],x[k],3);

    }
}

}
