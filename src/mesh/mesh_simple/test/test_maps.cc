#include <iostream>
#include "stdlib.h"
#include "math.h"
#include "UnitTest++.h"
#include "../Mesh_simple.hh"

#include <AmanziComm.hh>
#include "Epetra_SerialComm.h"
#include "GenerationSpec.hh"

SUITE (MeshSimple) {
TEST(MAPS) {
  
  using namespace std;

  auto comm = Amanzi::getDefaultComm();

  Amanzi::AmanziMesh::Mesh_simple Mm(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1, 1, 1, comm); 

  double xc[] = { 2.0, 2.0, 2.0 };
  Mm.node_set_coordinates(7,xc);

  int expcellnodes[8] = {0,1,3,2,4,5,7,6};
  double expnodecoords[8][3] = {{0.0,0.0,0.0},
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


  CHECK_EQUAL(1,Mm.num_entities(Amanzi::AmanziMesh::CELL,Amanzi::AmanziMesh::Parallel_type::OWNED));
  CHECK_EQUAL(6,Mm.num_entities(Amanzi::AmanziMesh::FACE,Amanzi::AmanziMesh::Parallel_type::OWNED));
  CHECK_EQUAL(8,Mm.num_entities(Amanzi::AmanziMesh::NODE,Amanzi::AmanziMesh::Parallel_type::OWNED));

  vector<Amanzi::AmanziGeometry::Point> x(8);
  vector<Amanzi::AmanziMesh::Entity_ID> nodes(8);
  vector<Amanzi::AmanziMesh::Entity_ID> faces(6);
  
  for (auto i=0; i<Mm.num_entities(Amanzi::AmanziMesh::CELL,Amanzi::AmanziMesh::Parallel_type::OWNED); i++) {
    Mm.cell_get_nodes(i, &nodes);

    CHECK_EQUAL(8,nodes.size());
    CHECK_ARRAY_EQUAL(expcellnodes,nodes,8);
      
    for (int j=0; j<8; j++) {
      Mm.node_get_coordinates(nodes[j],&(x[j]));
      CHECK_ARRAY_EQUAL(expnodecoords[expcellnodes[j]],x[j],3);
    }

    Mm.cell_get_faces(i, &faces, true);
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
