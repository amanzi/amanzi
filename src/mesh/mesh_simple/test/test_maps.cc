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

  Amanzi::AmanziGeometry::Point xc(2.0, 2.0, 2.0);
  Mm.setNodeCoordinate(7,xc);

  std::set<int> expcellnodes{0,1,3,2,4,5,7,6};
  std::vector<Amanzi::AmanziGeometry::Point> expnodecoords;
  expnodecoords.emplace_back(Amanzi::AmanziGeometry::Point{0.0,0.0,0.0});
  expnodecoords.emplace_back(Amanzi::AmanziGeometry::Point{1.0,0.0,0.0});
  expnodecoords.emplace_back(Amanzi::AmanziGeometry::Point{0.0,1.0,0.0});
  expnodecoords.emplace_back(Amanzi::AmanziGeometry::Point{1.0,1.0,0.0});
  expnodecoords.emplace_back(Amanzi::AmanziGeometry::Point{0.0,0.0,1.0});
  expnodecoords.emplace_back(Amanzi::AmanziGeometry::Point{1.0,0.0,1.0});
  expnodecoords.emplace_back(Amanzi::AmanziGeometry::Point{0.0,1.0,1.0});
  expnodecoords.emplace_back(Amanzi::AmanziGeometry::Point{2.0,2.0,2.0});

  int expfacenodes[6][4] = {{0,1,3,2},
                            {4,5,7,6},
                            {0,1,5,4},
                            {2,3,7,6},
                            {0,2,6,4},
                            {1,3,7,5}};


  CHECK_EQUAL(1,Mm.getNumEntities(Amanzi::AmanziMesh::Entity_kind::CELL,Amanzi::AmanziMesh::Parallel_type::OWNED));
  CHECK_EQUAL(6,Mm.getNumEntities(Amanzi::AmanziMesh::Entity_kind::FACE,Amanzi::AmanziMesh::Parallel_type::OWNED));
  CHECK_EQUAL(8,Mm.getNumEntities(Amanzi::AmanziMesh::Entity_kind::NODE,Amanzi::AmanziMesh::Parallel_type::OWNED));

  vector<Amanzi::AmanziGeometry::Point> x(8);
  vector<Amanzi::AmanziMesh::Entity_ID> nodes(8);
  vector<Amanzi::AmanziMesh::Entity_ID> faces(6);

  for (auto i=0; i<Mm.getNumEntities(Amanzi::AmanziMesh::Entity_kind::CELL,Amanzi::AmanziMesh::Parallel_type::OWNED); i++) {
    Mm.getCellNodes(i, nodes);

    CHECK_EQUAL(8,nodes.size());
    std::set<int> nodeset;
    for (const auto& n : nodes) nodeset.insert(n);
    CHECK(expcellnodes == nodeset);

    std::vector<Amanzi::AmanziGeometry::Point> node_coords;
    for (int j=0; j<8; j++) {
      auto nc = Mm.getNodeCoordinate(nodes[j]);
      // check unique
      CHECK(node_coords.end() == std::find(node_coords.begin(), node_coords.end(), nc));
      node_coords.emplace_back(nc);
      // check in the expected list
      CHECK(expnodecoords.end() != std::find(expnodecoords.begin(), expnodecoords.end(), nc));
    }


    Mm.getCellFaces(i, faces);
    for (int j=0; j<6; j++) {
      Amanzi::AmanziMesh::Entity_ID_List fnodes;

      Mm.getFaceNodes(faces[j],fnodes);
      CHECK_ARRAY_EQUAL(expfacenodes[faces[j]],fnodes,4);

      x = Mm.getFaceCoordinates(faces[j]);
	        
      for (int k=0; k<4; k++) {
        CHECK_ARRAY_EQUAL(expnodecoords[expfacenodes[faces[j]][k]],x[k],3);
      }
    }

    x = Mm.getCellCoordinates(i);
    std::vector<Amanzi::AmanziGeometry::Point> cellcoords;
    CHECK_EQUAL(8,x.size());
    for (int k = 0; k < 8; k++) {
      CHECK(expnodecoords.end() != std::find(expnodecoords.begin(), expnodecoords.end(), x[k]));
      CHECK(cellcoords.end() == std::find(cellcoords.begin(), cellcoords.end(), x[k]));
      cellcoords.push_back(x[k]);
    }
  }
}

}
