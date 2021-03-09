
#include "UnitTest++.h"
#include "Teuchos_RCP.hpp"

#include <AmanziComm.hh>
#include "Epetra_SerialComm.h"

#include "../Mesh_simple.hh"

TEST(NUMBERING) {
  auto comm = Amanzi::getDefaultComm();

  double expnodecoords[8][3] = {{0.0,0.0,0.0},
                                {1.0,0.0,0.0},
                                {0.0,1.0,0.0},
                                {1.0,1.0,0.0},
                                {0.0,0.0,1.0},
                                {1.0,0.0,1.0},
                                {0.0,1.0,1.0},
                                {1.0,1.0,1.0}};
  int expcellnodes[8] = {0,1,3,2,4,5,7,6};
  int expfacenodes[6][4] = {{0,1,3,2},
                            {4,5,7,6},
                            {0,1,5,4},
                            {2,3,7,6},
                            {0,2,6,4},
                            {1,3,7,5}};

  int expcellfaces[6] = {2,5,3,4,0,1};
  int expfacedirs[6] = {1,1,-1,-1,-1,1};
                              
  int expcelledges[12] = {0,1,3,2,4,5,7,6,8,9,11,10};

  // Create a single-cell mesh;
  auto plist = Teuchos::rcp(new Teuchos::ParameterList("my simple mesh with edges"));
  plist->set("request edges", true);
  auto mesh = Teuchos::rcp(new Amanzi::AmanziMesh::Mesh_simple(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1, 1, 1, comm, Teuchos::null, plist));
  
  // Write node coordinates
  Amanzi::AmanziGeometry::Point x;
  for (Amanzi::AmanziMesh::Entity_ID j = 0; j < 8; ++j) {
    x = mesh->getNodeCoordinate(j);
    CHECK_ARRAY_EQUAL(expnodecoords[j],x,3);
  }
  
  // Write face-node connectivity
  Amanzi::AmanziMesh::Entity_ID_List fnode;
  for (Amanzi::AmanziMesh::Entity_ID j = 0; j < 6; ++j) {
    mesh->getFaceNodes(j, fnode);
    CHECK_EQUAL(4,fnode.size());
    CHECK_ARRAY_EQUAL(expfacenodes[j],fnode,fnode.size());
  }
  
  // Write cell-node connectivity
  Amanzi::AmanziMesh::Entity_ID_List cnode;
  mesh->getCellNodes(0, cnode);
  CHECK_EQUAL(8,cnode.size());
  std::set<int> expcellnode_set, cnode_set;
  for (const auto& n : expcellnodes) expcellnode_set.insert(n);
  for (const auto& n : cnode) cnode_set.insert(n);
  CHECK(expcellnode_set == cnode_set);
  
  // Write cell face-node connectivity
  // Amanzi::AmanziMesh::Entity_ID cface[6];
  // int fdir[6];
  Amanzi::AmanziMesh::Entity_ID_List cface;
  std::vector<int> fdir;
  mesh->getCellFacesAndDirs(0,cface,&fdir);
  CHECK_ARRAY_EQUAL(expcellfaces,cface,6);
  CHECK_ARRAY_EQUAL(expfacedirs,fdir,6);
  
  // Write cell-edge connectivity
  Amanzi::AmanziMesh::Entity_ID_List cedge;
  mesh->getCellEdges(0, cedge);
  CHECK_EQUAL(12,cedge.size());

  std::set<int> expcelledge_set, cedge_set;
  for (const auto& n : expcelledges) expcelledge_set.insert(n);
  for (const auto& n : cedge) cedge_set.insert(n);
  CHECK(expcelledge_set == cedge_set);
}
