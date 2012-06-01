
#include "UnitTest++.h"
#include "Teuchos_RCP.hpp"
#include <Epetra_Comm.h>
#include <Epetra_MpiComm.h>
#include "Epetra_SerialComm.h"

#include "../Mesh_simple.hh"

// #include "State.hpp"

TEST(NUMBERING) {

#ifdef HAVE_MPI
  Epetra_MpiComm *comm = new Epetra_MpiComm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm *comm = new Epetra_SerialComm();
#endif

  // Create a single-cell mesh;
  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh(new Amanzi::AmanziMesh::Mesh_simple(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1, 1, 1, comm));
  
  //  State S(1,mesh);
  
  //  std::string gmvfile = "out.gmv";
  //  S.write_gmv(gmvfile);
  
  // Write node coordinates
  std::cout << "Node coordinates..." << std::endl;
  double x[3];
  for (Amanzi::AmanziMesh::Entity_ID j = 0; j < 8; ++j) {
    std::cout << j << ":";
    mesh->node_to_coordinates(j, x, x+3);
    for (int i = 0; i < 3; ++i) std::cout << " " << x[i];
    std::cout << std::endl;
  }
  std::cout << std::endl;
  
  // Write face-node connectivity
  Amanzi::AmanziMesh::Entity_ID fnode[8];
  std::cout << "Face node connectivity..." << std::endl;
  for (Amanzi::AmanziMesh::Entity_ID j = 0; j < 6; ++j) {
    mesh->face_to_nodes(j, fnode, fnode+4);
    std::cout << j << ":";
    for (int i = 0; i < 4; ++i) std::cout << " " << fnode[i];
    std::cout << std::endl;
  }
  std::cout << std::endl;
  
  // Write cell-node connectivity
  Amanzi::AmanziMesh::Entity_ID cnode[8];
  std::cout << "Cell node connectivity..." << std::endl;
  for (Amanzi::AmanziMesh::Entity_ID j = 0; j < 1; ++j) {
    mesh->cell_to_nodes(j, cnode, cnode+8);
    std::cout << j << ":";
    for (int i = 0; i < 8; ++i) std::cout << " " << cnode[i];
    std::cout << std::endl;
  }
  std::cout << std::endl;
  
  // Write cell face-node connectivity
  //  Amanzi::AmanziMesh::Entity_ID cface[6];
  //  int fdir[6];
  Amanzi::AmanziMesh::Entity_ID_List cface;
  std::vector<int> fdir;
  std::cout << "Cell " << 0 << " faces (relative orientation)..." << std::endl;
  mesh->cell_get_faces_and_dirs(0,&cface,&fdir);
  for (int j = 0; j < 6; ++j) std::cout << j << ": " << cface[j] << "(" << fdir[j] << ")" << std::endl;
  
}
