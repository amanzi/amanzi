/*
  The discretization component of Amanzi.
  License: BSD
  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "UnitTest++.h"

#include "Mesh.hh"
#include "MeshFactory.hh"
#include "MeshAudit.hh"
#include "Point.hh"

#include "MFD3D_CrouzeixRaviart.hh"
#include "Tensor.hh"


/* **************************************************************** */
TEST(CROUZEIX_RAVIART) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "Test: High-order Crouzeix Raviart in 2D" << std::endl;
  Epetra_MpiComm *comm = new Epetra_MpiComm(MPI_COMM_WORLD);

  MeshFactory meshfactory(comm);
  meshfactory.preference(FrameworkPreference({MSTK}));
  Teuchos::RCP<Mesh> mesh = meshfactory(0.0, 0.0, 1.0, 1.0, 1, 1); 
 
  MFD3D_CrouzeixRaviart mfd(mesh);

  int nfaces = 4, cell = 0;
  DenseMatrix N, R, Ac, G;

  Tensor T(2, 1);
  T(0, 0) = 1.0;

  mfd.H1consistencyHO(cell, 1, T, N, R, Ac, G);

  printf("Mass matrix for cell %3d\n", cell);
  Ac.PrintMatrix("%8.4f ");

  // verify SPD propery
  for (int i=0; i<nfaces; i++) CHECK(Ac(i, i) > 0.0);

  // verify exact integration property
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;
  mesh->cell_get_faces_and_dirs(cell, &faces, &dirs);
    
  double xi, yi, xj, yj;
  double vxx = 0.0, vxy = 0.0, volume = mesh->cell_volume(cell); 
  for (int i = 0; i < nfaces; i++) {
    int f1 = faces[i];
    for (int j = 0; j < nfaces; j++) {
      int f2 = faces[j];

      xi = mesh->face_centroid(f1)[0];
      yi = mesh->face_centroid(f1)[1];
      xj = mesh->face_centroid(f2)[0];

      vxx += Ac(i, j) * xi * xj;
      vxy += Ac(i, j) * yi * xj;
    }
  }

  // CHECK_CLOSE(T(0,0) * volume, vxx, 1e-10);
  // CHECK_CLOSE(T(0,0) * volume, vxy, 1e-10);

  delete comm;
}


