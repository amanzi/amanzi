/*
The transport component of the Amanzi code, serial unit tests.
License: BSD
Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>

#include "UnitTest++.h"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

#include "Mesh_simple.hh"
#include "MeshAudit.hh"

#include "mfd3d.hpp"
#include "tensor.hpp"


/* **************************************************************** */
TEST(DARCY_MASS) {
  using namespace Teuchos;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "============ TEST DARCY: MASS MATRIX =====================" << endl;
#ifdef HAVE_MPI
  Epetra_MpiComm *comm = new Epetra_MpiComm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm *comm = new Epetra_SerialComm();
#endif

  int num_components = 3;
  RCP<Mesh> mesh = rcp(new Mesh_simple(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1, 2, 3, comm)); 
 
  MFD3D mfd(mesh);

  int nfaces = 6, cell = 0;
  Tensor T(3, 1);
  T(0, 0) = 1;

  Teuchos::SerialDenseMatrix<int, double> N(nfaces, 3);
  Teuchos::SerialDenseMatrix<int, double> Mc(nfaces, nfaces);
  Teuchos::SerialDenseMatrix<int, double> M(nfaces, nfaces);

  int ok = mfd.L2_consistency(cell, T, N, Mc);
  mfd.stability_scalar(cell, N, Mc, M);

  printf("Mass matrix for cell %3d\n", cell);
  for (int i=0; i<6; i++) {
    for (int j=0; j<6; j++ ) printf("%8.3f ", M(i, j)); 
    printf("\n");
  }

  delete comm;
}


/* **************************************************************** */
TEST(DARCY_INVERSE_MASS) {
  using namespace Teuchos;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "============ TEST DARCY: INVERSE MASS MATRIX ==============" << endl;
#ifdef HAVE_MPI
  Epetra_MpiComm *comm = new Epetra_MpiComm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm *comm = new Epetra_SerialComm();
#endif

  int num_components = 3;
  RCP<Mesh> mesh = rcp(new Mesh_simple(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1, 2, 3, comm)); 
 
  MFD3D mfd(mesh);

  int nfaces = 6, cell = 0;
  Tensor T(3, 1);
  T(0, 0) = 1;

  Teuchos::SerialDenseMatrix<int, double> R(nfaces, 3);
  Teuchos::SerialDenseMatrix<int, double> Wc(nfaces, nfaces);
  Teuchos::SerialDenseMatrix<int, double> W(nfaces, nfaces);

  int ok = mfd.L2_consistency_inverse(cell, T, R, Wc);
  mfd.stability_scalar(cell, R, Wc, W);

  printf("Inverse of mass matrix for cell %3d\n", cell);
  for (int i=0; i<6; i++) {
    for (int j=0; j<6; j++ ) printf("%8.3f ", W(i, j)); 
    printf("\n");
  }

  delete comm;
}



