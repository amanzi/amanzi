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
  Teuchos::RCP<const Amanzi::AmanziGeometry::GeometricModel> gm;
  Teuchos::RCP<Mesh> mesh = meshfactory(0.0, 0.0, 1.0, 1.0, 1, 1, gm, true, true); 
 
  MFD3D_CrouzeixRaviart mfd(mesh);

  int cell(0);
  DenseMatrix N, R, Ac, A1, Ak, G;

  Tensor T(2, 1);
  T(0, 0) = 1.0;

  // 1st-order scheme
  mfd.StiffnessMatrix(cell, T, A1);

  // 1st-order scheme (new algorithm)
  mfd.StiffnessMatrixHO(cell, 1, T, Ak);

  printf("Stiffness matrix for order = 1\n");
  Ak.PrintMatrix("%8.4f ");

  A1 -= Ak;
  CHECK(A1.NormInf() <= 1e-10);

  // 2nd-order scheme (new algorithm)
  mfd.H1consistencyHO(cell, 2, T, N, R, Ac, G);
  mfd.StiffnessMatrixHO(cell, 2, T, Ak);

  printf("Stiffness matrix for order = 2\n");
  Ak.PrintMatrix("%8.4f ");

  // verify SPD propery
  int nrows = Ak.NumRows();
  for (int i = 0; i < nrows; i++) CHECK(Ak(i, i) > 0.0);

  // verify exact integration property
  DenseMatrix G1(G);
  G1.Multiply(N, R, true);
  //std::cout << "R="<< R << std::endl;
  //std::cout << "N=" << N << std::endl;
  //std::cout << G1 << std::endl;
    
  delete comm;
}


