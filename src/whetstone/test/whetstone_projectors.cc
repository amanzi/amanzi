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
#include "ProjectorH1.hh"


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
  DenseMatrix N, R, G, A1, Ak;

  Tensor T(2, 1);
  T(0, 0) = 1.0;

  // 1st-order scheme
  mfd.StiffnessMatrix(cell, T, A1);

  // 1st-order scheme (new algorithm)
  mfd.StiffnessMatrixHO(cell, 1, T, R, G, Ak);

  printf("Stiffness matrix for order = 1\n");
  Ak.PrintMatrix("%8.4f ");

  A1 -= Ak;
  CHECK(A1.NormInf() <= 1e-10);

  // 2nd-order scheme (new algorithm)
  mfd.H1consistencyHO(cell, 2, T, N, R, G, Ak);
  mfd.StiffnessMatrixHO(cell, 2, T, R, G, Ak);

  printf("Stiffness matrix for order = 2\n");
  Ak.PrintMatrix("%8.4f ");

  // verify SPD propery
  int nrows = Ak.NumRows();
  for (int i = 0; i < nrows; i++) CHECK(Ak(i, i) > 0.0);

  // verify exact integration property
  DenseMatrix G1(G);
  G1.Multiply(N, R, true);
  G1(0, 0) = 1.0;
  G1.Inverse();
  G1(0, 0) = 0.0;
  G1 -= G;
  CHECK(G1.NormInf() <= 1e-10);
    
  delete comm;
}


/* **************************************************************** */
TEST(HARMINIC_PROJECTORS) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "\nTest: High-order projectors in 2D" << std::endl;
  Epetra_MpiComm *comm = new Epetra_MpiComm(MPI_COMM_WORLD);

  MeshFactory meshfactory(comm);
  meshfactory.preference(FrameworkPreference({MSTK}));
  Teuchos::RCP<const Amanzi::AmanziGeometry::GeometricModel> gm;
  Teuchos::RCP<Mesh> mesh = meshfactory(0.0, 0.0, 1.2, 1.1, 1, 1, gm, true, true); 
 
  int cell(0);
  AmanziGeometry::Point p0(0.0, 0.0);
  VectorPolynomial uc;
  std::vector<VectorPolynomial> vf(4);

  // test zero cell deformation
  for (int n = 0; n < 4; ++n) {
    vf[n].resize(2);
    for (int i = 0; i < 2; ++i) {
      vf[n][i].Reshape(2, 1, true);
    }
  }

  ProjectorH1 projector(mesh);
  projector.HarmonicPk_Cell(cell, p0, 2, vf, uc);

  CHECK(uc[0].NormMax() < 1e-12 && uc[1].NormMax() < 1e-12);

  // test re-location of the right-top corner to (2,3)
  p0 = AmanziGeometry::Point(0.2, 0.475);

  vf[1][0](1, 1) = 0.8 / 1.1; 
  vf[1][1](1, 1) = 1.9 / 1.1; 

  vf[2][0](1, 0) = 0.8 / 1.2; 
  vf[2][1](1, 0) = 1.9 / 1.2; 

  projector.HarmonicPk_Cell(cell, p0, 2, vf, uc);
  std::cout << uc[0] << std::endl;
  std::cout << uc[1] << std::endl;

  CHECK(fabs(uc[0](2, 0)) < 1.0e-12 && uc[0](2, 2) < 1.0e-12);
  CHECK(fabs(uc[1](2, 0)) < 1.0e-12 && uc[1](2, 2) < 1.0e-12);

  delete comm;
}
