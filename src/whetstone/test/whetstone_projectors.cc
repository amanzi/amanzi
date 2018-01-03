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
  // Teuchos::RCP<Mesh> mesh = meshfactory(0.0, 0.0, 1.2, 1.1, 1, 1, gm, true, true); 
  Teuchos::RCP<Mesh> mesh = meshfactory("test/one_pentagon.exo", gm, true, true); 
 
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
  A1.PrintMatrix("%8.4f ");

  A1 -= Ak;
  CHECK(A1.NormInf() <= 1e-10);

  // high-order scheme (new algorithm)
  for (int k = 2; k < 4; ++k) {
    mfd.H1consistencyHO(cell, k, T, N, R, G, Ak);
    mfd.StiffnessMatrixHO(cell, k, T, R, G, Ak);

    printf("Stiffness matrix for order = %d\n", k);
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
  }
 
  delete comm;
}


/* **************************************************************** */
TEST(HARMONIC_PROJECTORS_SQUARE) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "\nTest: High-order projectors for square (linear deformation)" << std::endl;
  Epetra_MpiComm *comm = new Epetra_MpiComm(MPI_COMM_WORLD);

  MeshFactory meshfactory(comm);
  meshfactory.preference(FrameworkPreference({MSTK}));
  Teuchos::RCP<const Amanzi::AmanziGeometry::GeometricModel> gm;
  Teuchos::RCP<Mesh> mesh = meshfactory(0.0, 0.0, 1.2, 1.1, 1, 1, gm, true, true); 
 
  int cell(0);
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
  projector.HarmonicPk_Cell(cell, 2, vf, uc);

  CHECK(uc[0].NormMax() < 1e-12 && uc[1].NormMax() < 1e-12);

  // test linear deformation
  for (int n = 0; n < 4; ++n) {
    for (int i = 0; i < 2; ++i) {
      vf[n][i](0, 0) = 1.0;
      vf[n][i](1, 0) = 2.0;
      vf[n][i](1, 1) = 3.0;
    }
  }
  
  projector.HarmonicP1_Cell(cell, vf, uc);
  std::cout << uc[0] << std::endl;

  uc[0] -= vf[0][0];
  uc[1] -= vf[0][1];
  CHECK(uc[0].NormMax() < 1e-12 && uc[1].NormMax() < 1e-12);

  for (int k = 2; k < 4; ++k) {
    projector.HarmonicPk_Cell(cell, k, vf, uc);
    std::cout << uc[1] << std::endl;

    uc[0] -= vf[0][0];
    uc[1] -= vf[0][1];
    CHECK(uc[0].NormMax() < 1e-12 && uc[1].NormMax() < 1e-12);
  }

  // test re-location of the right-top corner to (2,3)
  std::cout << "Test: High-order projectors for square (bilinear deformation)" << std::endl;
  for (int n = 0; n < 4; ++n) vf[n].PutScalar(0.0);
  vf[1][0](1, 1) = 0.8 / 1.1; 
  vf[1][1](1, 1) = 1.9 / 1.1; 

  vf[2][0](1, 0) = 0.8 / 1.2; 
  vf[2][1](1, 0) = 1.9 / 1.2; 

  projector.HarmonicPk_Cell(cell, 2, vf, uc);

  std::cout << uc[0] << std::endl;
  std::cout << uc[1] << std::endl;
  auto p = AmanziGeometry::Point(1.2, 1.1);
  CHECK(fabs(uc[0].Value(p) - 0.8) < 1e-12 &&
        fabs(uc[1].Value(p) - 1.9) < 1e-12);

  // add additive constant deformation 
  std::cout << "Test: High-order projectors for square (2nd bilinear deformation)" << std::endl;
  for (int n = 0; n < 4; ++n) vf[n][0](0, 0) = 1.0;
  projector.HarmonicPk_Cell(cell, 2, vf, uc);
  std::cout << uc[0] << std::endl;

  projector.HarmonicPk_Cell(cell, 3, vf, uc);
  std::cout << uc[0] << std::endl;

  delete comm;
}


/* **************************************************************** */
TEST(HARMONIC_PROJECTORS_POLYGON) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "\nTest: High-order projectors for pentagon (linear deformation)" << std::endl;
  Epetra_MpiComm *comm = new Epetra_MpiComm(MPI_COMM_WORLD);

  MeshFactory meshfactory(comm);
  meshfactory.preference(FrameworkPreference({MSTK}));
  Teuchos::RCP<const Amanzi::AmanziGeometry::GeometricModel> gm;
  Teuchos::RCP<Mesh> mesh = meshfactory("test/one_pentagon.exo", gm, true, true);
  // Teuchos::RCP<Mesh> mesh = meshfactory("test/one_quad.exo", gm, true, true);
 
  int cell(0), nfaces(5);
  VectorPolynomial uc;
  std::vector<VectorPolynomial> vf(nfaces);

  // test linear deformation
  for (int n = 0; n < nfaces; ++n) {
    vf[n].resize(2);
    for (int i = 0; i < 2; ++i) {
      vf[n][i].Reshape(2, 1, true);
      vf[n][i](0, 0) = 1.0;
      vf[n][i](1, 0) = 2.0;
      vf[n][i](1, 1) = 3.0;
    }
  }
  
  // -- old scheme
  ProjectorH1 projector(mesh);
  projector.HarmonicP1_Cell(cell, vf, uc);
  std::cout << uc[0] << std::endl;

  uc[0] -= vf[0][0];
  uc[1] -= vf[0][1];
  CHECK(uc[0].NormMax() < 1e-12 && uc[1].NormMax() < 1e-12);

  // -- new scheme (k=1)
  for (int k = 1; k < 4; ++k) {
    projector.HarmonicPk_Cell(cell, k, vf, uc);
    std::cout << uc[0] << std::endl;

    uc[0] -= vf[0][0];
    uc[1] -= vf[0][1];
    CHECK(uc[0].NormMax() < 1e-12 && uc[1].NormMax() < 1e-12);
  }

  // test quadratic deformation
  std::cout << "\nTest: High-order projectors for pentagon (quadratic deformation)" << std::endl;
  for (int n = 0; n < nfaces; ++n) {
    for (int i = 0; i < 2; ++i) {
      vf[n][i].Reshape(2, 2, false);
      vf[n][i](2, 0) = 4.0;
      vf[n][i](2, 1) = 5.0;
      vf[n][i](2, 2) = -4.0;
    }
  }

  for (int k = 2; k < 4; ++k) {
    projector.HarmonicPk_Cell(cell, k, vf, uc);
    std::cout << uc[0] << std::endl;
    uc[0] -= vf[0][0];
    uc[1] -= vf[0][1];
    CHECK(uc[0].NormMax() < 1e-10 && uc[1].NormMax() < 1e-10);
  }

  // test cubic deformation
  std::cout << "\nTest: High-order projectors for pentagon (cubic deformation)" << std::endl;
  for (int n = 0; n < nfaces; ++n) {
    for (int i = 0; i < 2; ++i) {
      vf[n][i].Reshape(2, 3, false);
      vf[n][i](3, 0) = 2.0;
      vf[n][i](3, 1) = -6.0;
      vf[n][i](3, 2) = -6.0;
      vf[n][i](3, 3) = 2.0;
    }
  }

  for (int k = 3; k < 4; ++k) {
    projector.HarmonicPk_Cell(cell, k, vf, uc);
    std::cout << vf[0][0] << std::endl;
    std::cout << uc[0] << std::endl;
    uc[0] -= vf[0][0];
    uc[1] -= vf[0][1];
    CHECK(uc[0].NormMax() < 1e-12 && uc[1].NormMax() < 1e-12);
  }

  delete comm;
}
