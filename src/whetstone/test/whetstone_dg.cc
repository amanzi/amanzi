/*
  WhetStone

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <cstdlib>
#include <cmath>
#include <iostream>

// TPLs
#include "Teuchos_RCP.hpp"
#include "UnitTest++.h"

// Amanzi
#include "Mesh.hh"
#include "MeshFactory.hh"

// Amanzi::WhetStone
#include "DenseMatrix.hh"
#include "DG_Modal.hh"
#include "MeshMaps_VEM.hh"
#include "Polynomial.hh"


/* ****************************************************************
* Test of Taylor polynomials
**************************************************************** */
TEST(DG_TAYLOR_POLYNOMIALS) {
  using namespace Amanzi;
  using namespace Amanzi::WhetStone;

  // polynomials in two dimentions
  Polynomial p(2, 3);

  int i(0);
  for (auto it = p.begin(); it.end() <= p.end(); ++it) {
    const int* index = it.multi_index();
    CHECK(index[0] >= 0 && index[1] >= 0);

    int pos = p.PolynomialPosition(index);
    CHECK(pos == i++);

    int m = p.MonomialPosition(index);
    p.monomials(index[0] + index[1]).coefs()[m] = pos;
  }
  std::cout << p << std::endl; 
  CHECK(p.size() == 10);

  // re-define polynomials
  p.Reshape(2, 4);
  std::cout << p << std::endl; 
  CHECK(p.size() == 15);

  Polynomial p_tmp(p);
  p.Reshape(2, 2);
  std::cout << "Reshaping last polynomial\n" << p << std::endl; 
  CHECK(p.size() == 6);

  // operations with polynomials
  AmanziGeometry::Point xy(1.0, 0.0);
  double val = p.Value(xy) + p_tmp.Value(xy);

  p += p_tmp;
  CHECK(p.size() == 15);
  CHECK_CLOSE(p.Value(xy), val, 1e-12);

  // polynomials in 3D
  Polynomial q(3, 3);

  i = 0;
  for (auto it = q.begin(); it.end() <= q.end(); ++it) {
    const int* index = it.multi_index();
    CHECK(index[0] >= 0 && index[1] >= 0 && index[2] >= 0);

    int pos = q.PolynomialPosition(index);
    CHECK(pos == i++);

    int m = q.MonomialPosition(index);
    q.monomials(index[0] + index[1] + index[2]).coefs()[m] = pos;
  }
  std::cout << q << std::endl; 
  CHECK(q.size() == 20);
  Polynomial q_orig(q);

  // reshape polynomials
  q.Reshape(3, 2);
  Polynomial q1(q), q2(q), q3(q);
  std::cout << "Reshaping last 3D polynomial\n" << q << std::endl; 
  CHECK(q.size() == 10);

  q.Reshape(3, 3);
  std::cout << "Reshaping last 3D polynomial, (name q)\n" << q << std::endl; 
  CHECK(q.size() == 20);

  // ring operations with polynomials
  AmanziGeometry::Point xyz(1.0, 2.0, 3.0);
  val = q1.Value(xyz);
  q1 *= q2;
  CHECK_CLOSE(q1.Value(xyz), val * val, 1e-10);

  val = q2.Value(xyz);
  q2 *= q2;
  CHECK_CLOSE(q2.Value(xyz), val * val, 1e-10);

  Polynomial q4 = q2 - q3 * q3;
  CHECK_CLOSE(q4.Value(xyz), 0.0, 1e-10);

  // derivatives
  std::vector<Polynomial> grad;
  q_orig.Gradient(grad);
  std::cout << "Gradient of a polynomial:\n" << grad[0] << grad[1] << grad[2] << std::endl;

  // change coordinate system
  AmanziGeometry::Point origin(0.5, 0.3, 0.2);
  q.Reshape(3, 2);
  val = q.Value(xyz);
  q.ChangeOrigin(origin);
  std::cout << "Changed origin of polynomial q\n" << q << std::endl; 
  CHECK_CLOSE(val, q.Value(xyz), 1e-10);
}


/* ****************************************************************
* Test of DG mass matrices: K is tensor
**************************************************************** */
TEST(DG_MASS_MATRIX) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "Test: DG mass matrices (tensors)" << std::endl;
  Epetra_MpiComm *comm = new Epetra_MpiComm(MPI_COMM_WORLD);

  MeshFactory meshfactory(comm);
  meshfactory.preference(FrameworkPreference({MSTK}));
  Teuchos::RCP<Mesh> mesh = meshfactory(0.0, 0.0, 0.5, 0.5, 1, 1); 

  DenseMatrix M;
  Tensor T(2, 1);
  T(0, 0) = 1.0;

  for (int k = 0; k < 3; k++) {
    DG_Modal dg(k, mesh);

    dg.MassMatrix(0, T, M);
    int nk = M.NumRows();

    printf("Mass matrix for order=%d\n", k);
    for (int i = 0; i < nk; i++) {
      for (int j = 0; j < nk; j++ ) printf("%9.5f ", M(i, j)); 
      printf("\n");
    }

    double area = mesh->cell_volume(0);
    for (int i = 0; i < nk; ++i) {
      CHECK_CLOSE(M(i, i), area, 1e-12);
    }
  }

  delete comm;
}


/* ****************************************************************
* Test of DG mass matrices: K is polynomial
**************************************************************** */
TEST(DG_MASS_MATRIX_POLYNOMIAL) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "\nTest: DG mass matrices (polynomials)" << std::endl;
  Epetra_MpiComm *comm = new Epetra_MpiComm(MPI_COMM_WORLD);

  MeshFactory meshfactory(comm);
  meshfactory.preference(FrameworkPreference({MSTK}));
  // Teuchos::RCP<Mesh> mesh = meshfactory(0.0, 0.0, 1.0, 1.0, 1, 1); 
  Teuchos::RCP<Mesh> mesh = meshfactory("test/one_cell2.exo");
 
  double tmp, integral[3];
  DenseMatrix A;

  for (int k = 0; k < 3; k++) {
    DG_Modal dg(2, mesh);

    Polynomial u(2, k);
    u.monomials(0).coefs()[0] = 1.0;
    u.monomials(k).coefs()[0] = 1.0;

    dg.MassMatrixPoly(0, u, A);
    int nk = A.NumRows();

    printf("Mass matrix for polynomial coefficient: nk=2, pk=%d\n", k);
    for (int i = 0; i < nk; i++) {
      for (int j = 0; j < nk; j++ ) printf("%8.4f ", A(i, j)); 
      printf("\n");
    }

    // TEST1: accuracy (gradient should be rescaled)
    DenseVector v(nk), av(nk);
    const AmanziGeometry::Point& xc = mesh->cell_centroid(0);

    v.PutScalar(0.0);
    v(0) = xc[0] + 2 * xc[1];
    v(1) = 1.0 / 1.6501110800;
    v(2) = 2.0 / 2.6871118178;
    
    A.Multiply(v, av, false);
    v.Dot(av, &tmp);
    integral[k] = tmp;
  }
  CHECK_CLOSE(20.2332916667, integral[0], 1e-10);
  CHECK(integral[0] < integral[1]);

  delete comm;
}


/* ****************************************************************
* Test of DG advection matrices on a face
**************************************************************** */
TEST(DG_ADVECTION_MATRIX_FACE) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "\nTest: DG advection matrices on faces" << std::endl;
  Epetra_MpiComm *comm = new Epetra_MpiComm(MPI_COMM_WORLD);

  MeshFactory meshfactory(comm);
  meshfactory.preference(FrameworkPreference({MSTK}));
  Teuchos::RCP<Mesh> mesh = meshfactory(0.0, 0.0, 1.0, 1.0, 2, 2); 
 
  for (int k = 0; k < 2; k++) {
    DG_Modal dg(k, mesh);
    DenseMatrix A0, A1;

    Polynomial un(2, 0);
    un.monomials(0).coefs()[0] = 1.0;

    // TEST1: constant u
    dg.FluxMatrixPoly(1, un, A0, false);

    printf("Advection matrix (face-based) for order=%d  u.n=1\n", k);
    int nk = A0.NumRows();
    for (int i = 0; i < nk; i++) {
      for (int j = 0; j < nk; j++ ) printf("%8.4f ", A0(i, j)); 
      printf("\n");
    }

    // TEST2: add zero gradient to polynomial un
    un.Reshape(2, 1);
    dg.FluxMatrixPoly(1, un, A1, false);

    A1 -= A0;
    CHECK_CLOSE(0.0, A1.NormInf(), 1e-12);

    // TEST3: nonzero linear component polynomial un
    un.monomials(1).coefs()[0] = 1.0;

    dg.FluxMatrixPoly(1, un, A1, false);

    printf("Advection matrix (face-based) for order=%d u.n=1+x\n", k);
    for (int i = 0; i < nk; i++) {
      for (int j = 0; j < nk; j++ ) printf("%8.4f ", A1(i, j)); 
      printf("\n");
    }
  }

  delete comm;
}


/* ****************************************************************
* Test of DG advection matrices in a cell
**************************************************************** */
TEST(DG_ADVECTION_MATRIX_CELL) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "\nTest: DG advection matrices in cells" << std::endl;
  Epetra_MpiComm *comm = new Epetra_MpiComm(MPI_COMM_WORLD);

  MeshFactory meshfactory(comm);
  meshfactory.preference(FrameworkPreference({MSTK}));
  // Teuchos::RCP<Mesh> mesh = meshfactory(0.0, 0.0, 1.0, 1.0, 2, 2); 
  Teuchos::RCP<Mesh> mesh = meshfactory("test/one_quad.exo"); 

  for (int k = 0; k < 3; k++) {
    DG_Modal dg(k, mesh);
    dg.set_basis(WhetStone::TAYLOR_BASIS_SIMPLE);

    DenseMatrix A0, A1;

    VectorPolynomial u(2);
    for (int i = 0; i < 2; ++i) u[i].Reshape(2, 2);

    // TEST1: constant u
    u[0].monomials(0).coefs()[0] = 1.0;
    u[1].monomials(0).coefs()[0] = 1.0;
    dg.AdvectionMatrixPoly(0, u, A0, false);

    printf("Advection matrix (cell-based) for order=%d u=(1,1)\n", k);
    int nk = A0.NumRows();
    for (int i = 0; i < nk; i++) {
      for (int j = 0; j < nk; j++ ) printf("%10.6f ", A0(i, j)); 
      printf("\n");
    }

    // TEST2: linear u
    u[0].monomials(1).coefs()[0] = 1.0;
    u[0].monomials(1).coefs()[1] = 1.0;
    dg.AdvectionMatrixPoly(0, u, A0, false);

    printf("Advection matrix (cell-based) for order=%d u=(1+x+y,1), f(x,y)=2+x+3y\n", k);
    nk = A0.NumRows();
    for (int i = 0; i < nk; i++) {
      for (int j = 0; j < nk; j++ ) printf("%10.6f ", A0(i, j)); 
      printf("\n");
    }

    // accuracy test for functions 1+x and 1+x
    DenseVector v1(nk), v2(nk), v3(nk);
    if (k > 0) {
      const AmanziGeometry::Point& xc = mesh->cell_centroid(0);
      v1.PutScalar(0.0);
      v1(0) = 2 + xc[0] + 3 * xc[1];
      v1(1) = 1.0;
      v1(2) = 3.0;
      v2 = v1;
 
      A0.Multiply(v1, v3, false);
      double integral(v2 * v3);
      printf("  inner product = %10.6f\n", integral);
      printf("  centroid = %10.6f %10.6f\n", xc[0], xc[1]);

      CHECK_CLOSE(integral, 1891.0 / 48.0, 1e-12);
    }

    // TEST3: quadratic u
    u[1].monomials(2).coefs()[0] = 1.0;
    dg.AdvectionMatrixPoly(0, u, A0, false);

    printf("Advection matrix (cell-based) for order=%d u=(1+x+y,1+x^2)\n", k);
    nk = A0.NumRows();
    for (int i = 0; i < nk; i++) {
      for (int j = 0; j < nk; j++ ) printf("%10.6f ", A0(i, j)); 
      printf("\n");
    }

    if (k > 0) {
      A0.Multiply(v1, v3, false);
      double integral(v2 * v3);
      printf("  inner product = %10.6f\n", integral);
    }
  }

  delete comm;
}

 
/* ****************************************************************
* Test of polynomial least-square approximation
**************************************************************** */
TEST(DG_MAP_APPROXIMATION_CELL) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "\nTest: Polynomial approximation of map in cells." << std::endl;
  Epetra_MpiComm *comm = new Epetra_MpiComm(MPI_COMM_WORLD);

  MeshFactory meshfactory(comm);
  meshfactory.preference(FrameworkPreference({MSTK}));
  Teuchos::RCP<Mesh> mesh = meshfactory("test/one_cell2.exo");

  // extract polygon from the mesh
  Entity_ID_List nodes;
  AmanziGeometry::Point xv;
  std::vector<AmanziGeometry::Point> x1;

  mesh->cell_get_nodes(0, &nodes);

  for (int i = 0; i < nodes.size(); ++i) {
    mesh->node_get_coordinates(nodes[i], &xv);
    x1.push_back(xv);
  }

  // test identity map
  MeshMaps_VEM maps(mesh, mesh);
  VectorPolynomial u;

  maps.LeastSquareFit(1, x1, x1, u);
  for (int i = 0; i < 2; ++i) {
    CHECK_CLOSE(u[i].monomials(0).coefs()[0], 0.0, 1e-12);
    CHECK_CLOSE(u[i].monomials(1).coefs()[i], 1.0, 1e-12);
    CHECK_CLOSE(u[i].monomials(1).coefs()[1 - i], 0.0, 1e-12);
  }

  // test linear map
  std::vector<AmanziGeometry::Point> x2(x1);
  AmanziGeometry::Point shift(0.1, 0.2);
  for (int i = 0; i < nodes.size(); ++i) {
    x2[i] += shift;
  }

  maps.LeastSquareFit(1, x1, x2, u);
  for (int i = 0; i < 2; ++i) {
    CHECK_CLOSE(u[i].monomials(0).coefs()[0], shift[i], 1e-12);
    CHECK_CLOSE(u[i].monomials(1).coefs()[i], 1.0, 1e-12);
    CHECK_CLOSE(u[i].monomials(1).coefs()[1 - i], 0.0, 1e-12);
  }

  // test rotation map
  double s(std::sin(0.3)), c(std::cos(0.3));
  for (int i = 0; i < nodes.size(); ++i) {
    x2[i][0] = c * x1[i][0] - s * x1[i][1];
    x2[i][1] = s * x1[i][0] + c * x1[i][1];
  }

  maps.LeastSquareFit(1, x1, x2, u);
  for (int i = 0; i < 2; ++i) {
    CHECK_CLOSE(u[i].monomials(0).coefs()[0], 0.0, 1e-12);
    CHECK_CLOSE(u[i].monomials(1).coefs()[i], c, 1e-12);
  }
  CHECK_CLOSE(u[0].monomials(1).coefs()[1], -s, 1e-12);
  CHECK_CLOSE(u[1].monomials(1).coefs()[0],  s, 1e-12);

  // test non-linear deformation map
  x1.clear();
  x1.push_back(AmanziGeometry::Point(-0.5, -0.5));
  x1.push_back(AmanziGeometry::Point( 0.5, -0.5));
  x1.push_back(AmanziGeometry::Point(-0.5,  0.5));
  x1.push_back(AmanziGeometry::Point( 0.5,  0.5));

  x2 = x1;
  x2[3] += AmanziGeometry::Point(0.1, 0.1);

  maps.LeastSquareFit(1, x1, x2, u);
  std::cout << u[0] << " " << u[1] << std::endl;

  CHECK_CLOSE(0.025, u[0].monomials(0).coefs()[0], 1e-12);

  delete comm;
}

