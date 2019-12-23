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
#include "UnitTest++.h"

// Amanzi::WhetStone
#include "Monomial.hh"
#include "Polynomial.hh"
#include "SpaceTimePolynomial.hh"
#include "VectorObjectsUtils.hh"


/* ****************************************************************
* Test of Taylor polynomials
**************************************************************** */
TEST(DG_TAYLOR_POLYNOMIALS) {
  using namespace Amanzi;
  using namespace Amanzi::WhetStone;

  // polynomials in two dimentions
  Polynomial p(2, 3);

  int i(0);
  for (auto it = p.begin(); it < p.end(); ++it) {
    const int* index = it.multi_index();
    CHECK(index[0] >= 0 && index[1] >= 0);

    int pos = PolynomialPosition(2, index);
    CHECK(pos == i++);

    int m = MonomialSetPosition(2, index);
    p(index[0] + index[1], m) = pos;
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
  for (auto it = q.begin(); it < q.end(); ++it) {
    const int* index = it.multi_index();
    CHECK(index[0] >= 0 && index[1] >= 0 && index[2] >= 0);

    int pos = PolynomialPosition(3, index);
    CHECK(pos == i++);

    int m = MonomialSetPosition(3, index);
    q(index[0] + index[1] + index[2], m) = pos;
  }
  std::cout << "Original polynomial\n" << q << std::endl; 
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
  auto grad = Gradient(q_orig);
  std::cout << "Gradient of a polynomial:\n" << grad << std::endl;
 
  Polynomial lp = q_orig.Laplacian();
  std::cout << "Laplacian of original polynomial:\n" << lp << std::endl;

  q4 = Divergence(grad) - lp; 
  CHECK_CLOSE(0.0, q4.NormInf(), 1e-12);

  // change origin of coordinate system
  AmanziGeometry::Point origin(0.5, 0.3, 0.2);
  q.Reshape(3, 2);
  val = q.Value(xyz);

  q.ChangeOrigin(origin);
  std::cout << "Changed origin of polynomial q\n" << q << std::endl; 
  CHECK_CLOSE(val, q.Value(xyz), 1e-10);

  // assignement small to large
  q1.Reshape(2, 3, true);
  q2.Reshape(2, 2, true);
  q1 = q2;

  // trace of a 2D polynomial on a line (x0, tau1)
  Polynomial p2d(2, 3);

  for (auto it = p2d.begin(); it < p2d.end(); ++it) {
    const int* index = it.multi_index();
    int pos = PolynomialPosition(2, index);

    int m = it.MonomialSetOrder();
    int k = it.MonomialSetPosition();
    p2d(m, k) = pos;
  }
  AmanziGeometry::Point x0(0.0, 0.0), v1(1.0, 1.0);
  std::vector<AmanziGeometry::Point> tau1;
  tau1.push_back(v1);

  Polynomial p1d(p2d);
  p1d.ChangeCoordinates(x0, tau1);
  std::cout << "2D coordinate change: tau[0]=" << tau1[0] << std::endl;
  std::cout << p2d << p1d << std::endl;

  p2d.set_origin(AmanziGeometry::Point(1.0, 1.0));
  p1d = p2d;
  p1d.ChangeCoordinates(x0, tau1);
  std::cout << p2d << p1d << std::endl;

  // trace of a 3D polynomial on a plane (y0, tau2)
  Polynomial p3d(3, 3);

  for (auto it = p3d.begin(); it < p3d.end(); ++it) {
    int pos = it.PolynomialPosition();
    p3d(pos) = pos;
  }
  AmanziGeometry::Point y0(0.0, 0.0, 0.0);
  std::vector<AmanziGeometry::Point> tau2;
  tau2.push_back(AmanziGeometry::Point(1.0, 1.5, 0.0));
  tau2.push_back(AmanziGeometry::Point(2.0,-1.2, 0.0));

  p2d = p3d;
  p2d.ChangeCoordinates(y0, tau2);
  std::cout << "3D coordinate change: tau[0]=" << tau2[0] << "  tau[1]=" << tau2[1] << std::endl;
  std::cout << p3d << p2d << std::endl;
}


/* ****************************************************************
* Test of space-time polynomials
**************************************************************** */
TEST(DG_SPACE_TIME_POLYNOMIALS) {
  using namespace Amanzi;
  using namespace Amanzi::WhetStone;

  std::cout << "Space-time polynomials..." << std::endl; 
  int d(2);
  Polynomial p0(d, 0), p1(d, 1), p2(d, 2);
  p0(0) = 1.0;

  p1(0) = 1.0;
  p1(1) = 2.0;
  p1(2) = 3.0;

  p2(0) = 1.0;
  p2(1) = 2.0;
  p2(5) = 3.0;

  SpaceTimePolynomial stp1(d, 2);
  stp1[0] = p0;
  stp1[1] = p1;
  stp1[2] = p2;

  auto stp2 = stp1;
  auto stp3 = stp1 * stp2;
  auto stp4 = stp2 * stp1;
  stp4 -= stp3;
  CHECK_CLOSE(0.0, stp4.NormInf(), 1e-14);

  stp1 *= 2;
  stp2 = stp1 + stp2;
  // std::cout << stp2 << std::endl;
}


/* ****************************************************************
* Test of vector iterators
**************************************************************** */
TEST(VECTOR_POLYNOMIAL_ITERATOR) {
  using namespace Amanzi;
  using namespace Amanzi::WhetStone;

  // polynomials in two dimentions
  std::cout << "\nVector polynomial iterator..." << std::endl; 
  VectorPolynomial p(2, 2, 3);

  int i(0);
  for (auto it = p.begin(); it < p.end(); ++it) {
    int pos = it.VectorPolynomialPosition();
    CHECK(pos == i++);
    CHECK(it.MonomialSetOrder() <= 3);
    CHECK(it.PolynomialPosition() == pos % 10);
    CHECK(it.VectorComponent() == (pos / 10));
  }
}


/* ****************************************************************
* Test of vector decomposition in 2D
**************************************************************** */
TEST(VECTOR_POLYNOMIAL_DECOMPOSITON_2D) {
  using namespace Amanzi;
  using namespace Amanzi::WhetStone;

  // polynomials in two dimentions
  std::cout << "\nVector polynomial decomposition (2D)..." << std::endl; 

  VectorPolynomial q(2, 2, 2);
  Polynomial p1, p2;

  for (auto it = q.begin(); it < q.end(); ++it) {
    int k = it.VectorComponent();
    int n = it.PolynomialPosition();
    q[k](n) = n + 1;

    VectorDecomposition2DRot(q, p1, p2);

    // verify 
    VectorPolynomial x(2, 2, 1);
    x[0](1) = 1.0;
    x[1](2) = 1.0;

    auto p4 = x;
    for (int i = 0; i < 2; ++i) p4[i] *= p2;

    auto p3 = Rot2D(p1) + p4;
    p3 -= q;
    CHECK_CLOSE(0.0, p3.NormInf(), 1e-12);

    // re-group coefficients
    auto v = ExpandCoefficients(q);
    v.Regroup(6, 10);
  }
}


/* ****************************************************************
* Test of vector decomposition in 3D
**************************************************************** */
TEST(VECTOR_POLYNOMIAL_DECOMPOSITON_3D) {
  using namespace Amanzi;
  using namespace Amanzi::WhetStone;

  // polynomials in two dimentions
  std::cout << "\nVector polynomial decomposition (3D)..." << std::endl; 

  VectorPolynomial q(3, 3, 3), p1;
  Polynomial p2;

  for (auto it = q.begin(); it < q.end(); ++it) {
    int k = it.VectorComponent();
    int n = it.PolynomialPosition();
    Monomial mono(3, it.multi_index(), n + 1);
    mono.set_origin(AmanziGeometry::Point(3));

    VectorDecomposition3DCurl(mono, k, p1, p2);

    // verify vector decomposition
    VectorPolynomial x(3, 3, 1);
    x[0](1) = 1.0;
    x[1](2) = 1.0;
    x[2](3) = 1.0;

    auto p4 = x;
    for (int i = 0; i < 3; ++i) p4[i] *= p2;

    auto p5 = p1 ^ x;
    auto p6 = Curl3D(p5);
    auto p3 = p6 + p4;
    p3[k] -= Polynomial(mono);
    CHECK_CLOSE(0.0, p3.NormInf(), 1e-12);
  }
}


/* ****************************************************************
* Test of matrix form of curl operator
**************************************************************** */
TEST(CURL3D_MATRIX_FORM) {
  using namespace Amanzi;
  using namespace Amanzi::WhetStone;

  // polynomials in two dimentions
  std::cout << "\nMatrix representation of 3D culr operator..." << std::endl; 

  VectorPolynomial q(3, 3, 2);

  for (auto it = q.begin(); it < q.end(); ++it) {
    int k = it.VectorComponent();
    int n = it.PolynomialPosition();
    q[k](n) = k + n;
  }

  auto p = Curl3D(q);

  auto A = Curl3DMatrix(3, 2);
  auto vq = ExpandCoefficients(q);
  auto vp = ExpandCoefficients(p);

  DenseVector tmp(vp);
  A.Multiply(vq, tmp, false);
  vp -= tmp;

  double err, norm;
  vp.NormInf(&err);
  tmp.NormInf(&norm);

  CHECK_CLOSE(0.0, err, 1e-12);
  // PrintMatrix(A);
}

