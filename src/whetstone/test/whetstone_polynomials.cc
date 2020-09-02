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
#include "Polynomial.hh"
#include "SpaceTimePolynomial.hh"
#include "SplinePolynomial.hh"
#include "VectorObjects.hh"


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

  // trace of a 2D polynomial
  Polynomial p2d(2, 3);

  for (auto it = p2d.begin(); it < p2d.end(); ++it) {
    const int* index = it.multi_index();
    int pos = PolynomialPosition(2, index);

    int m = it.MonomialSetOrder();
    int k = it.MonomialSetPosition();
    p2d(m, k) = pos;
  }
  AmanziGeometry::Point x0(0.0, 0.0), v1(1.0, 1.0);
  std::vector<AmanziGeometry::Point> tau;
  tau.push_back(v1);

  Polynomial p1d(p2d);
  p1d.ChangeCoordinates(x0, tau);
  std::cout << "tau[0]=" << tau[0] << std::endl;
  std::cout << "Before ChangeCoordinates: " << p2d << std::endl;
  std::cout << "After ChangeCoordinates: " << p1d << std::endl;

  p2d.set_origin(AmanziGeometry::Point(1.0, 1.0));
  p1d = p2d;
  p1d.ChangeCoordinates(x0, tau);
  std::cout << "Before ChangeCoordinates: " << p2d << std::endl;
  std::cout << "After ChangeCoordinates: " << p1d << std::endl;

  // assignement small to large
  q1.Reshape(2, 3, true);
  q2.Reshape(2, 2, true);
  q1 = q2;
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
* Test of splines
**************************************************************** */
TEST(SPLINE_POLYNOMIALS) {
  using namespace Amanzi;
  using namespace Amanzi::WhetStone;

  std::cout << "Splines.." << std::endl; 

  // polynomail x^3
  double xp;
  AmanziGeometry::Point x0(1), x1(1);
  x0[0] = 0.0;
  x1[0] = 1.0;
  xp = 0.2;
  SplineCubic sp0;
  sp0.Setup(x0[0], 0.0, 0.0, x1[0], 1.0, 3.0);
  CHECK_CLOSE(0.008, sp0.Value(xp), 1e-15);

  // general cubic polymonial
  x0[0] = 0.5;
  Polynomial p1(1, 3);
  for (int n = 0; n < 4; ++n) p1(n) = n + 1;
  p1.set_origin(x0);

  auto g1 = Gradient(p1);

  SplineCubic sp1;
  sp1.Setup(x0[0], p1.Value(x0), g1[0].Value(x0),
            x1[0], p1.Value(x1), g1[0].Value(x1));

  p1 -= sp1.poly();
  CHECK_CLOSE(0.0, p1.NormInf(), 1e-15);

  // exterior linear interpolant
  SplineExteriorLinear sp2;
  sp2.Setup(1.0, 1.0, 2.0,  2.0, 3.0, 1.5);
  CHECK_CLOSE(0.0, sp2.Value(0.5), 1e-15);
  CHECK_CLOSE(3.3, sp2.Value(2.2), 1e-15);
}

