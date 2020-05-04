/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Operations with space-time polynomials of type p0 + t p1 + t^2 p2
  where p0, p1, and p2 are polynomials with possibly different origins.
  This class implements ring algebra for such polynomials. If origins
  are different, operations fails.
*/

#ifndef AMANZI_WHETSTONE_SPACE_TIME_POLYNOMIAL_HH_
#define AMANZI_WHETSTONE_SPACE_TIME_POLYNOMIAL_HH_

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <vector>

#include "Point.hh"
#include "Polynomial.hh"

namespace Amanzi {
namespace WhetStone {

class SpaceTimePolynomial {
 public:
  SpaceTimePolynomial() : d_(0), order_(-1), size_(-1) {};
  SpaceTimePolynomial(int d, int order);
  SpaceTimePolynomial(const SpaceTimePolynomial& poly);
  ~SpaceTimePolynomial() {};

  // reshape polynomial and erase (optionally) memory
  void reshape(int d, int order, bool reset = false);

  // modifiers
  // -- polynomial is recalculated
  void set_origin(const AmanziGeometry::Point& origin) {
    for (int i = 0; i < size_; ++i) coefs_[i].set_origin(origin);
  }
  void ChangeOrigin(const AmanziGeometry::Point& origin) {
    for (int i = 0; i < size_; ++i) coefs_[i].ChangeOrigin(origin);
  }

  // typical operations with polynomials
  // -- polynomial values
  double Value(const AmanziGeometry::Point& xp, double t) const;
  Polynomial<> Value(double t) const;

  // -- operators (extended ring algebra)
  SpaceTimePolynomial& operator+=(const SpaceTimePolynomial& poly);
  SpaceTimePolynomial& operator-=(const SpaceTimePolynomial& poly);
  SpaceTimePolynomial& operator*=(const SpaceTimePolynomial& poly);
  SpaceTimePolynomial& operator*=(double val);
 
  friend SpaceTimePolynomial operator+(const SpaceTimePolynomial& poly1, const SpaceTimePolynomial& poly2) {
    SpaceTimePolynomial tmp(poly1);
    return tmp += poly2;
  }

  friend SpaceTimePolynomial operator-(const SpaceTimePolynomial& poly1, const SpaceTimePolynomial& poly2) {
    SpaceTimePolynomial tmp(poly1);
    return tmp -= poly2;
  }

  friend SpaceTimePolynomial operator*(const SpaceTimePolynomial& poly1, const SpaceTimePolynomial& poly2) {
    SpaceTimePolynomial tmp(poly1);
    return tmp *= poly2;
  }

  friend SpaceTimePolynomial operator*(double val, const SpaceTimePolynomial& poly) {
    SpaceTimePolynomial tmp(poly);
    return tmp *= val;
  }
  friend SpaceTimePolynomial operator*(const SpaceTimePolynomial& poly, double val) {
    SpaceTimePolynomial tmp(poly);
    return tmp *= val;
  }

  // Change of coordinates:
  // --  x = xf + B * s
  void ChangeCoordinates(const AmanziGeometry::Point& xf,
                         const std::vector<AmanziGeometry::Point>& B) {
    for (int i = 1; i < size_; ++i) coefs_[i].ChangeCoordinates(xf, B);
  }
  // --  s = B^+ (x - xf)
  void InverseChangeCoordinates(const AmanziGeometry::Point& xf,
                                const std::vector<AmanziGeometry::Point>& B) {
    for (int i = 1; i < size_; ++i) coefs_[i].InverseChangeCoordinates(xf, B);
  }

  // access
  int dimension() const { return d_; }
  int order() const { return order_; }
  int size() const { return size_; }

  // -- one-index access
  Polynomial<>& operator[](int i) { return coefs_[i]; }
  const Polynomial<>& operator[](int i) const { return coefs_[i]; }

  // norms
  double normInf() const {
    double val(0.0);
    for (int i = 0; i < size_; ++i) std::max(val, coefs_[i].normInf());
    return val;
  }

  // output 
  friend std::ostream& operator << (std::ostream& os, const SpaceTimePolynomial& p);

 protected:
  int d_, order_, size_;
  std::vector<Polynomial<>> coefs_;
};

}  // namespace WhetStone
}  // namespace Amanzi

#endif

