/*
  WhetStone, version 2.1
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Operations with polynomials of type p(x - x0) where x0 is a
  constant vector. This class implements ring algebra for polynomials 
  and a few useful transformations. See also class VectorPolynomial.
*/

#ifndef AMANZI_WHETSTONE_POLYNOMIAL_HH_
#define AMANZI_WHETSTONE_POLYNOMIAL_HH_

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <vector>

#include "Point.hh"

#include "DenseVector.hh"
#include "MonomialSet.hh"
#include "PolynomialIterator.hh"

namespace Amanzi {
namespace WhetStone {

class Polynomial {
 public:
  Polynomial() : d_(0), order_(-1), size_(0) {};
  Polynomial(int d, int order);
  Polynomial(int d, const int* multi_index, double factor);

  // reshape polynomial with erase (optionally) memory
  void Reshape(int d, int order, bool reset = false);

  // initialization options
  // -- reset all coefficients to a scalar
  void PutScalar(double val);
  // -- set polynomial coefficients from a vector.
  //    The vector size should match that of polynomial.
  void SetPolynomialCoefficients(const DenseVector& coefs);
  // -- copy polynomial coefficients to a vector. 
  //    The vector is resized to accomodate data.
  void GetPolynomialCoefficients(DenseVector& coefs) const;

  // change the coordinate system
  // -- without changing polynomial
  void set_origin(const AmanziGeometry::Point& origin) { origin_ = origin; }
  // -- polynomial is recalculated
  void ChangeOrigin(const AmanziGeometry::Point& origin);

  // typical operations with polynomials
  // -- polynomial values
  double Value(const AmanziGeometry::Point& xp) const;
  // -- polynomial norms
  double NormMax() const;

  // -- operators (ring algebra)
  Polynomial& operator+=(const Polynomial& poly);
  Polynomial& operator-=(const Polynomial& poly);
  Polynomial& operator*=(const Polynomial& poly);
  Polynomial& operator*=(double val);
 
  friend Polynomial operator+(const Polynomial& poly1, const Polynomial& poly2) {
    Polynomial tmp(poly1);
    return tmp += poly2;
  }

  friend Polynomial operator-(const Polynomial& poly1, const Polynomial& poly2) {
    Polynomial tmp(poly1);
    return tmp -= poly2;
  }

  friend Polynomial operator*(const Polynomial& poly1, const Polynomial& poly2) {
    Polynomial tmp(poly1);
    return tmp *= poly2;
  }

  friend Polynomial operator*(double val, const Polynomial& poly) {
    Polynomial tmp(poly);
    return tmp *= val;
  }
  friend Polynomial operator*(const Polynomial& poly, double val) {
    Polynomial tmp(poly);
    return tmp *= val;
  }

  // multi index defines both monomial order and current monomial
  // whenever possible, use faster Iterator's functions 
  int MonomialSetPosition(const int* multi_index) const;
  int PolynomialPosition(const int* multi_index) const;

  // iterator starts with constant term for correct positioning
  PolynomialIterator begin() const { PolynomialIterator it(d_); return it.begin(); }
  int end() const { return order_; }

  // Change of coordinates:
  // --  x = xf + B * s
  void ChangeCoordinates(const AmanziGeometry::Point& xf,
                         const std::vector<AmanziGeometry::Point>& B);
  // --  s = B^+ (x - xf)
  void InverseChangeCoordinates(const AmanziGeometry::Point& xf,
                                const std::vector<AmanziGeometry::Point>& B);

  // access
  int dimension() const { return d_; }
  int order() const { return order_; }
  int size() const { return size_; }
  const AmanziGeometry::Point& origin() const { return origin_; }

  double& operator()(int i, int j) { return coefs_[i](j); }
  const double& operator()(int i, int j) const { return coefs_[i](j); }

  MonomialSet& monomials(int i) { return coefs_[i]; }
  const MonomialSet& monomials(int i) const { return coefs_[i]; }

  // output 
  friend std::ostream& operator << (std::ostream& os, const Polynomial& p);

  // special non-member functions
  // -- Laplacian
  Polynomial Laplacian();

 private:
  int d_, order_, size_;
  AmanziGeometry::Point origin_;
  std::vector<MonomialSet> coefs_;
};

}  // namespace WhetStone
}  // namespace Amanzi

#endif

