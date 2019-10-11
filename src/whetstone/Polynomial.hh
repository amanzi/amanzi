/*
  WhetStone, Version 2.2
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
#include "PolynomialBase.hh"
#include "PolynomialIterator.hh"

namespace Amanzi {
namespace WhetStone {

// formard declarations
class Monomial;

int
PolynomialSpaceDimension(int d, int order);

class Polynomial : public PolynomialBase {
 public:
  Polynomial(){};
  Polynomial(int d, int order);
  Polynomial(int d, int order, const DenseVector& coefs);
  Polynomial(int d, const int* multi_index, double factor);
  Polynomial(const Polynomial& poly);
  Polynomial(const Monomial& mono);
  Polynomial(int d, int order, const std::vector<AmanziGeometry::Point>& xyz,
             const DenseVector& values);

  // reshape polynomial with erase (optionally) memory
  void Reshape(int d, int order, bool reset = false);

  // initialization options
  // -- reset all coefficients to a scalar
  void PutScalar(double val) { coefs_.PutScalar(val); }

  // modifiers
  // -- polynomial is recalculated
  void ChangeOrigin(const AmanziGeometry::Point& origin);
  // -- polynomial is created from monomial by changing its origin
  Polynomial
  ChangeOrigin(const Monomial& mono, const AmanziGeometry::Point& origin);

  // typical operations with polynomials
  // -- polynomial values
  virtual double Value(const AmanziGeometry::Point& xp) const override;
  // -- get all polynomial coefficients
  virtual DenseVector ExpandCoefficients() const override { return coefs_; }
  // -- polynomial norms (we use 'inf' instead of 'max' for uniformity)
  double NormInf() const
  {
    double tmp;
    coefs_.NormInf(&tmp);
    return tmp;
  }

  // -- operators (ring algebra)
  Polynomial& operator+=(const Polynomial& poly);
  Polynomial& operator-=(const Polynomial& poly);
  Polynomial& operator*=(const Polynomial& poly);
  Polynomial& operator*=(double val);

  friend Polynomial operator+(const Polynomial& poly1, const Polynomial& poly2)
  {
    Polynomial tmp(poly1);
    return tmp += poly2;
  }

  friend Polynomial operator-(const Polynomial& poly1, const Polynomial& poly2)
  {
    Polynomial tmp(poly1);
    return tmp -= poly2;
  }

  friend Polynomial operator*(const Polynomial& poly1, const Polynomial& poly2)
  {
    Polynomial tmp(poly1);
    return tmp *= poly2;
  }

  friend Polynomial operator*(double val, const Polynomial& poly)
  {
    Polynomial tmp(poly);
    return tmp *= val;
  }
  friend Polynomial operator*(const Polynomial& poly, double val)
  {
    Polynomial tmp(poly);
    return tmp *= val;
  }

  // iterator starts with constant term
  PolynomialIterator begin(int k0 = 0) const
  {
    PolynomialIterator it(d_);
    return it.begin(k0);
  }
  // returns an iterator referring to the past-the-end term in the polynomial
  PolynomialIterator end() const
  {
    PolynomialIterator it(d_);
    return it.begin(order_ + 1);
  }

  // Change of coordinates:
  // --  x = xf + B * s
  void ChangeCoordinates(const AmanziGeometry::Point& xf,
                         const std::vector<AmanziGeometry::Point>& B);
  // --  s = B^+ (x - xf)
  void InverseChangeCoordinates(const AmanziGeometry::Point& xf,
                                const std::vector<AmanziGeometry::Point>& B);

  // -- one-index access
  double& operator()(int i) { return coefs_(i); }
  const double& operator()(int i) const { return coefs_(i); }

  // -- two-index access
  double& operator()(int i, int j)
  {
    return coefs_(PolynomialSpaceDimension(d_, i - 1) + j);
  }
  const double& operator()(int i, int j) const
  {
    return coefs_(PolynomialSpaceDimension(d_, i - 1) + j);
  }

  // output
  friend std::ostream& operator<<(std::ostream& os, const Polynomial& p);

  // specialized member functions
  // -- Laplacian
  Polynomial Laplacian();
};

} // namespace WhetStone
} // namespace Amanzi

#endif
