/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Monomial is a simple polynomial that is defined completely by
  multi_index, origin, and coefficient. Example of second-order 
  monomial is  m(x,y) = 2 (x-x0) * (y-y0) where (x0, y0) is the
  origin.
*/

#ifndef AMANZI_WHETSTONE_MONOMIAL_HH_
#define AMANZI_WHETSTONE_MONOMIAL_HH_

#include <cstdlib>

#include "Point.hh"
#include "PolynomialBase.hh"

namespace Amanzi {
namespace WhetStone {

class Monomial : public PolynomialBase {
 public:
  Monomial() { coefs_.Reshape(1); coefs_(0) = 0.0; }
  Monomial(int d, const int* multi_index, double coef);
  ~Monomial() {};

  // typical operations with monomials
  // -- polynomial values
  virtual double Value(const AmanziGeometry::Point& xp) const override;
  // -- get all polynomial coefficients
  virtual DenseVector ExpandCoefficients() const override; 

  // access
  const int* multi_index() const { return multi_index_; }

 private:
  int multi_index_[3];
};
 
}  // namespace WhetStone
}  // namespace Amanzi

#endif

