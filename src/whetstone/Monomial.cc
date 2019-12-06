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

#include <cstdlib>

#include "Monomial.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Constructor from multi-index
****************************************************************** */
Monomial::Monomial(int d, const int* multi_index, double coef)
{
  d_ = d;
  order_ = 0;
  for (int i = 0; i < d_; ++i) {
    multi_index_[i] = multi_index[i];
    order_ += multi_index_[i];
  }
  coefs_.Reshape(1);
  coefs_(0) = coef;
}


/* ******************************************************************
* Return ordered list of all polynomial coefficients of monomial.
****************************************************************** */
DenseVector Monomial::ExpandCoefficients() const
{
  int size = PolynomialSpaceDimension(d_, order_);
  DenseVector coefs(size);
  coefs.PutScalar(0.0);

  int l = PolynomialPosition(d_, multi_index_);
  coefs(l) = coefs_(0);

  return coefs;
}


/* ******************************************************************
* Calculate monomial value
****************************************************************** */
double Monomial::Value(const AmanziGeometry::Point& xp) const
{
  double tmp = coefs_(0);
  if (tmp != 0.0 && order_ > 0) {
    for (int i = 0; i < d_; ++i) {
      tmp *= std::pow(xp[i] - origin_[i], multi_index_[i]);
    }
  }
  return tmp;
}

}  // namespace WhetStone
}  // namespace Amanzi


