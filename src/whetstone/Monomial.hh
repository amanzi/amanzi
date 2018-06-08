/*
  WhetStone, version 2.1
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

namespace Amanzi {
namespace WhetStone {

class Monomial {
 public:
  Monomial() : d_(0), order_(-1), coef_(0.0) {};
  Monomial(int d, const int* multi_index, double coef) : d_(d), coef_(coef) {
    order_ = 0;
    for (int i = 0; i < d_; ++i) {
      multi_index_[i] = multi_index[i];
      order_ += multi_index_[i];
    }
  }
  ~Monomial() {};

  // modifiers
  void set_origin(const AmanziGeometry::Point& origin) { origin_ = origin; }

  // typical operations with monomials
  virtual double Value(const AmanziGeometry::Point& xp) const override;

  // -- polynomial norms
  // access
  int dimension() const { return d_; }
  int order() const { return order_; }
  const int* multi_index() const { return multi_index_; }
  const AmanziGeometry::Point& origin() const { return origin_; }

  double& coef() { return coef_; } 
  const double& coef() const { return coef_; } 

 private:
  int d_, order_;
  int multi_index_[3];
  double coef_;
  AmanziGeometry::Point origin_;
};
 

/* ******************************************************************
* Calculate monomial value
****************************************************************** */
double Monomial::Value(const AmanziGeometry::Point& xp) const
{
  double tmp = coef_;
  if (tmp != 0.0) {
    for (int i = 0; i < d_; ++i) {
      tmp *= std::pow(xp[i] - origin_[i], multi_index_[i]);
    }
  }

  return tmp;
}


}  // namespace WhetStone
}  // namespace Amanzi

#endif

