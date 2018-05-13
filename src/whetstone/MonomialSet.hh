/*
  WhetStone, version 2.1
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Monomial set is the ordered collection of homogeneous polynomials 
  of the same order. In Amanzi, a 'regular' monomial is defined 
  with coefficient one. For example, regular monomials of order two 
  in 2D are x^2, xy and y^2.
*/

#ifndef AMANZI_WHETSTONE_MONOMIAL_SET_HH_
#define AMANZI_WHETSTONE_MONOMIAL_SET_HH_

#include <cstdlib>
#include <vector>

#include "PolynomialIterator.hh"

namespace Amanzi {
namespace WhetStone {

class MonomialSet {
 public:
  MonomialSet() : d_(0), order_(-1) {};
  MonomialSet(int d, int order) : d_(d), order_(order), it_(d) {
    int nk = (order_ == 0) ? 1 : d_;
    for (int i = 1; i < order_; ++i) {
      nk *= d_ + i;
      nk /= i + 1;
    }
    coefs_.resize(nk, 0.0);
  }
  ~MonomialSet() {};

  // iterators
  PolynomialIterator& begin() const { return it_.begin(order_); }
  int end() const { return order_; }

  // reset all coefficients to a scalar
  void PutScalar(double val) {
    for (auto it = coefs_.begin(); it != coefs_.end(); ++it) *it = val;
  }

  // access
  int order() const { return order_; }
  int size() const { return coefs_.size(); }

  std::vector<double>& coefs() { return coefs_; } 
  const std::vector<double>& coefs() const { return coefs_; } 

 private:
  // direct memory access operations
  double& operator()(int i) { return coefs_[i]; }
  const double& operator()(int i) const { return coefs_[i]; }

  friend class Polynomial;

 protected:
  mutable PolynomialIterator it_; 

 private:
  int d_, order_;
  std::vector<double> coefs_;
};
 
}  // namespace WhetStone
}  // namespace Amanzi

#endif

