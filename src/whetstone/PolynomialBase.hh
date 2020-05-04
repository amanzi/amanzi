/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Virtual helper class for monomials and polynomials. It is used
  mainly in numerical integrators to access special properties of
  polynomials.
*/

#ifndef AMANZI_WHETSTONE_POLYNOMIAL_BASE_HH_
#define AMANZI_WHETSTONE_POLYNOMIAL_BASE_HH_

#include "Point.hh"

#include "DenseVector.hh"
#include "WhetStoneFunction.hh"

namespace Amanzi {
namespace WhetStone {

template<class MEMSPACE = DefaultHostMemorySpace>
class PolynomialBase : public WhetStoneFunction {
 public:
  PolynomialBase() : d_(0), order_(-1), size_(0) {};
  PolynomialBase(int d, int order) : d_(d), order_(order), origin_(d) {};
  virtual ~PolynomialBase() {};

  // convert to regular vector
  virtual DenseVector<MEMSPACE> ExpandCoefficients() const = 0;

  // modifiers
  void set_origin(const AmanziGeometry::Point& origin) { origin_ = origin; }

  // access
  int dimension() const { return d_; }
  int order() const { return order_; }
  int size() const { return size_; }
  const AmanziGeometry::Point& origin() const { return origin_; }
  const DenseVector<MEMSPACE>& coefs() const { return coefs_; }

 protected:
  int d_, order_, size_;
  AmanziGeometry::Point origin_;
  DenseVector<MEMSPACE> coefs_;
};


// calculate dimension of monomial space of given order 
inline
int MonomialSpaceDimension(int d, int order)
{
  int nk = (order == 0) ? 1 : d;
  for (int i = 1; i < order; ++i) {
    nk *= d + i;
    nk /= i + 1;
  }
  return nk;
}

// calculate dimension of polynomial space of given order 
// we assume that space of negative order has dimension zero
inline
int PolynomialSpaceDimension(int d, int order)
{
  if (order < 0) return 0;

  int nk = order + 1;
  for (int i = 1; i < d; ++i) {
    nk *= order + i + 1;
    nk /= i + 1;
  }
  return nk;
}


// Position in set of same-order monomials defined by multi_index.
// Multi index defines both monomial order and monomial position.
inline
int MonomialSetPosition(int d, const int* multi_index)
{
  int m(multi_index[1]);
  if (d == 3) {
    int n = multi_index[1] + multi_index[2];
    m = n * (n + 1) / 2 + multi_index[2];
  }
  return m;
}


// Position in polynomial defined by multi_index.
inline
int PolynomialPosition(int d, const int* multi_index)
{
  // current monomial order
  int k = 0;
  for (int i = 0; i < d; ++i) k += multi_index[i];

  // space of low-order monomials
  int nk = PolynomialSpaceDimension(d, k - 1);

  return nk + MonomialSetPosition(d, multi_index);
}


// calculate order of polynomial space from given dspace dimension 
inline
int PolynomialSpaceOrder(int d, int nk)
{
  int order = -1;
  while (PolynomialSpaceDimension(d, order) != nk) order++;
  return order;
}

} // namespace WhetStone
} // namespace Amanzi

#endif
