/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  WhetStone, Version 2.2
  Release name: naka-to.

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

template <class MEMSPACE = DefaultHostMemorySpace>
class Monomial : public PolynomialBase<MEMSPACE> {
  using PolynomialBase<MEMSPACE>::d_;
  using PolynomialBase<MEMSPACE>::order_;
  using PolynomialBase<MEMSPACE>::size_;
  using PolynomialBase<MEMSPACE>::origin_;
  using PolynomialBase<MEMSPACE>::coefs_;

 public:
  Monomial()
  {
    coefs_.reshape(1);
    coefs_(0) = 0.0;
  }
  Monomial(int d, const int* multi_index, double coef)
  {
    d_ = d;
    order_ = 0;
    for (int i = 0; i < d_; ++i) {
      multi_index_[i] = multi_index[i];
      order_ += multi_index_[i];
    }
    coefs_.reshape(1);
    coefs_(0) = coef;
  }
  ~Monomial(){};

  // typical operations with monomials
  // -- polynomial values
  virtual double Value(const AmanziGeometry::Point& xp) const override;
  // -- get all polynomial coefficients
  virtual DenseVector<MEMSPACE> ExpandCoefficients() const override;

  // access
  const int* multi_index() const { return multi_index_; }

 private:
  int multi_index_[3];
};


/* ******************************************************************
 * Return ordered list of all polynomial coefficients of monomial.
 ****************************************************************** */
template <typename MEMSPACE>
inline DenseVector<MEMSPACE>
Monomial<MEMSPACE>::ExpandCoefficients() const
{
  int size = PolynomialSpaceDimension(d_, order_);
  DenseVector<> coefs(size);
  coefs.putScalar(0.0);

  int l = PolynomialPosition(d_, multi_index_);
  coefs(l) = coefs_(0);

  return coefs;
}


/* ******************************************************************
 * Calculate monomial value
 ****************************************************************** */
template <typename MEMSPACE>
inline double
Monomial<MEMSPACE>::Value(const AmanziGeometry::Point& xp) const
{
  double tmp = coefs_(0);
  if (tmp != 0.0) {
    for (int i = 0; i < d_; ++i) { tmp *= std::pow(xp[i] - origin_[i], multi_index_[i]); }
  }
  return tmp;
}


} // namespace WhetStone
} // namespace Amanzi

#endif
