/*
  WhetStone, version 2.1
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Operations with polynomials.
*/

#include "Polynomial.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Constructor
****************************************************************** */
Polynomial::Polynomial(int d, int order) : 
    d_(d), order_(order), it_(d), origin_(d)
{
  size_ = 0;
  coefs_.resize(order_ + 1);
  for (int i = 0; i <= order_; ++i) {
    coefs_[i] = Monomial(d_, i);
    size_ += coefs_[i].size();
  }
}


/* ******************************************************************
* Re-shape polynomial
* NOTE: case d > d_ can be treated more intelligently.
****************************************************************** */
void Polynomial::Reshape(int d, int order, bool reset)
{
  if (d_ != d) {
    d_ = d;
    order_ = order;
    it_.set_dimension(d);
    origin_ = AmanziGeometry::Point(d);

    coefs_.clear();
    coefs_.resize(order_ + 1);
    for (int i = 0; i <= order_; ++i) {
      coefs_[i] = Monomial(d_, i);
    }
  } else if (order_ != order) {
    coefs_.resize(order + 1);
    for (int i = order_ + 1; i <= order; ++i) {
      coefs_[i] = Monomial(d_, i);
    }
    order_ = order;
  }

  size_ = 0;
  for (int i = 0; i <= order_; ++i) {
    size_ += coefs_[i].size();
  }

  if (reset) Reset();
}


/* ******************************************************************
* Implemented ring algebra operations.
* NOTE: implementation is order independent.
****************************************************************** */
Polynomial& Polynomial::operator+=(const Polynomial& poly)
{
  ASSERT(d_ == poly.dimension());  // FIXME
  ASSERT(origin_ == poly.origin());

  int order_max = std::max(order_, poly.order());
  Reshape(d_, order_max);

  for (int i = 0; i <= poly.order(); ++i) {
    auto& mono = coefs_[i].coefs();
    for (int k = 0; k < mono.size(); ++k) {
      mono[k] += poly.monomials(i).coefs()[k];
    }
  }

  return *this;
}


Polynomial& Polynomial::operator-=(const Polynomial& poly)
{
  ASSERT(d_ == poly.dimension());  // FIXME
  ASSERT(origin_ == poly.origin());

  int order_max = std::max(order_, poly.order());
  Reshape(d_, order_max);

  for (int i = 0; i <= poly.order(); ++i) {
    auto& mono = coefs_[i].coefs();
    for (int k = 0; k < mono.size(); ++k) {
      mono[k] -= poly.monomials(i).coefs()[k];
    }
  }

  return *this;
}


Polynomial& Polynomial::operator*=(const Polynomial& poly)
{
  ASSERT(d_ == poly.dimension());  // FIXME
  ASSERT(origin_ == poly.origin());

  Polynomial arg1(*this);
  const Polynomial* arg2 = &poly;
  if (this == arg2) arg2 = &arg1; 

  int order_prod = order_ * poly.order();
  Reshape(d_, order_prod, true);

  int index[3];
  for (auto it1 = arg1.begin(); it1.end() <= arg1.end(); ++it1) {
    const int* idx1 = it1.multi_index();
    int k1 = it1.MonomialOrder();
    int m1 = it1.MonomialPosition();
    double val1 = arg1.monomials(k1).coefs()[m1];

    for (auto it2 = arg2->begin(); it2.end() <= arg2->end(); ++it2) {
      const int* idx2 = it2.multi_index();
      int k2 = it2.MonomialOrder();
      int m2 = it2.MonomialPosition();
      double val2 = arg2->monomials(k2).coefs()[m2];

      int n(0);
      for (int i = 0; i < d_; ++i) {
        index[i] = idx1[i] + idx2[i];
        n += index[i];
      }
      int l = MonomialPosition(index);
      coefs_[n](l) += val1 * val2;
    }
  }

  return *this;
}


/* ******************************************************************
* Compute gradient of a polynomial
***************************************************************** */
void Polynomial::Gradient(std::vector<Polynomial>& grad) const
{
  grad.resize(d_);

  int order = std::max(0, order_ - 1);
  for (int i = 0; i < d_; ++i) {
    grad[i].Reshape(d_, order, true);
  }

  int index[3];
  for (auto it = begin(); it.end() <= end(); ++it) {
    int k = it.MonomialOrder();
    if (k > 0) {
      const int* idx = it.multi_index();
      int m = it.MonomialPosition();
      double val = monomials(k).coefs()[m];

      for (int i = 0; i < d_; ++i) {
        for (int j = 0; j < d_; ++j) index[j] = idx[j];

        if (index[i] > 0) {
          index[i]--;
          m = MonomialPosition(index);
          grad[i].monomials(k - 1)(m) = val * idx[i];
        }
      }
    }
  }
}


/* ******************************************************************
* Rebase polynomial to different origin.
****************************************************************** */
void Polynomial::ChangeOrigin(const AmanziGeometry::Point& origin)
{
  ASSERT(order_ < 2);  // FIXME

  if (order_ == 1) {
    AmanziGeometry::Point shift(origin - origin_);
    for (int i = 0; i < d_; ++i) {
      coefs_[0](0) += coefs_[1](i) * shift[i];
    }
  }
  origin_ = origin;   
}


/* ******************************************************************
* Resets all coefficients to zero
****************************************************************** */
void Polynomial::Reset()
{
  for (int i = 0; i <= order_; ++i) {
    std::vector<double>& tmp = coefs_[i].coefs();
    for (auto it = tmp.begin(); it != tmp.end(); ++it) *it = 0.0;
  }
}


/* ******************************************************************
* Calculate polynomial value
****************************************************************** */
double Polynomial::Value(const AmanziGeometry::Point& xp)
{
  double sum(0.0);

  for (auto it = begin(); it.end() <= end(); ++it) {
    const int* index = it.multi_index();
    int l = MonomialPosition(index);

    int k(0);
    double tmp(1.0);
    for (int i = 0; i < d_; ++i) {
      tmp *= std::pow(xp[i] - origin_[i], index[i]);
      k += index[i];
    }
    sum += tmp * coefs_[k](l);
  }

  return sum;
}



/* ******************************************************************
* Position in monomial defined by multi_index. 2D algorithm
****************************************************************** */
int Polynomial::MonomialPosition(const int* multi_index) const
{
  int m(multi_index[1]);
  if (d_ == 3) {
    int n = multi_index[1] + multi_index[2];
    m = n * (n + 1) / 2 + multi_index[2];
  }
  return m;
}


/* ******************************************************************
* Position in polynomial defined by multi_index. 2D algorithm.
****************************************************************** */
int Polynomial::PolynomialPosition(const int* multi_index) const
{
  // current monomial order
  int k = 0;
  for (int i = 0; i < d_; ++i) k += multi_index[i];

  // sum of previous monomials
  int nk = (d_ == 2) ? k * (k + 1) / 2 : k * (k + 1) * (k + 2) / 6; 

  return nk + MonomialPosition(multi_index);
}


/* ******************************************************************
* Fancy I/O
****************************************************************** */
std::ostream& operator << (std::ostream& os, const Polynomial& p)
{
  for (auto it = p.begin(); it.end() <= p.end(); ++it) {
    int k = it.MonomialOrder();
    int m = it.MonomialPosition();
    double val = p.monomials(k).coefs()[m];

    if (m == 0) os << "k=" << k << "  coefs:";
    if (m > 0 && val >= 0) os << " +";
    os << " " << val << " ";

    const int* index = it.multi_index();
    if (index[0] == 1) os << "x";
    if (index[0] > 1)  os << "x^" << index[0];
    if (index[1] == 1) os << "y";
    if (index[1] > 1)  os << "y^" << index[1];
    if (index[2] == 1) os << "z";
    if (index[2] > 1)  os << "z^" << index[2];

    if (m == p.monomials(k).size() - 1) os << std::endl;
  } 
  os << "origin: " << p.origin() << std::endl;
  return os;
}


}  // namespace WhetStone
}  // namespace Amanzi


