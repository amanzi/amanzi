/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Operations with space-time polynomials of special structure.
*/

#include <cmath>

#include "Polynomial.hh"
#include "SpaceTimePolynomial.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Constructor of a valid zero polynomial.
****************************************************************** */
SpaceTimePolynomial::SpaceTimePolynomial(int d, int order)
{
  d_ = d;
  order_ = order;
  size_ = order + 1;

  Polynomial<> tmp(d, 0);
  coefs_.resize(size_, tmp);
}


/* ******************************************************************
* Copy constructor.
****************************************************************** */
SpaceTimePolynomial::SpaceTimePolynomial(const SpaceTimePolynomial& poly)
{
  d_ = poly.dimension();
  order_ = poly.order();
  size_ = poly.size();

  coefs_.resize(size_);
  for (int i = 0; i < size_; ++i) coefs_[i] = poly[i];
}


/* ******************************************************************
* Re-shape polynomial
* NOTE: case d > d_ can be treated more intelligently.
****************************************************************** */
void SpaceTimePolynomial::reshape(int d, int order, bool reset)
{
  if (d_ != d) {
    d_ = d;
    order_ = order;
    size_ = order + 1;

    Polynomial<> tmp(d, 0);
    coefs_.resize(size_, tmp);
  } else if (order_ != order) {
    int size = size_;
    order_ = order;
    size_ = order_ + 1;

    coefs_.resize(size_);

    if (reset) { 
      Polynomial<> tmp(d, 0);
      for (int i = 0; i < size_; ++i) coefs_[i] = tmp;
    } else {
      Polynomial<> tmp(d, 0);
      for (int i = size; i < size_; ++i) coefs_[i] = tmp;
    }
  } else if (reset) {
    Polynomial<> tmp(d, 0);
    for (int i = 0; i < size_; ++i) coefs_[i] = tmp;
  }
}


/* ******************************************************************
* Implemented ring algebra operations.
****************************************************************** */
SpaceTimePolynomial& SpaceTimePolynomial::operator+=(const SpaceTimePolynomial& poly)
{
  AMANZI_ASSERT(d_ == poly.dimension());

  int order = poly.order();
  if (order_ < order) reshape(d_, order);
  for (int i = 0; i < poly.size(); ++i) coefs_[i] += poly[i];

  return *this;
}


SpaceTimePolynomial& SpaceTimePolynomial::operator-=(const SpaceTimePolynomial& poly)
{
  AMANZI_ASSERT(d_ == poly.dimension());

  int order = poly.order();
  if (order_ < order) reshape(d_, order);
  for (int i = 0; i < poly.size(); ++i) coefs_[i] -= poly[i];

  return *this;
}


SpaceTimePolynomial& SpaceTimePolynomial::operator*=(const SpaceTimePolynomial& poly)
{
  AMANZI_ASSERT(d_ == poly.dimension());

  int order = poly.order();
  int order_prod = order_ + order;
  SpaceTimePolynomial product(d_, order_prod);
  product.set_origin(poly[0].origin());

  for (int i = 0; i < size_; ++i) {
    for (int j = 0; j < poly.size(); ++j) {
      product[i + j] += coefs_[i] * poly[j];
    }
  }

  *this = product;
  return *this;
}


SpaceTimePolynomial& SpaceTimePolynomial::operator*=(double val) {
  for (int i = 0; i < size_; ++i) coefs_[i] *= val;
  return *this;
}


/* ******************************************************************
* Calculate polynomial value at a given point. 
****************************************************************** */
double SpaceTimePolynomial::Value(const AmanziGeometry::Point& xp, double t) const
{
  double sum(coefs_[0](0)), tmp(t);
  for (int i = 1; i < size_; ++i) {
    sum += coefs_[i].Value(xp) * tmp;
    tmp *= t;
  }
  
  return sum;
}


/* ******************************************************************
* Calculate polynomial value at a given time point.
****************************************************************** */
Polynomial<> SpaceTimePolynomial::Value(double t) const
{
  double tmp(t);
  auto poly = coefs_[0];

  for (int i = 1; i < size_; ++i) {
    poly += coefs_[i] * tmp;
    tmp *= t;
 }

  return poly;
}


/* ******************************************************************
* Fancy I/O
****************************************************************** */
std::ostream& operator << (std::ostream& os, const SpaceTimePolynomial& p)
{
  int d = p.dimension();
  os << "space-time polynomial: order=" << p.order() << " d=" << d
     << " time_terms=" << p.size() << std::endl;
  for (int i = 0; i < p.size(); ++i) {
    os << "t^" << i << " * " << p[i];
  } 
  return os;
}

}  // namespace WhetStone
}  // namespace Amanzi


