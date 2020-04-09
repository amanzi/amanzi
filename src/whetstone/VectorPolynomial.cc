/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Operations with vector of polynomials of type p(x - x0) where x0
  could be different for each vector component.
*/

#include <vector>

#include "Point.hh"
#include "MatrixPolynomial.hh"
#include "VectorPolynomial.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
 * Trivial constructor: each component is polynomial=0
 ****************************************************************** */
VectorPolynomial::VectorPolynomial(int d, int size) : d_(d)
{
  polys_.resize(size);
  for (int i = 0; i < size; ++i) polys_[i].Reshape(d_, 0, true);
}


/* ******************************************************************
 * Constructor: each component is zero polynomial of given order
 ****************************************************************** */
VectorPolynomial::VectorPolynomial(int d, int size, int order) : d_(d)
{
  polys_.resize(size);
  for (int i = 0; i < size; ++i) polys_[i].Reshape(d_, order, true);
}


/* ******************************************************************
 * Converter-type constructor
 ****************************************************************** */
VectorPolynomial::VectorPolynomial(const Polynomial& p)
{
  d_ = p.dimension();
  polys_.resize(1);
  polys_[0] = p;
}


/* ******************************************************************
 * Re-shape polynomials
 ****************************************************************** */
void
VectorPolynomial::Reshape(int d, int m, int order, bool reset)
{
  d_ = d;

  polys_.resize(m);
  for (int i = 0; i < m; ++i) { polys_[i].Reshape(d, order, reset); }
}


/* ******************************************************************
 * Reset all coefficients to thesame number
 ****************************************************************** */
void
VectorPolynomial::PutScalar(double val)
{
  for (int i = 0; i < size(); ++i) { polys_[i].PutScalar(val); }
}


/* ******************************************************************
 * Calculate value at a point
 ****************************************************************** */
DenseVector
VectorPolynomial::Value(const AmanziGeometry::Point& xp) const
{
  int n = polys_.size();
  DenseVector val(n);

  for (int i = 0; i < n; ++i) { val(i) = polys_[i].Value(xp); }

  return val;
}


/* ******************************************************************
 * Ring algebra
 ****************************************************************** */
VectorPolynomial&
VectorPolynomial::operator+=(const VectorPolynomial& vp)
{
  for (int i = 0; i < polys_.size(); ++i) { polys_[i] += vp[i]; }
  return *this;
}

VectorPolynomial&
VectorPolynomial::operator-=(const VectorPolynomial& vp)
{
  for (int i = 0; i < polys_.size(); ++i) { polys_[i] -= vp[i]; }
  return *this;
}


/* ******************************************************************
 * Set same origin for all polynomials without modyfying them
 ****************************************************************** */
void
VectorPolynomial::set_origin(const AmanziGeometry::Point& origin)
{
  for (int i = 0; i < size(); ++i) { polys_[i].set_origin(origin); }
}


/* ******************************************************************
 * Change all polynomials to new same origin
 ****************************************************************** */
void
VectorPolynomial::ChangeOrigin(const AmanziGeometry::Point& origin)
{
  for (int i = 0; i < size(); ++i) { polys_[i].ChangeOrigin(origin); }
}


/* ******************************************************************
 * Ring algebra
 ****************************************************************** */
double
VectorPolynomial::NormInf() const
{
  double tmp(0.0);
  for (int i = 0; i < polys_.size(); ++i) {
    tmp = std::max(tmp, polys_[i].NormInf());
  }
  return tmp;
}

} // namespace WhetStone
} // namespace Amanzi
