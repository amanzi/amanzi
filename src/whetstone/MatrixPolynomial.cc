/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Operations with matrix of polynomials of type p(x - x0) where x0
  could be different for each matrix entry.
*/

#include <vector>

#include "Point.hh"
#include "MatrixPolynomial.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
 * Trivial constructor: each component is polynomial=0
 ****************************************************************** */
MatrixPolynomial::MatrixPolynomial(int d, int m, int n, int order)
  : d_(d), m_(m), n_(n), order_(order)
{
  polys_.resize(m_);
  for (int i = 0; i < m_; ++i) {
    polys_[i].resize(n_);
    for (int j = 0; j < n_; ++j) { polys_[i][j].Reshape(d_, order, true); }
  }
}


/* ******************************************************************
 * Re-shape polynomials
 ****************************************************************** */
void
MatrixPolynomial::Reshape(int d, int m, int n, int order, bool reset)
{
  d_ = d;
  m_ = m;
  n_ = n;

  polys_.resize(m_);
  for (int i = 0; i < m_; ++i) {
    polys_[i].resize(n_);
    for (int j = 0; j < n_; ++j) { polys_[i][j].Reshape(d, order, reset); }
  }
}


/* ******************************************************************
 * Reset all coefficients to thesame number
 ****************************************************************** */
void
MatrixPolynomial::PutScalar(double val)
{
  for (int i = 0; i < m_; ++i)
    for (int j = 0; j < n_; ++j) polys_[i][j].PutScalar(val);
}


/* ******************************************************************
 * Calculate value at a point
 ****************************************************************** */
DenseMatrix
MatrixPolynomial::Value(const AmanziGeometry::Point& xp) const
{
  DenseMatrix val(m_, n_);

  for (int i = 0; i < m_; ++i)
    for (int j = 0; j < n_; ++j) val(i, j) = polys_[i][j].Value(xp);

  return val;
}


/* ******************************************************************
 * Matrix-vector operations
 ***************************************************************** */
void
MatrixPolynomial::Multiply(const VectorPolynomial& v, VectorPolynomial& av,
                           bool transpose)
{
  if (!transpose) {
    av.resize(m_);

    for (int i = 0; i < m_; ++i) {
      av[i] = polys_[i][0] * v[0];

      for (int k = 1; k < n_; ++k) { av[i] += polys_[i][k] * v[k]; }
    }
  } else {
    av.resize(n_);

    for (int i = 0; i < n_; ++i) {
      av[i] = polys_[0][i] * v[0];

      for (int k = 1; k < m_; ++k) { av[i] += polys_[k][i] * v[k]; }
    }
  }
}

void
MatrixPolynomial::Multiply(const DenseVector& v, VectorPolynomial& av,
                           bool transpose)
{
  if (!transpose) {
    av.resize(m_);
    for (int i = 0; i < m_; ++i) {
      av[i] = polys_[i][0] * v(0);

      for (int k = 1; k < n_; ++k) { av[i] += polys_[i][k] * v(k); }
    }
  } else {
    av.resize(n_);
    for (int i = 0; i < n_; ++i) {
      av[i] = polys_[0][i] * v(0);

      for (int k = 1; k < m_; ++k) { av[i] += polys_[k][i] * v(k); }
    }
  }
}

void
MatrixPolynomial::Multiply(const AmanziGeometry::Point& p, VectorPolynomial& av,
                           bool transpose)
{
  int d(p.dim());
  AMANZI_ASSERT(NumCols() == d);

  if (!transpose) {
    av.resize(m_);
    for (int i = 0; i < m_; ++i) {
      av[i] = polys_[i][0] * p[0];

      for (int k = 1; k < d; ++k) { av[i] += polys_[i][k] * p[k]; }
    }
  } else {
    av.resize(d);
    for (int i = 0; i < d; ++i) {
      av[i] = polys_[0][i] * p[0];

      for (int k = 1; k < d; ++k) { av[i] += polys_[k][i] * p[k]; }
    }
  }
}


/* ******************************************************************
 * Ring algebra
 ****************************************************************** */
MatrixPolynomial&
MatrixPolynomial::operator+=(const MatrixPolynomial& mp)
{
  for (int i = 0; i < m_; ++i)
    for (int j = 0; j < n_; ++j) polys_[i][j] += mp(i, j);

  return *this;
}

MatrixPolynomial&
MatrixPolynomial::operator-=(const MatrixPolynomial& mp)
{
  for (int i = 0; i < m_; ++i)
    for (int j = 0; j < n_; ++j) polys_[i][j] -= mp(i, j);

  return *this;
}


/* ******************************************************************
 * Set same origin for all polynomials without modyfying them
 ****************************************************************** */
void
MatrixPolynomial::set_origin(const AmanziGeometry::Point& origin)
{
  for (int i = 0; i < m_; ++i)
    for (int j = 0; j < n_; ++j) polys_[i][j].set_origin(origin);
}


/* ******************************************************************
 * Change all polynomials to new same origin
 ****************************************************************** */
void
MatrixPolynomial::ChangeOrigin(const AmanziGeometry::Point& origin)
{
  for (int i = 0; i < m_; ++i)
    for (int j = 0; j < n_; ++j) polys_[i][j].ChangeOrigin(origin);
}


/* ******************************************************************
 * Ring algebra
 ****************************************************************** */
double
MatrixPolynomial::NormInf() const
{
  double tmp(0.0);
  for (int i = 0; i < m_; ++i)
    for (int j = 0; j < n_; ++j) tmp = std::max(tmp, polys_[i][j].NormInf());

  return tmp;
}

} // namespace WhetStone
} // namespace Amanzi
