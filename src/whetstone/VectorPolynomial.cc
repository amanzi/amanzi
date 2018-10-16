/*
  WhetStone, version 2.1
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
* Reset all coefficients to thesame number
****************************************************************** */
void VectorPolynomial::PutScalar(double val)
{
  for (int i = 0; i < size(); ++i) {
    polys_[i].PutScalar(val);
  }
}


/* ******************************************************************
* Calculate value at a point 
****************************************************************** */
DenseVector VectorPolynomial::Value(const AmanziGeometry::Point& xp) const
{
  int n = polys_.size();
  DenseVector val(n);

  for (int i = 0; i < n; ++i) {
    val(i) = polys_[i].Value(xp);
  }

  return val;
}


/* ******************************************************************
* Matrix-vector operations
***************************************************************** */
void VectorPolynomial::Multiply(const std::vector<std::vector<Polynomial> >& A, 
                                const VectorPolynomial& v, bool transpose)
{
  int nrows(A.size());
  int ncols(v.size());

  if (!transpose) {
    resize(nrows);

    for (int i = 0; i < nrows; ++i) {
      polys_[i] = A[i][0] * v[0];

      for (int k = 1; k < ncols; ++k) {
        polys_[i] += A[i][k] * v[k];
      }
    }
  } else {
    resize(ncols);

    for (int i = 0; i < ncols; ++i) {
      polys_[i] = A[0][i] * v[0];

      for (int k = 1; k < nrows; ++k) {
        polys_[i] += A[k][i] * v[k];
      }
    }
  }
}


void VectorPolynomial::Multiply(const std::vector<std::vector<Polynomial> >& A, 
                                const AmanziGeometry::Point& p, bool transpose)
{
  int d(p.dim());
  AMANZI_ASSERT(A.size() == d);

  resize(d);
  if (!transpose) {
    for (int i = 0; i < d; ++i) {
      polys_[i] = A[i][0] * p[0];

      for (int k = 1; k < d; ++k) {
        polys_[i] += A[i][k] * p[k];
      }
    }
  } else {
    for (int i = 0; i < d; ++i) {
      polys_[i] = A[0][i] * p[0];

      for (int k = 1; k < d; ++k) {
        polys_[i] += A[k][i] * p[k];
      }
    }
  }
}


/* ******************************************************************
* Ring algebra
****************************************************************** */
VectorPolynomial& VectorPolynomial::operator+=(const VectorPolynomial& vp)
{
  for (int i = 0; i < polys_.size(); ++i) {
    polys_[i] += vp[i];
  }
  return *this;
}

VectorPolynomial& VectorPolynomial::operator-=(const VectorPolynomial& vp)
{
  for (int i = 0; i < polys_.size(); ++i) {
    polys_[i] -= vp[i];
  }
  return *this;
}


/* ******************************************************************
* Set same origin for all polynomials without modyfying them 
****************************************************************** */
void VectorPolynomial::set_origin(const AmanziGeometry::Point& origin)
{
  for (int i = 0; i < size(); ++i) {
    polys_[i].set_origin(origin);
  }
}


/* ******************************************************************
* Change all polynomials to new same origin
****************************************************************** */
void VectorPolynomial::ChangeOrigin(const AmanziGeometry::Point& origin)
{
  for (int i = 0; i < size(); ++i) {
    polys_[i].ChangeOrigin(origin);
  }
}


/* ******************************************************************
* Ring algebra
****************************************************************** */
double VectorPolynomial::NormInf() const
{
  double tmp(0.0);
  for (int i = 0; i < polys_.size(); ++i) {
    tmp = std::max(tmp, polys_[i].NormInf());
  }
  return tmp;
}

}  // namespace WhetStone
}  // namespace Amanzi


