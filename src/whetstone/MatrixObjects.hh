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

  Operations with matrix of objects that admit a ring algebra.
*/

#ifndef AMANZI_WHETSTONE_MATRIX_OBJECTS_HH_
#define AMANZI_WHETSTONE_MATRIX_OBJECTS_HH_

#include <vector>

#include "Point.hh"

#include "DenseMatrix.hh"
#include "VectorObjects.hh"

namespace Amanzi {
namespace WhetStone {

template <class T>
class MatrixObjects {
 public:
  MatrixObjects() : d_(0), m_(0), n_(0){};
  MatrixObjects(int d, int m, int n, int order) : d_(d), m_(m), n_(n), order_(order)
  {
    polys_.resize(m_);
    for (int i = 0; i < m_; ++i) {
      polys_[i].resize(n_);
      for (int j = 0; j < n_; ++j) polys_[i][j].Reshape(d_, order, true);
    }
  }
  ~MatrixObjects(){};

  // reshape polynomial with erase (optionally) memory
  void Reshape(int d, int m, int n, int order, bool reset = false);

  // minimal set of vector operations
  int NumRows() const { return m_; }
  int NumCols() const { return n_; }

  T& operator()(int i, int j) { return polys_[i][j]; }
  const T& operator()(int i, int j) const { return polys_[i][j]; }

  // typical operations with polynomials
  void PutScalar(double val)
  {
    for (int i = 0; i < m_; ++i)
      for (int j = 0; j < n_; ++j) polys_[i][j].PutScalar(val);
  }
  double NormInf() const;

  // ring algebra
  template <typename U>
  MatrixObjects& operator*=(const U val)
  {
    for (int i = 0; i < m_; ++i)
      for (int j = 0; j < n_; ++j) polys_[i][j] *= val;
    return *this;
  }

  MatrixObjects<T>& operator+=(const MatrixObjects<T>& mp);
  MatrixObjects<T>& operator-=(const MatrixObjects<T>& mp);

  friend MatrixObjects<T> operator+(const MatrixObjects<T>& m1, const MatrixObjects<T>& m2)
  {
    MatrixObjects<T> tmp(m1);
    return tmp += m2;
  }

  friend MatrixObjects<T> operator-(const MatrixObjects<T>& m1, const MatrixObjects<T>& m2)
  {
    MatrixObjects<T> tmp(m1);
    return tmp -= m2;
  }

  friend MatrixObjects<T> operator*(double val, const MatrixObjects<T>& m)
  {
    MatrixObjects<T> tmp(m);
    return tmp *= val;
  }

  friend MatrixObjects<T> operator*(const MatrixObjects<T>& m, double val)
  {
    MatrixObjects tmp(m);
    return tmp *= val;
  }

  friend VectorObjects<T> operator*(const MatrixObjects<T>& m, const AmanziGeometry::Point& p)
  {
    VectorObjects<T> v;
    m.Multiply(p, v, false);
    return v;
  }

  // change the coordinate system
  // change the coordinate system
  void set_origin(const AmanziGeometry::Point& origin);
  void ChangeOrigin(const AmanziGeometry::Point& origin);

  // simple operations with vector polynomials
  // -- value
  DenseMatrix Value(const AmanziGeometry::Point& xp) const;

  // -- matrix-vector products
  void Multiply(const VectorObjects<T>& v, VectorObjects<T>& av, bool transpose);
  void Multiply(const DenseVector& v, VectorObjects<T>& av, bool transpose);
  void Multiply(const AmanziGeometry::Point& p, VectorObjects<T>& av, bool transpose) const;

  // output
  friend std::ostream& operator<<(std::ostream& os, const MatrixObjects<T>& poly)
  {
    os << "Matrix of Objects: " << poly.NumRows() << " x " << poly.NumCols() << std::endl;
    for (int i = 0; i < poly.NumRows(); ++i)
      for (int j = 0; j < poly.NumCols(); ++j) os << "i=" << i << ", j=" << j << " " << poly(i, j);
    return os;
  }

 private:
  int d_, m_, n_, order_;
  std::vector<std::vector<T>> polys_;
};


// used types
typedef MatrixObjects<Polynomial> MatrixPolynomial;
typedef MatrixObjects<SpaceTimePolynomial> MatrixSpaceTimePolynomial;


/* ******************************************************************
* Re-shape polynomials
****************************************************************** */
template <class T>
void
MatrixObjects<T>::Reshape(int d, int m, int n, int order, bool reset)
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
* Calculate value at a point
****************************************************************** */
template <class T>
DenseMatrix
MatrixObjects<T>::Value(const AmanziGeometry::Point& xp) const
{
  DenseMatrix val(m_, n_);

  for (int i = 0; i < m_; ++i)
    for (int j = 0; j < n_; ++j) val(i, j) = polys_[i][j].Value(xp);

  return val;
}


/* ******************************************************************
* Ring algebra
****************************************************************** */
template <class T>
MatrixObjects<T>&
MatrixObjects<T>::operator+=(const MatrixObjects<T>& mp)
{
  for (int i = 0; i < m_; ++i)
    for (int j = 0; j < n_; ++j) polys_[i][j] += mp(i, j);

  return *this;
}

template <class T>
MatrixObjects<T>&
MatrixObjects<T>::operator-=(const MatrixObjects<T>& mp)
{
  for (int i = 0; i < m_; ++i)
    for (int j = 0; j < n_; ++j) polys_[i][j] -= mp(i, j);

  return *this;
}


/* ******************************************************************
* Matrix-vector operations
***************************************************************** */
template <class T>
void
MatrixObjects<T>::Multiply(const VectorObjects<T>& v, VectorObjects<T>& av, bool transpose)
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


template <class T>
void
MatrixObjects<T>::Multiply(const DenseVector& v, VectorObjects<T>& av, bool transpose)
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


template <class T>
void
MatrixObjects<T>::Multiply(const AmanziGeometry::Point& p,
                           VectorObjects<T>& av,
                           bool transpose) const
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
* Set same origin for all polynomials without modyfying them
****************************************************************** */
template <class T>
void
MatrixObjects<T>::set_origin(const AmanziGeometry::Point& origin)
{
  for (int i = 0; i < m_; ++i)
    for (int j = 0; j < n_; ++j) polys_[i][j].set_origin(origin);
}


/* ******************************************************************
* Change all polynomials to new same origin
****************************************************************** */
template <class T>
void
MatrixObjects<T>::ChangeOrigin(const AmanziGeometry::Point& origin)
{
  for (int i = 0; i < m_; ++i)
    for (int j = 0; j < n_; ++j) polys_[i][j].ChangeOrigin(origin);
}


/* ******************************************************************
* Maximum norm
****************************************************************** */
template <class T>
double
MatrixObjects<T>::NormInf() const
{
  double tmp(0.0);
  for (int i = 0; i < m_; ++i)
    for (int j = 0; j < n_; ++j) tmp = std::max(tmp, polys_[i][j].NormInf());

  return tmp;
}


// non-member functions
// -- gradient
template <class T>
MatrixObjects<T>
Gradient(const VectorObjects<T>& p)
{
  int d = p[0].dimension(), n = p.size();
  MatrixObjects<T> grad(d, d, n, 0);
  for (int i = 0; i < n; ++i) {
    auto tmp = Gradient(p[i]);
    for (int k = 0; k < d; ++k) grad(i, k) = tmp[k];
  }
  return grad;
}

} // namespace WhetStone
} // namespace Amanzi

#endif
