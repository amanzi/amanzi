/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Operations with matrix of polynomials of type p(x - x0) where x0
  could be different for each entry.
*/

#ifndef AMANZI_WHETSTONE_MATRIX_POLYNOMIAL_HH_
#define AMANZI_WHETSTONE_MATRIX_POLYNOMIAL_HH_

#include <vector>

#include "Point.hh"

#include "DenseMatrix.hh"
#include "Polynomial.hh"
#include "VectorPolynomial.hh"
#include "Tensor.hh"

namespace Amanzi {
namespace WhetStone {

class MatrixPolynomial {
 public:
  MatrixPolynomial() : d_(0), m_(0), n_(0){};
  MatrixPolynomial(int d, int m, int n, int order);
  ~MatrixPolynomial(){};

  // reshape polynomial with erase (optionally) memory
  void Reshape(int d, int m, int n, int order, bool reset = false);

  // minimal set of vector operations
  int NumRows() const { return m_; }
  int NumCols() const { return n_; }

  Polynomial& operator()(int i, int j) { return polys_[i][j]; }
  const Polynomial& operator()(int i, int j) const { return polys_[i][j]; }

  // typical operations with polynomials
  void PutScalar(double val);
  double NormInf() const;

  // ring algebra
  template <typename Type>
  MatrixPolynomial& operator*=(Type val)
  {
    for (int i = 0; i < m_; ++i)
      for (int j = 0; j < n_; ++j) polys_[i][j] *= val;
    return *this;
  }

  MatrixPolynomial& operator+=(const MatrixPolynomial& mp);
  MatrixPolynomial& operator-=(const MatrixPolynomial& mp);

  friend MatrixPolynomial
  operator+(const MatrixPolynomial& vp1, const MatrixPolynomial& vp2)
  {
    MatrixPolynomial tmp(vp1);
    return tmp += vp2;
  }

  friend MatrixPolynomial
  operator-(const MatrixPolynomial& vp1, const MatrixPolynomial& vp2)
  {
    MatrixPolynomial tmp(vp1);
    return tmp -= vp2;
  }

  template <typename Type>
  friend MatrixPolynomial operator*(const Type& val, const MatrixPolynomial& vp)
  {
    MatrixPolynomial tmp(vp);
    return tmp *= val;
  }

  template <typename Type>
  friend MatrixPolynomial operator*(const MatrixPolynomial& vp, const Type& val)
  {
    MatrixPolynomial tmp(vp);
    return tmp *= val;
  }

  // change the coordinate system
  void set_origin(const AmanziGeometry::Point& origin);
  void ChangeOrigin(const AmanziGeometry::Point& origin);

  // typical operations with vector polynomials
  // -- value
  DenseMatrix Value(const AmanziGeometry::Point& xp) const;

  // -- matrix-vector products
  void
  Multiply(const VectorPolynomial& v, VectorPolynomial& av, bool transpose);
  void Multiply(const DenseVector& v, VectorPolynomial& av, bool transpose);
  void Multiply(const AmanziGeometry::Point& p, VectorPolynomial& av,
                bool transpose);

  // output
  friend std::ostream&
  operator<<(std::ostream& os, const MatrixPolynomial& poly)
  {
    os << "Matrix Polynomial: " << poly.NumRows() << " x " << poly.NumCols()
       << std::endl;
    for (int i = 0; i < poly.NumRows(); ++i)
      for (int j = 0; j < poly.NumCols(); ++j)
        os << "i=" << i << ", j=" << j << " " << poly(i, j);
    return os;
  }

 private:
  int d_, m_, n_, order_;
  std::vector<std::vector<Polynomial>> polys_;
};

} // namespace WhetStone
} // namespace Amanzi

#endif
