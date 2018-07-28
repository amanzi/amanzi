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

#ifndef AMANZI_WHETSTONE_VECTORPOLYNOMIAL_HH_
#define AMANZI_WHETSTONE_VECTORPOLYNOMIAL_HH_

#include <vector>

#include "Point.hh"

#include "DenseVector.hh"
#include "Polynomial.hh"
#include "Tensor.hh"

namespace Amanzi {
namespace WhetStone {

class VectorPolynomial {
 public:
  VectorPolynomial() : d_(0) {};
  VectorPolynomial(int d, int size);
  VectorPolynomial(int d, int size, int order);
  VectorPolynomial(const Polynomial& p);
  ~VectorPolynomial() {};

  // minimal set of vector operations
  int size() const { return polys_.size(); }
  void resize(int size) { polys_.resize(size); }

  Polynomial& operator[](int i) { return polys_[i]; }
  const Polynomial& operator[](int i) const { return polys_[i]; }

  // typical operations with polynomials
  void PutScalar(double val);
  double NormMax() const;

  // ring algebra
  template<typename Type>
  VectorPolynomial& operator*=(Type val) {
    for (int i = 0; i < polys_.size(); ++i) {
      polys_[i] *= val;
    }
    return *this;
  }

  VectorPolynomial& operator+=(const VectorPolynomial& vp);
  VectorPolynomial& operator-=(const VectorPolynomial& vp);

  friend VectorPolynomial operator+(const VectorPolynomial& vp1, const VectorPolynomial& vp2) {
    VectorPolynomial tmp(vp1);
    return tmp += vp2;
  }

  friend VectorPolynomial operator-(const VectorPolynomial& vp1, const VectorPolynomial& vp2) {
    VectorPolynomial tmp(vp1);
    return tmp -= vp2;
  }

  template<typename Type>
  friend VectorPolynomial operator*(const Type& val, const VectorPolynomial& vp) {
    VectorPolynomial tmp(vp);
    return tmp *= val;
  }

  template<typename Type>
  friend VectorPolynomial operator*(const VectorPolynomial& vp, const Type& val) {
    VectorPolynomial tmp(vp);
    return tmp *= val;
  }

  // change the coordinate system
  void set_origin(const AmanziGeometry::Point& origin);
  void ChangeOrigin(const AmanziGeometry::Point& origin);

  // typical operations with vector polynomials
  // -- value
  DenseVector Value(const AmanziGeometry::Point& xp) const;

  // -- gradient of a polynomial
  void Gradient(const Polynomial p);

  // -- matrix-vector product A * v
  void Multiply(const std::vector<std::vector<Polynomial> >& A,
                const VectorPolynomial& v, bool transpose);

  // -- matrix-point product A * p
  void Multiply(const std::vector<std::vector<Polynomial> >& A,
                const AmanziGeometry::Point& p, bool transpose);

  // dot product v1 * p 
  friend Polynomial operator*(const VectorPolynomial& poly, const AmanziGeometry::Point& p) {
    int d(p.dim());
    Polynomial tmp(d, 0);
    tmp.set_origin(poly[0].origin());

    for (int i = 0; i < d; ++i) {
      tmp += poly[i] * p[i];
    }
    return tmp;
  }

  // dot product T * v1 
  friend VectorPolynomial operator*(const Tensor& T, const VectorPolynomial& poly) {
    int d(T.dimension());
    VectorPolynomial tmp(d, d, 0);

    if (T.rank() == 1) {
      tmp = poly;
      tmp *= T(0, 0);
    } else {
      tmp.set_origin(poly[0].origin());

      for (int i = 0; i < d; ++i) {
        for (int j = 0; j < d; ++j) {
          tmp[i] += T(i, j) * poly[j];
        }
      }
    }
    return tmp;
  }

  // dot product v1 * v2
  friend Polynomial operator*(const VectorPolynomial& v1, const VectorPolynomial& v2) {
    AMANZI_ASSERT(v1.size() == v2.size());

    int d(v1[0].dimension());
    Polynomial tmp(d, 0);
    tmp.set_origin(v1[0].origin());

    for (int i = 0; i < v1.size(); ++i) {
      tmp += v1[i] * v2[i];
    }
    return tmp;
  }

  // output 
  friend std::ostream& operator << (std::ostream& os, const VectorPolynomial& poly) {
    os << "Vector Polynomial (size=" << poly.size() << "):" << std::endl;
    for (int i = 0; i < poly.size(); ++i) {
      os << "i=" << i << " " << poly[i];
    }
    return os;
  }

 private:
  int d_;
  std::vector<Polynomial> polys_;
};

// non-member functions
// --divergence
inline
Polynomial Divergence(const VectorPolynomial vp) 
{
  int d = vp[0].dimension();
  AMANZI_ASSERT(d == vp.size());

  int order = vp[0].order();
  order = std::max(0, order - 1);

  Polynomial div(d, order);
  div.set_origin(vp[0].origin());

  int index[3];
  for (int i = 0; i < d; ++i) {
    for (auto it = vp[i].begin(); it < vp[i].end(); ++it) {
      int k = it.MonomialSetOrder();
      if (k > 0) {
        const int* idx = it.multi_index();
        for (int j = 0; j < d; ++j) index[j] = idx[j];

        if (index[i] > 0) {
          int m = it.MonomialSetPosition();
          double val = vp[i](k, m);

          index[i]--;
          m = vp[i].MonomialSetPosition(index);
          div(k - 1, m) += val * idx[i];
        }
      }
    }
  }

  return div;
}

}  // namespace WhetStone
}  // namespace Amanzi

#endif

