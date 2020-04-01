/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Operations with vector of objects that admit ring algebra.
  Avalaible values for template are Polynomial and SpaceTimePolynomial.
*/ 

#ifndef AMANZI_WHETSTONE_VECTOR_OBJECTS_HH_
#define AMANZI_WHETSTONE_VECTOR_OBJECTS_HH_

#include <vector>

#include "Point.hh"

#include "DenseVector.hh"
#include "Polynomial.hh"
#include "SpaceTimePolynomial.hh"
#include "Tensor.hh"

namespace Amanzi {
namespace WhetStone {

template<class T>
class VectorObjects {
 public:
  VectorObjects() : d_(0) {};
  VectorObjects(int d, int size) : d_(d) {
    polys_.resize(size);
    for (int i = 0; i < size; ++i) polys_[i].reshape(d_, 0, true);
  }
  VectorObjects(int d, int size, int order) : d_(d) {
    polys_.resize(size);
    for (int i = 0; i < size; ++i) polys_[i].reshape(d_, order, true);
  }
  VectorObjects(const T& p) : d_(p.dimension()) {
    polys_.resize(1);
    polys_[0] = p;
  }
  ~VectorObjects() {};

  // reshape polynomial with erase (optionally) memory
  void reshape(int d, int m, int order, bool reset = false) {
    d_ = d;
    polys_.resize(m);
    for (int i = 0; i < m; ++i) polys_[i].reshape(d, order, reset);
  }

  // minimal set of vector operations
  int size() const { return polys_.size(); }
  void resize(int size) { polys_.resize(size); }

  T& operator[](int i) { return polys_[i]; }
  const T& operator[](int i) const { return polys_[i]; }

  // typical operations with polynomials
  void putScalar(double val) {
    for (int i = 0; i < size(); ++i) polys_[i].putScalar(val);
  }
  double NormInf() const {
    double tmp(0.0);
    for (int i = 0; i <size(); ++i) tmp = std::max(tmp, polys_[i].NormInf());
    return tmp;
  }

  // ring algebra
  VectorObjects<T>& operator*=(double val) {
    for (int i = 0; i < polys_.size(); ++i) polys_[i] *= val;
    return *this;
  }

  VectorObjects<T>& operator+=(const VectorObjects<T>& vp) {
    for (int i = 0; i < polys_.size(); ++i) polys_[i] += vp[i];
    return *this;
  }
  VectorObjects<T>& operator-=(const VectorObjects<T>& vp) {
    for (int i = 0; i < polys_.size(); ++i) polys_[i] -= vp[i];
    return *this;
  }

  friend VectorObjects<T> operator+(const VectorObjects<T>& vp1, const VectorObjects<T>& vp2) {
    VectorObjects<T> tmp(vp1);
    return tmp += vp2;
  }

  friend VectorObjects<T> operator-(const VectorObjects<T>& vp1, const VectorObjects<T>& vp2) {
    VectorObjects<T> tmp(vp1);
    return tmp -= vp2;
  }

  friend VectorObjects<T> operator*(double val, const VectorObjects<T>& vp) {
    VectorObjects<T> tmp(vp);
    return tmp *= val;
  }

  friend VectorObjects<T> operator*(const VectorObjects<T>& vp, double val) {
    VectorObjects<T> tmp(vp);
    return tmp *= val;
  }

  // change the coordinate system
  void set_origin(const AmanziGeometry::Point& origin) {
    for (int i = 0; i < size(); ++i) polys_[i].set_origin(origin);
  }
  void ChangeOrigin(const AmanziGeometry::Point& origin) {
    for (int i = 0; i < size(); ++i) polys_[i].ChangeOrigin(origin);
  }

  // typical operations with vector polynomials
  // -- value
  DenseVector Value(const AmanziGeometry::Point& xp) const {
    int n = polys_.size();
    DenseVector val(n);
    for (int i = 0; i < n; ++i) val(i) = polys_[i].Value(xp);
    return val;
  }

  // dot product v1 * p 
  friend T operator*(const VectorObjects<T>& poly, const AmanziGeometry::Point& p) {
    AMANZI_ASSERT(poly.size() == p.dim());

    T tmp(poly[0] * p[0]);

    for (int i = 1; i < p.dim(); ++i) {
      tmp += poly[i] * p[i];
    }
    return tmp;
  }

  // dot product T * v1 
  template<class MEMSPACE>
  friend VectorObjects<T> operator*(const Tensor<MEMSPACE>& K, const VectorObjects<T>& poly) {
    int d(K.dimension());
    VectorObjects<T> tmp(d, d, 0);

    if (K.rank() == 1) {
      tmp = poly;
      tmp *= K(0, 0);
    } else {
      tmp.set_origin(poly[0].origin());

      for (int i = 0; i < d; ++i) {
        for (int j = 0; j < d; ++j) {
          tmp[i] += K(i, j) * poly[j];
        }
      }
    }
    return tmp;
  }

  // dot product v1 * v2
  friend T operator*(const VectorObjects<T>& v1, const VectorObjects<T>& v2) {
    AMANZI_ASSERT(v1.size() == v2.size());

    T tmp(v1[0] * v2[0]);

    for (int i = 1; i < v1.size(); ++i) {
      tmp += v1[i] * v2[i];
    }
    return tmp;
  }

  // output 
  friend std::ostream& operator << (std::ostream& os, const VectorObjects<T>& poly) {
    os << "Vector Object (length=" << poly.size() << "):" << std::endl;
    for (int i = 0; i < poly.size(); ++i) {
      os << "i=" << i << " " << poly[i];
    }
    return os;
  }

 private:
  int d_;
  std::vector<T> polys_;
};


// used types
typedef VectorObjects<Polynomial> VectorPolynomial;
typedef VectorObjects<SpaceTimePolynomial> VectorSpaceTimePolynomial;


// non-member functions
VectorPolynomial Gradient(const Polynomial& p);
VectorSpaceTimePolynomial Gradient(const SpaceTimePolynomial& p);

Polynomial Divergence(const VectorObjects<Polynomial>& vp);
VectorPolynomial GradientOnUnitSphere(const Polynomial& poly, int k);

// project gradient of the given polynomial on unit sphere using
// the Taylor expansion with k terms
VectorPolynomial GradientOnUnitSphere(const Polynomial& poly, int k);

}  // namespace WhetStone
}  // namespace Amanzi

#endif

