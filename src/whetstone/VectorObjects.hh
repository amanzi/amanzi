/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Vector of objects that admit a ring algebra. Current used values 
  for the template parameter are Polynomial and SpaceTimePolynomial.
  In addition to the conventional ring algebra operations, a few
  additional functions are need from the template class.
*/ 

#ifndef AMANZI_WHETSTONE_VECTOR_OBJECTS_HH_
#define AMANZI_WHETSTONE_VECTOR_OBJECTS_HH_

#include <vector>

#include "Point.hh"

#include "DenseVector.hh"
#include "Polynomial.hh"
#include "SpaceTimePolynomial.hh"
#include "Tensor.hh"
#include "VectorObjectsIterator.hh"

namespace Amanzi {
namespace WhetStone {

template<class T>
class VectorObjects {
 public:
  VectorObjects() : d_(0) {};
  VectorObjects(int d, int m) : d_(d) {
    polys_.resize(m);
    for (int i = 0; i < m; ++i) polys_[i].Reshape(d_, 0, true);
  }
  VectorObjects(int d, int m, int order) : d_(d) {
    polys_.resize(m);
    for (int i = 0; i < m; ++i) polys_[i].Reshape(d_, order, true);
  }
  VectorObjects(const T& p) : d_(p.dimension()) {
    polys_.resize(1);
    polys_[0] = p;
  }
  ~VectorObjects() {};

  // reshape polynomial and zero-out (optionally) memory
  void Reshape(int d, int m, int order, bool reset = false) {
    d_ = d;
    polys_.resize(m);
    for (int i = 0; i < m; ++i) polys_[i].Reshape(d, order, reset);
  }

  // minimal set of vector operations
  int NumRows() const { return polys_.size(); }
  void resize(int m) { polys_.resize(m); }

  T& operator[](int i) { return polys_[i]; }
  const T& operator[](int i) const { return polys_[i]; }

  // typical operations with polynomials
  void PutScalar(double val) {
    for (int i = 0; i < NumRows(); ++i) polys_[i].PutScalar(val);
  }
  double NormInf() const {
    double tmp(0.0);
    for (int i = 0; i < NumRows(); ++i) tmp = std::max(tmp, polys_[i].NormInf());
    return tmp;
  }

  // ring algebra
  VectorObjects<T>& operator*=(double val) {
    for (int i = 0; i < NumRows(); ++i) polys_[i] *= val;
    return *this;
  }

  VectorObjects<T>& operator+=(const VectorObjects<T>& vp) {
    for (int i = 0; i < NumRows(); ++i) polys_[i] += vp[i];
    return *this;
  }
  VectorObjects<T>& operator-=(const VectorObjects<T>& vp) {
    for (int i = 0; i < NumRows(); ++i) polys_[i] -= vp[i];
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

  // algebra of 2D and 3D vectors
  friend VectorObjects<T> operator^(const VectorObjects<T>& vp1, const VectorObjects<T>& vp2) {
    int d = vp1[0].dimension();

    AMANZI_ASSERT(d == 3);
    AMANZI_ASSERT(vp1.NumRows() == d);
    AMANZI_ASSERT(vp2.NumRows() == d);
    AMANZI_ASSERT(vp1[0].origin() == vp2[0].origin());

    VectorObjects<T> tmp(d, d, 0);
    tmp[0] = vp1[1] * vp2[2] - vp1[2] * vp2[1];
    tmp[1] = vp1[2] * vp2[0] - vp1[0] * vp2[2];
    tmp[2] = vp1[0] * vp2[1] - vp1[1] * vp2[0];
    tmp.set_origin(vp1[0].origin());

    return tmp;
  }

  friend VectorObjects<T> operator^(const VectorObjects<T>& vp1, const AmanziGeometry::Point& p2) {
    int d = vp1[0].dimension();

    AMANZI_ASSERT(d == 3);
    AMANZI_ASSERT(vp1.NumRows() == d);

    VectorObjects<T> tmp(d, d, 0);
    tmp[0] = vp1[1] * p2[2] - vp1[2] * p2[1];
    tmp[1] = vp1[2] * p2[0] - vp1[0] * p2[2];
    tmp[2] = vp1[0] * p2[1] - vp1[1] * p2[0];
    tmp.set_origin(vp1[0].origin());

    return tmp;
  }

  friend VectorObjects<T> operator^(const AmanziGeometry::Point& p1, const VectorObjects<T>& vp2) {
    return vp2 ^ p1;
  }

  // iterators require specialization
  VectorObjectsIterator<T> begin() {
    int order = (NumRows() > 0) ? polys_[0].order() : 0;
    auto it = VectorObjectsIterator<T>(d_, NumRows(), order);
    return it.begin();
  }

  VectorObjectsIterator<T> end() {
    int order = (NumRows() > 0) ? polys_[0].order() : 0;
    auto it = VectorObjectsIterator<T>(d_, NumRows(), order);
    return it.end();
  }

  // change of the coordinate system
  void set_origin(const AmanziGeometry::Point& origin) {
    for (int i = 0; i < NumRows(); ++i) polys_[i].set_origin(origin);
  }
  void ChangeOrigin(const AmanziGeometry::Point& origin) {
    for (int i = 0; i < NumRows(); ++i) polys_[i].ChangeOrigin(origin);
  }

  // typical operations with vector polynomials
  // -- value
  DenseVector Value(const AmanziGeometry::Point& xp) const {
    int m = NumRows();
    DenseVector val(m);
    for (int i = 0; i < m; ++i) val(i) = polys_[i].Value(xp);
    return val;
  }

  // generaization of the ring algebra
  // -- product vp * p 
  friend VectorObjects<T> operator*(const VectorObjects<T>& vp, const T& poly) {
    VectorObjects<T> tmp(vp);
    for (int i = 0; i < tmp.NumRows(); ++i) tmp[i] *= poly;
    return tmp;
  }

  // -- dot product vp * p
  friend T operator*(const VectorObjects<T>& vp, const AmanziGeometry::Point& p) {
    AMANZI_ASSERT(vp.NumRows() == p.dim());

    T tmp(vp[0] * p[0]);
    for (int i = 1; i < p.dim(); ++i) tmp += vp[i] * p[i];
    return tmp;
  }

  // -- product T * vp
  friend VectorObjects<T> operator*(const Tensor& K, const VectorObjects<T>& vp) {
    int d(K.dimension());
    VectorObjects<T> tmp(d, d, 0);

    if (K.rank() == 1) {
      tmp = vp;
      tmp *= K(0, 0);
    } else {
      tmp.set_origin(vp[0].origin());

      for (int i = 0; i < d; ++i) {
        for (int j = 0; j < d; ++j) {
          tmp[i] += K(i, j) * vp[j];
        }
      }
    }
    return tmp;
  }

  // -- dot product v1 * v2
  friend T operator*(const VectorObjects<T>& v1, const VectorObjects<T>& v2) {
    AMANZI_ASSERT(v1.NumRows() == v2.NumRows());

    T tmp(v1[0] * v2[0]);
    for (int i = 1; i < v1.NumRows(); ++i) tmp += v1[i] * v2[i];
    return tmp;
  }

  // output 
  friend std::ostream& operator << (std::ostream& os, const VectorObjects<T>& poly) {
    os << "Vector Object (length=" << poly.NumRows() << "):" << std::endl;
    for (int i = 0; i < poly.NumRows(); ++i) {
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

}  // namespace WhetStone
}  // namespace Amanzi

#endif

