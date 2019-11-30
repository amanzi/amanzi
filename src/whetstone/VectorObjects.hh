/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Vector of objects that admit a ring algebra. Possible values 
  for the template parameter are Monomial, Polynomial and
  SpaceTimePolynomial.
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

  // reshape polynomial with erase (optionally) memory
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

  // change the coordinate system
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

  // dot product v1 * p 
  friend T operator*(const VectorObjects<T>& poly, const AmanziGeometry::Point& p) {
    AMANZI_ASSERT(poly.NumRows() == p.dim());

    T tmp(poly[0] * p[0]);

    for (int i = 1; i < p.dim(); ++i) {
      tmp += poly[i] * p[i];
    }
    return tmp;
  }

  // dot product T * v1 
  friend VectorObjects<T> operator*(const Tensor& K, const VectorObjects<T>& poly) {
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
    AMANZI_ASSERT(v1.NumRows() == v2.NumRows());

    T tmp(v1[0] * v2[0]);

    for (int i = 1; i < v1.NumRows(); ++i) {
      tmp += v1[i] * v2[i];
    }
    return tmp;
  }

  // vector decompositions
  // -- curl-based
  //  3D: q_k = curl(p_k ^ x) + x . p_{k-1}
  //  2D: q_k =  rot(p_{k+1}) + x . p_{k-1}
  void VectorDecomposition3DCurl(VectorObjects<T>& p1, T& p2);
  void VectorDecomposition2DRot(T& p1, T& p2);
  // -- grad based
  //  3D: q_k = grad(p_{k+1}) + x ^ p_{k-1}
  //  2D: q_k = grad(p_{k+1}) + x*. p_{k-1}, x* = (-y, x)
  void VectorDecomposition3DGrad(T& p1, VectorObjects<T>& p2) {};
  void VectorDecomposition2DGrad(T& p1, T& p2) {};

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


// non-member functions
// -- differential operators
VectorPolynomial Gradient(const Polynomial& p);
VectorSpaceTimePolynomial Gradient(const SpaceTimePolynomial& p);

VectorPolynomial Curl3D(const VectorPolynomial& p);
Polynomial Curl2D(const VectorPolynomial& p);
VectorPolynomial Rot2D(const Polynomial& p);

Polynomial Divergence(const VectorObjects<Polynomial>& vp);

// -- algebra
VectorPolynomial operator^(const VectorPolynomial& p1, const VectorPolynomial& p2);

// project gradient of the given polynomial on unit sphere using
// the Taylor expansion with k terms
VectorPolynomial GradientOnUnitSphere(const Polynomial& poly, int k);


// Specializations

/* ******************************************************************
* 3D vector decomposition: q_k = curl(p_k ^ x) + x . p_{k-1}
****************************************************************** */
template<>
inline
void VectorPolynomial::VectorDecomposition3DCurl(VectorPolynomial& p1, Polynomial& p2)
{
  // reshape output
  p1 = *this;

  int d = p1[0].dimension();
  int order(0);
  for (int k = 0; k < d; ++k) order = std::max(order, p1[k].order() - 1);
  p2.Reshape(d, order, true);
  p2.set_origin(p1[0].origin());

  // calculate decomposition for each monomial of each component
  int idx[3];
  for (int k = 0; k < d; ++k) {
    for (auto it = p1[k].begin(); it < p1[k].end(); ++it) {
      int n = it.PolynomialPosition();
      const int* index = it.multi_index();

      double a(2.0);
      for (int i = 0; i < d; ++i) a += index[i];

      double coef = p1[k](n);
      p1[k](n) = coef / a;

      if (index[k] > 0) {
        for (int i = 0; i < d; ++i) idx[i] = index[i];
        idx[k]--;
        
        int l = PolynomialPosition(d, idx);
        p2(l) += coef * index[k] / a;
      }
    }
  }
}


/* ******************************************************************
* 2D vector decomposition: q_k = tor(p_{k+1}) + x . p_{k-1}
****************************************************************** */
template<>
inline
void VectorPolynomial::VectorDecomposition2DRot(Polynomial& p1, Polynomial& p2)
{
  // reshape output
  int d = polys_[0].dimension();
  int order(0);
  for (int k = 0; k < d; ++k) order = std::max(order, polys_[k].order() - 1);

  p1.Reshape(d, order + 2, true);
  p1.set_origin(polys_[0].origin());

  p2.Reshape(d, order, true);
  p2.set_origin(polys_[0].origin());

  // calculate decomposition for each monomial of each component
  int idx[3];
  for (int k = 0; k < d; ++k) {
    for (auto it = polys_[k].begin(); it < polys_[k].end(); ++it) {
      int n = it.PolynomialPosition();
      const int* index = it.multi_index();

      double a(1.0);
      for (int i = 0; i < d; ++i) {
        a += index[i];
        idx[i] = index[i];
      }

      idx[1 - k]++;
      int l = PolynomialPosition(d, idx);

      double coef = polys_[k](n);
      if (k == 0)
        p1(l) += coef / a;
      else
        p1(l) -= coef / a;
      
      if (index[k] > 0) {
        idx[1 - k]--;
        idx[k]--;
        
        int l = PolynomialPosition(d, idx);
        p2(l) += coef * index[k] / a;
      }
    }
  }
}

}  // namespace WhetStone
}  // namespace Amanzi

#endif

