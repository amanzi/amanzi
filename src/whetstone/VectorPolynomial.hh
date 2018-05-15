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
#include "Polynomial.hh"

namespace Amanzi {
namespace WhetStone {

class VectorPolynomial {
 public:
  VectorPolynomial() : d_(0) {};
  VectorPolynomial(int d, int size);
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
  VectorPolynomial& operator*=(double val);
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

  friend VectorPolynomial operator*(double val, const VectorPolynomial& vp) {
    VectorPolynomial tmp(vp);
    return tmp *= val;
  }

  // complex constructions
  // -- gradient of a polynomial
  void Gradient(const Polynomial p);

  // -- matrix-vector product A * v
  void Multiply(const std::vector<std::vector<Polynomial> >& A,
                const VectorPolynomial& v, bool transpose);

  // -- matrix-point product A * p
  void Multiply(const std::vector<std::vector<Polynomial> >& A,
                const AmanziGeometry::Point& p, bool transpose);

  // dot product 
  friend Polynomial operator*(const VectorPolynomial& poly, const AmanziGeometry::Point& p) {
    int d(p.dim());
    Polynomial tmp(d, 0);
    tmp.set_origin(poly[0].origin());

    for (int i = 0; i < d; ++i) {
      tmp += poly[i] * p[i];
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

 private:
  int d_;
  std::vector<Polynomial> polys_;
};

}  // namespace WhetStone
}  // namespace Amanzi

#endif

