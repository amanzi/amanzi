/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Operations with polynomials.
*/

#include <cmath>

#include "DenseMatrix.hh"
#include "DenseVector.hh"
#include "Monomial.hh"
#include "Polynomial.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Constructor of zero polynomial.
****************************************************************** */
Polynomial::Polynomial(int d, int order) : PolynomialBase(d, order)
{
  size_ = PolynomialSpaceDimension(d_, order_);
  coefs_.Reshape(size_);
  coefs_.PutScalar(0.0);
}


/* ******************************************************************
* Constructor from a given vector.
****************************************************************** */
Polynomial::Polynomial(int d, int order, const DenseVector& coefs)
  : PolynomialBase(d, order)
{
  size_ = PolynomialSpaceDimension(d_, order_);
  AMANZI_ASSERT(size_ == coefs.NumRows());
  coefs_ = coefs;
}


/* ******************************************************************
* Constructor of a polynomial with a single term:
*    p(x) = factor * (x)^multi_index
****************************************************************** */
Polynomial::Polynomial(int d, const int* multi_index, double factor)
   : PolynomialBase(d, 0)
{
  for (int i = 0; i < d_; ++i) order_ += multi_index[i];

  size_ = PolynomialSpaceDimension(d_, order_);
  coefs_.Reshape(size_);
  coefs_.PutScalar(0.0);

  int l = PolynomialPosition(d_, multi_index);
  coefs_(l) = factor;
}


/* ******************************************************************
* Copy constructor.
****************************************************************** */
Polynomial::Polynomial(const Polynomial& poly)
{
  d_ = poly.dimension();
  order_ = poly.order();
  origin_ = poly.origin();
  size_ = poly.size();
  coefs_ = poly.coefs();
}


/* ******************************************************************
* Constructor of a polynomial from a monomial.
****************************************************************** */
Polynomial::Polynomial(const Monomial& mono)
{
  d_ = mono.dimension();
  order_ = mono.order();
  origin_ = mono.origin();

  size_ = PolynomialSpaceDimension(d_, order_);
  coefs_.Reshape(size_);
  coefs_.PutScalar(0.0);

  int l = PolynomialPosition(d_, mono.multi_index());
  coefs_(l) = mono.coefs()(0);
}


/* ******************************************************************
* Complex constructor based on minimum set of given points and value.
****************************************************************** */
Polynomial::Polynomial(int d, int order,
                       const std::vector<AmanziGeometry::Point>& xyz, 
                       const DenseVector& values)
  : PolynomialBase(d, order)
{
  d_ = d;
  order_ = order;
  size_ = PolynomialSpaceDimension(d, order);

  AMANZI_ASSERT(size_ == xyz.size());
  AMANZI_ASSERT(size_ == values.NumRows());

  coefs_.Reshape(size_);

  if (order == 0) {
    coefs_(0) = values(0);
  } else {
    // evaluate basis functions at given points
    DenseMatrix psi(size_, size_);

    for (auto it = begin(); it < end(); ++it) {
      int i = it.PolynomialPosition();
      const int* idx = it.multi_index();

      for (int n = 0; n < size_; ++n) {
        double val(1.0);
        for (int k = 0; k < d_; ++k) {
          val *= std::pow(xyz[n][k], idx[k]);
        }
        psi(n, i) = val;
      }
    }
      
    // form linear system
    DenseMatrix A(size_, size_);
    DenseVector b(size_);

    A.Multiply(psi, psi, true);
    psi.Multiply(values, b, true);

    // solver linear systems
    A.Inverse();
    A.Multiply(b, coefs_, false);
  }
}


/* ******************************************************************
* Re-shape polynomial
* NOTE: case d > d_ can be treated more intelligently.
****************************************************************** */
void Polynomial::Reshape(int d, int order, bool reset)
{
  if (d_ != d) {
    d_ = d;
    order_ = order;
    origin_ = AmanziGeometry::Point(d);

    size_ = PolynomialSpaceDimension(d_, order_);
    coefs_.Reshape(size_);
    coefs_.PutScalar(0.0);
  } else if (order_ != order) {
    int size = size_;

    order_ = order;
    size_ = PolynomialSpaceDimension(d_, order_);
    coefs_.Reshape(size_);

    if (reset) { 
      coefs_.PutScalar(0.0);
    } else {
      double* data = coefs_.Values();
      for (int i = size; i < size_; ++i) data[i] = 0.0;
    }
  } else if (reset) {
    coefs_.PutScalar(0.0);
  }
}


/* ******************************************************************
* Implemented ring algebra operations.
* NOTE: implementation is order independent.
****************************************************************** */
Polynomial& Polynomial::operator+=(const Polynomial& poly)
{
  AMANZI_ASSERT(d_ == poly.dimension());  // FIXME
  AMANZI_ASSERT(origin_ == poly.origin());

  int order = poly.order();
  if (order_ < order) Reshape(d_, order);
  for (int i = 0; i < poly.size(); ++i) coefs_(i) += poly(i);

  return *this;
}


Polynomial& Polynomial::operator-=(const Polynomial& poly)
{
  AMANZI_ASSERT(d_ == poly.dimension());  // FIXME
  AMANZI_ASSERT(origin_ == poly.origin());

  int order = poly.order();
  if (order_ < order) Reshape(d_, order);
  for (int i = 0; i < poly.size(); ++i) coefs_(i) -= poly(i);

  return *this;
}


Polynomial& Polynomial::operator*=(const Polynomial& poly)
{
  Polynomial tmp(*this);
  *this = tmp * poly;
  return *this;
}


Polynomial& Polynomial::operator*=(double val) {
  coefs_ *= val;
  return *this;
}


/* ******************************************************************
* Rebase polynomial to different origin.
****************************************************************** */
void Polynomial::ChangeOrigin(const AmanziGeometry::Point& origin)
{
  AmanziGeometry::Point shift(origin - origin_);

  if (order_ == 1) {
    for (int i = 0; i < d_; ++i) {
      coefs_(0) += coefs_(i + 1) * shift[i];
    }
  } else if (order_ > 1) {
    if (AmanziGeometry::L22(shift) == 0.0) return;

    // create powers (x_i - o_i)^k
    // FIXME: cost could be reduced using split Cnk * a^k
    std::vector<std::vector<DenseVector> > powers(d_);

    for (int i = 0; i < d_; ++i) {
      powers[i].resize(order_ + 1);

      for (int k = 0; k <= order_; ++k) {
        DenseVector& dv = powers[i][k];
        dv.Reshape(k + 1);

        int cnk(1);
        double val(1.0), a(shift[i]);
        for (int n = 0; n <= k; ++n) {
          dv(k - n) = val * cnk;
          cnk *= k - n;
          cnk /= n + 1;
          val *= a;
        }
      }
    }

    // iterate over polynomial and sum up products
    Polynomial rebased(d_, order_);
    for (auto it = begin(); it < end(); ++it) {
      int k = it.MonomialSetOrder();
      int m = it.PolynomialPosition();
      double coef = coefs_(m);
      if (coef == 0.0) continue;

      const int* index = it.multi_index();
      // product of (x-a)^i (y-b)^j (z-c)^k
      int idx[3];
      Polynomial tmp(d_, k);
      for (int i0 = 0; i0 <= index[0]; ++i0) {
        idx[0] = i0;
        double coef0 = powers[0][index[0]](i0); 

        for (int i1 = 0; i1 <= index[1]; ++i1) {
          idx[1] = i1;
          double coef1 = powers[1][index[1]](i1); 

          if (d_ == 2) {
            int pos = PolynomialPosition(d_, idx);
            tmp(pos) = coef * coef0 * coef1;
          } else {
            for (int i2 = 0; i2 <= index[2]; ++i2) {
              idx[2] = i2;
              double coef2 = powers[2][index[2]](i2); 

              int pos = PolynomialPosition(d_, idx);
              tmp(pos) = coef * coef0 * coef1 * coef2;
            }
          }
        }
      }
      rebased += tmp;
    }
    
    *this = rebased;
  }
  origin_ = origin;   
}


/* ******************************************************************
* Rebase monomial to different origin. 
****************************************************************** */
Polynomial Polynomial::ChangeOrigin(
    const Monomial& mono, const AmanziGeometry::Point& origin)
{
  int d = mono.dimension();
  int order = mono.order();
  double coef = mono.coefs()(0);

  Polynomial poly(d, order);

  if (order == 0) {
    poly(0) = coef;
  }
  else if (order > 0) {
    AmanziGeometry::Point shift(origin - mono.origin());
    const int* index = mono.multi_index();

    // create powers (x_i - o_i)^k
    std::vector<DenseVector> powers(d);

    for (int i = 0; i < d; ++i) {
      powers[i].Reshape(index[i] + 1);

      int cnk(1);
      double val(1.0), a(shift[i]);
      for (int n = 0; n <= index[i]; ++n) {
        powers[i](index[i] - n) = val * cnk;
        cnk *= index[i] - n;
        cnk /= n + 1;
        val *= a;
      }
    }

    // product of (x-a)^k (y-b)^l (z-c)^m
    int idx[3];
    for (int i0 = 0; i0 <= index[0]; ++i0) {
      idx[0] = i0;
      double coef0 = powers[0](i0); 

      for (int i1 = 0; i1 <= index[1]; ++i1) {
        idx[1] = i1;
        double coef1 = powers[1](i1); 

        if (d == 2) {
          int pos = MonomialSetPosition(d_, idx);
          poly(i0 + i1, pos) = coef * coef0 * coef1;
        } else {
          for (int i2 = 0; i2 <= index[2]; ++i2) {
            idx[2] = i2;
            double coef2 = powers[2](i2); 

            int pos = MonomialSetPosition(d_, idx);
            poly(i0 + i1 + i2, pos) = coef * coef0 * coef1 * coef2;
          }
        }
      }
    }
  }

  poly.set_origin(origin);
  return poly;
}


/* ******************************************************************
* Calculate polynomial value at a given point. 
****************************************************************** */
double Polynomial::Value(const AmanziGeometry::Point& xp) const
{
  double sum(coefs_(0));
  AmanziGeometry::Point dx(xp - origin_);
  
  if (order_ > 0) {
    for (int i = 0; i < d_; ++i) {
      sum += dx[i] * coefs_(i + 1);
    }
  }

  for (auto it = begin(2); it < end(); ++it) {
    int n = it.PolynomialPosition();
    double tmp = coefs_(n);

    if (tmp != 0.0) {
      const int* index = it.multi_index();
      for (int i = 0; i < d_; ++i) {
        tmp *= std::pow(dx[i], index[i]);
      }
      sum += tmp;
    }
  }

  return sum;
}


/* ******************************************************************
* Change of coordinates: x = x0 + B * s 
* Note: the resulting polynomial is centered at the new local origin.
****************************************************************** */
void Polynomial::ChangeCoordinates(
    const AmanziGeometry::Point& x0, const std::vector<AmanziGeometry::Point>& B)
{
  int dnew = B.size();
  AMANZI_ASSERT(dnew > 0);

  // center polynomial at x0
  ChangeOrigin(x0);

  // populate new polynomial using different algorithms
  Polynomial tmp(dnew, order_);

  if (dnew == 1) {
    for (auto it = begin(); it < end(); ++it) {
      int m = it.MonomialSetOrder();
      int n = it.PolynomialPosition();

      double coef = coefs_(n);
      if (coef != 0.0) {
        const int* multi_index = it.multi_index();
        for (int i = 0; i < d_; ++i) {
          coef *= std::pow(B[0][i], multi_index[i]);  
        }
        tmp(m, 0) += coef;
      }
    }
  }

  else if (dnew == 2) {
    Polynomial poly0(dnew, 1), poly1(dnew, 1), poly2(dnew, 1);

    poly0(1) = B[0][0];
    poly0(2) = B[1][0];

    poly1(1) = B[0][1];
    poly1(2) = B[1][1];

    poly2(1) = B[0][2];
    poly2(2) = B[1][2];

    // lazy computation of powers of three binomials
    Polynomial one(dnew, 0);
    one(0) = 1.0;

    std::vector<Polynomial> power0, power1, power2;
    power0.push_back(one);
    power1.push_back(one);
    power2.push_back(one);

    for (auto it = begin(); it < end(); ++it) {
      double coef = coefs_(it.PolynomialPosition());

      if (coef != 0.0) {
        const int* idx = it.multi_index();

        for (int k = power0.size() - 1; k < idx[0]; ++k)
           power0.push_back(power0[k] * poly0);

        for (int k = power1.size() - 1; k < idx[1]; ++k)
           power1.push_back(power1[k] * poly1);

        for (int k = power2.size() - 1; k < idx[2]; ++k)
           power2.push_back(power2[k] * poly2);

        tmp += (coef * power0[idx[0]]) * power1[idx[1]] * power2[idx[2]];        
      }
    }
  }  
  
  *this = tmp;
}


/* ******************************************************************
* Inverse change of coordinates: s = B^+ (x - x0) for polynomial
* centered at the origin.
* Note: resulting polynomial is centered at x0.
****************************************************************** */
void Polynomial::InverseChangeCoordinates(
    const AmanziGeometry::Point& x0, const std::vector<AmanziGeometry::Point>& B)
{
  int dnew = x0.dim();

  // new polynomial will be centered at x0
  Polynomial tmp(dnew, order_);
  tmp.set_origin(x0);

  // populate new polynomial using different algorithms
  if (d_ == 1 && dnew == 2) {
    int i = (fabs(B[0][0]) > fabs(B[0][1])) ? 0 : 1;
    double scale = 1.0 / B[0][i];

    for (int m = 0; m < size_; ++m) {
      tmp(m, i * m) = coefs_(m) * std::pow(scale, m);
    }
  } else if (d_ == 1 && dnew == 3) {
    int i = (fabs(B[0][0]) > fabs(B[0][1])) ? 0 : 1;
    if (fabs(B[0][2]) > fabs(B[0][i])) i = 2;
    double scale = 1.0 / B[0][i];

    int multi_index[3] = {0, 0, 0};

    for (int n = 0; n < size_; ++n) {
      multi_index[i] = n;
      int m = PolynomialPosition(3, multi_index);
      tmp(m) = coefs_(n) * std::pow(scale, n);
    }
  } else if (d_ == 2) {
    // find monor with the largest determinant
    int i0(0), i1(1), i2;
    double det01, det02, det12, tmp01, tmp02, tmp12;
    
    det01 = B[0][0] * B[1][1] - B[0][1] * B[1][0];
    det02 = B[0][0] * B[1][2] - B[0][2] * B[1][0];
    det12 = B[0][1] * B[1][2] - B[0][2] * B[1][1];

    tmp01 = std::fabs(det01); 
    tmp02 = std::fabs(det02);
    tmp12 = std::fabs(det12);

    if (tmp02 >= std::max(tmp01, tmp12)) {
      i1 = 2;
      det01 = det02;
    } else if (tmp12 >= std::max(tmp01, tmp02)) {
      i0 = 1;
      i1 = 2;
      det01 = det12;
    }

    // invert indirectly the minor
    Polynomial poly0(dnew, 1), poly1(dnew, 1);
    poly0(1 + i0) = B[1][i1] / det01;
    poly0(1 + i1) =-B[1][i0] / det01;
    poly0.set_origin(x0);

    poly1(1 + i0) =-B[0][i1] / det01;
    poly1(1 + i1) = B[0][i0] / det01;
    poly1.set_origin(x0);

    // lazy computation of powers of two binomials
    Polynomial one(dnew, 0);
    one.set_origin(x0);
    one(0) = 1.0;

    std::vector<Polynomial> power0, power1;
    power0.push_back(one);
    power1.push_back(one);

    for (auto it = begin(); it < end(); ++it) {
      const int* idx = it.multi_index();
      double coef = coefs_(it.PolynomialPosition());

      if (coef != 0.0) {
        for (int k = power0.size() - 1; k < idx[0]; ++k)
           power0.push_back(power0[k] * poly0);

        for (int k = power1.size() - 1; k < idx[1]; ++k)
           power1.push_back(power1[k] * poly1);

        tmp += coef * (power0[idx[0]] * power1[idx[1]]);        
      }
    }
  } 
  
  *this = tmp;
}


/* ******************************************************************
* Laplacian operator
****************************************************************** */
Polynomial Polynomial::Laplacian()
{
  int order = std::max(0, order_ - 2);

  Polynomial tmp(d_, order);

  int index[3];
  for (auto it = begin(); it < end(); ++it) {
    int k = it.MonomialSetOrder();
    if (k > 1) {
      const int* idx = it.multi_index();
      int n = it.PolynomialPosition();
      double val = coefs_(n);

      for (int i = 0; i < d_; ++i) {
        for (int j = 0; j < d_; ++j) index[j] = idx[j];

        if (index[i] > 1) {
          index[i] -= 2;
          int m = MonomialSetPosition(d_, index);
          tmp(k - 2, m) += val * idx[i] * (idx[i] - 1);
        }
      }
    }
  }

  return tmp;
}


/* ******************************************************************
* Fancy I/O
****************************************************************** */
std::ostream& operator << (std::ostream& os, const Polynomial& p)
{
  int d = p.dimension();
  os << "polynomial: order=" << p.order() << " d=" << d
     << " size=" << p.size() << std::endl;
  for (auto it = p.begin(); it < p.end(); ++it) {
    int k = it.MonomialSetOrder();
    int m = it.MonomialSetPosition();
    double val = p(k, m);

    if (m == 0) os << "k=" << k << "  coefs:";
    if (m > 0 && val >= 0) os << " +";
    os << " " << val << " ";

    const int* index = it.multi_index();
    if (index[0] == 1) os << "x";
    if (index[0] > 1)  os << "x^" << index[0];
    if (index[1] == 1) os << "y";
    if (index[1] > 1)  os << "y^" << index[1];
    if (index[2] == 1) os << "z";
    if (index[2] > 1)  os << "z^" << index[2];

    if (m == MonomialSpaceDimension(d, k) - 1) os << std::endl;
  } 
  os << "origin: " << p.origin() << std::endl;
  return os;
}


/* ******************************************************************
* Friends: copy optimized product of polynomial
****************************************************************** */
Polynomial operator*(const Polynomial& poly1, const Polynomial& poly2)
{
  int d = poly1.dimension();
  AMANZI_ASSERT(d == poly2.dimension());  // FIXME
  AMANZI_ASSERT(poly1.origin() == poly2.origin());

  int order1 = poly1.order();
  int order2 = poly2.order();
  Polynomial tmp(d, order1 + order2);
  tmp.set_origin(poly1.origin());

  if (order1 == 0) {
    for (int n = 0; n < tmp.size(); ++n) {
      tmp(n) = poly1(0) * poly2(n);
    }
    return tmp; 
  }

  if (order2 == 0) {
    for (int n = 0; n < tmp.size(); ++n) {
      tmp(n) = poly2(0) * poly1(n);
    }
    return tmp; 
  }

  int index[3];
  for (auto it1 = poly1.begin(); it1 < poly1.end(); ++it1) {
    const int* idx1 = it1.multi_index();
    int n1 = it1.PolynomialPosition();
    double val1 = poly1(n1);
    if (val1 == 0.0) continue;

    for (auto it2 = poly2.begin(); it2 < poly2.end(); ++it2) {
      const int* idx2 = it2.multi_index();
      int n2 = it2.PolynomialPosition();
      double val2 = poly2(n2);

      for (int i = 0; i < d; ++i) {
        index[i] = idx1[i] + idx2[i];
      }
      int l = PolynomialPosition(d, index);
      tmp(l) += val1 * val2;
    }
  }

  return tmp;
}

}  // namespace WhetStone
}  // namespace Amanzi


