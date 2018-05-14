/*
  WhetStone, version 2.1
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Operations with polynomials.
*/

#include "Monomial.hh"
#include "Polynomial.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Constructor of zero polynomial.
****************************************************************** */
Polynomial::Polynomial(int d, int order) : 
    d_(d), order_(order), origin_(d)
{
  size_ = 0;
  coefs_.resize(order_ + 1);
  for (int i = 0; i <= order_; ++i) {
    int msize = MonomialSpaceDimension(d_, i);
    coefs_[i].Reshape(msize);
    coefs_[i].PutScalar(0.0);
    size_ += msize;
  }
}


/* ******************************************************************
* Constructor of a polynomial with a single term:
*    p(x) = factor * (x)^multi_index
****************************************************************** */
Polynomial::Polynomial(int d, const int* multi_index, double factor) : 
    d_(d), origin_(d)
{
  order_ = 0;
  for (int i = 0; i < d_; ++i) order_ += multi_index[i];

  size_ = 0;
  coefs_.resize(order_ + 1);
  for (int i = 0; i <= order_; ++i) {
    int msize = MonomialSpaceDimension(d_, i);
    coefs_[i].Reshape(msize);
    coefs_[i].PutScalar(0.0);
    size_ += msize;
  }

  int l = MonomialSetPosition(multi_index);
  coefs_[order_](l) = factor;
}


/* ******************************************************************
* Constructor of a polynomial from a monomial.
****************************************************************** */
Polynomial::Polynomial(const Monomial& mono)
{
  d_ = mono.dimension();
  order_ = mono.order();
  origin_ = mono.origin();

  size_ = 0;
  coefs_.resize(order_ + 1);
  for (int i = 0; i <= order_; ++i) {
    int msize = MonomialSpaceDimension(d_, i);
    coefs_[i].Reshape(msize);
    coefs_[i].PutScalar(0.0);
    size_ += msize;
  }

  int l = MonomialSetPosition(mono.multi_index());
  coefs_[order_](l) = mono.coef();
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

    coefs_.clear();
    coefs_.resize(order_ + 1);
    for (int i = 0; i <= order_; ++i) {
      int msize = MonomialSpaceDimension(d_, i);
      coefs_[i].Reshape(msize);
      coefs_[i].PutScalar(0.0);
    }
  } else if (order_ != order) {
    coefs_.resize(order + 1);
    for (int i = order_ + 1; i <= order; ++i) {
      int msize = MonomialSpaceDimension(d_, i);
      coefs_[i].Reshape(msize);
      coefs_[i].PutScalar(0.0);
    }
    order_ = order;
  }

  size_ = 0;
  for (int i = 0; i <= order_; ++i) {
    size_ += coefs_[i].NumRows();
  }

  if (reset) PutScalar(0.0);
}


/* ******************************************************************
* Implemented ring algebra operations.
* NOTE: implementation is order independent.
****************************************************************** */
Polynomial& Polynomial::operator+=(const Polynomial& poly)
{
  ASSERT(d_ == poly.dimension());  // FIXME
  ASSERT(origin_ == poly.origin());

  int order_max = std::max(order_, poly.order());
  Reshape(d_, order_max);

  for (int i = 0; i <= poly.order(); ++i) {
    coefs_[i] += poly.MonomialSet(i);
  }

  return *this;
}


Polynomial& Polynomial::operator-=(const Polynomial& poly)
{
  ASSERT(d_ == poly.dimension());  // FIXME
  ASSERT(origin_ == poly.origin());

  int order_max = std::max(order_, poly.order());
  Reshape(d_, order_max);

  for (int i = 0; i <= poly.order(); ++i) {
    coefs_[i] -= poly.MonomialSet(i);
  }

  return *this;
}


Polynomial& Polynomial::operator*=(const Polynomial& poly)
{
  ASSERT(d_ == poly.dimension());  // FIXME
  ASSERT(origin_ == poly.origin());

  Polynomial arg1(*this);
  const Polynomial* arg2 = &poly;
  if (this == arg2) arg2 = &arg1; 

  int order_prod = order_ + poly.order();
  Reshape(d_, order_prod, true);

  int index[3];
  for (auto it1 = arg1.begin(); it1.end() <= arg1.end(); ++it1) {
    const int* idx1 = it1.multi_index();
    int k1 = it1.MonomialSetOrder();
    int m1 = it1.MonomialSetPosition();
    double val1 = arg1(k1, m1);
    if (val1 == 0.0) continue;

    for (auto it2 = arg2->begin(); it2.end() <= arg2->end(); ++it2) {
      const int* idx2 = it2.multi_index();
      int k2 = it2.MonomialSetOrder();
      int m2 = it2.MonomialSetPosition();
      double val2 = arg2->operator()(k2, m2);

      int n(0);
      for (int i = 0; i < d_; ++i) {
        index[i] = idx1[i] + idx2[i];
        n += index[i];
      }
      int l = MonomialSetPosition(index);
      coefs_[n](l) += val1 * val2;
    }
  }

  return *this;
}


Polynomial& Polynomial::operator*=(double val)
{
  for (int i = 0; i <= order_; ++i) {
    coefs_[i] *= val;
  }

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
      coefs_[0](0) += coefs_[1](i) * shift[i];
    }
  } else if (order_ > 1) {
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
    for (auto it = begin(); it.end() <= end(); ++it) {
      int k = it.MonomialSetOrder();
      int m = it.MonomialSetPosition();
      double coef = coefs_[k](m);
      if (coef == 0.0) continue;

      const int* index = it.multi_index();
      // product of (x-a)^k (y-b)^l (z-c)^m
      int idx[3];
      Polynomial tmp(d_, k);
      for (int i0 = 0; i0 <= index[0]; ++i0) {
        idx[0] = i0;
        double coef0 = powers[0][index[0]](i0); 

        for (int i1 = 0; i1 <= index[1]; ++i1) {
          idx[1] = i1;
          double coef1 = powers[1][index[1]](i1); 

          if (d_ == 2) {
            int pos = MonomialSetPosition(idx);
            tmp(i0 + i1, pos) = coef * coef0 * coef1;
          } else {
            for (int i2 = 0; i2 <= index[2]; ++i2) {
              idx[2] = i2;
              double coef2 = powers[2][index[2]](i2); 

              int pos = MonomialSetPosition(idx);
              tmp(i0 + i1 + i2, pos) = coef * coef0 * coef1 * coef2;
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
  double coef = mono.coef();

  Polynomial poly(d, order);

  if (order == 0) {
    poly(0, 0) = coef;
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
          int pos = MonomialSetPosition(idx);
          poly(i0 + i1, pos) = coef * coef0 * coef1;
        } else {
          for (int i2 = 0; i2 <= index[2]; ++i2) {
            idx[2] = i2;
            double coef2 = powers[2](i2); 

            int pos = MonomialSetPosition(idx);
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
* Reset all coefficients to the same number
****************************************************************** */
void Polynomial::PutScalar(double val)
{
  for (int i = 0; i <= order_; ++i) {
    coefs_[i] = val;
  }
}


/* ******************************************************************
* Set polynomial coefficients to entries of the given vector.
****************************************************************** */
void Polynomial::SetPolynomialCoefficients(const DenseVector& coefs)
{
  ASSERT(size_ == coefs.NumRows());

  const double* data = coefs.Values();
  for (int k = 0; k <= order_; ++k) {
    int mk = coefs_[k].NumRows();
    double* it = coefs_[k].Values();

    for (int i = 0; i < mk; ++i) {
      *it = *data; 
      it++;
      data++;
    }
  }
}


/* ******************************************************************
* Copy polynomial coefficients to a vector. Vector is resized.
****************************************************************** */
void Polynomial::GetPolynomialCoefficients(DenseVector& coefs) const
{
  coefs.Reshape(size_);

  double* data = coefs.Values();
  for (int k = 0; k <= order_; ++k) {
    int mk = coefs_[k].NumRows();
    const double* it = coefs_[k].Values();

    for (int i = 0; i < mk; ++i) {
      *data = *it; 
      it++;
      data++;
    }
  }
}


/* ******************************************************************
* Calculate polynomial value
****************************************************************** */
double Polynomial::Value(const AmanziGeometry::Point& xp) const
{
  double sum(0.0);

  for (auto it = begin(); it.end() <= end(); ++it) {
    int k = it.MonomialSetOrder();
    int l = it.MonomialSetPosition();
    const int* index = it.multi_index();

    double tmp = coefs_[k](l);
    if (tmp != 0.0) {
      for (int i = 0; i < d_; ++i) {
        tmp *= std::pow(xp[i] - origin_[i], index[i]);
      }
      sum += tmp;
    }
  }

  return sum;
}


/* ******************************************************************
* Calculate polynomial maximum norm
****************************************************************** */
double Polynomial::NormMax() const
{
  double tmp(0.0);

  for (int k = 0; k <= order_; ++k) {
    tmp = std::max(tmp, coefs_[k].NormMax());
  }

  return tmp;
}


/* ******************************************************************
* Position in monomial defined by multi_index. 2D algorithm
****************************************************************** */
int Polynomial::MonomialSetPosition(const int* multi_index) const
{
  int m(multi_index[1]);
  if (d_ == 3) {
    int n = multi_index[1] + multi_index[2];
    m = n * (n + 1) / 2 + multi_index[2];
  }
  return m;
}


/* ******************************************************************
* Position in polynomial defined by multi_index. 2D algorithm.
****************************************************************** */
int Polynomial::PolynomialPosition(const int* multi_index) const
{
  // current monomial order
  int k = 0;
  for (int i = 0; i < d_; ++i) k += multi_index[i];

  // sum of previous monomials
  int nk = (d_ == 2) ? k * (k + 1) / 2 : k * (k + 1) * (k + 2) / 6; 

  return nk + MonomialSetPosition(multi_index);
}


/* ******************************************************************
* Change of coordinates: x = x0 + B * s 
* Note: resulting polynomial is centered at new origin.
****************************************************************** */
void Polynomial::ChangeCoordinates(
    const AmanziGeometry::Point& x0, const std::vector<AmanziGeometry::Point>& B)
{
  int dnew = B.size();
  ASSERT(dnew > 0);

  // center polynomial at x0
  ChangeOrigin(x0);

  // populate new polynomial using different algorithms
  Polynomial tmp(dnew, order_);
  for (auto it = begin(); it.end() <= end(); ++it) {
    const int* multi_index = it.multi_index();
    int m = it.MonomialSetOrder();
    int k = it.MonomialSetPosition();
    if (dnew == 1) {
      double coef = coefs_[m](k);
      for (int i = 0; i < d_; ++i) {
        coef *= std::pow(B[0][i], multi_index[i]);  
      }
      tmp(m, 0) += coef;
    }
  }  
  
  *this = tmp;
}


/* ******************************************************************
* Inverse change of coordinates: s = B^+ (x - x0) 
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
  if (d_ == 1) {
    int i = (fabs(B[0][0]) > fabs(B[0][1])) ? 0 : 1;
    double scale = 1.0 / B[0][i];

    for (auto it = begin(); it.end() <= end(); ++it) {
      int m = it.MonomialSetOrder();
      tmp(m, i * m) = coefs_[m](0) * std::pow(scale, m);
    }
  }  
  
  *this = tmp;
}


/* ******************************************************************
* Special operations: Laplacian
****************************************************************** */
Polynomial Polynomial::Laplacian()
{
  int order = std::max(0, order_ - 2);

  Polynomial tmp(d_, order);

  int index[3];
  for (auto it = begin(); it.end() <= end(); ++it) {
    int k = it.MonomialSetOrder();
    if (k > 1) {
      const int* idx = it.multi_index();
      int m = it.MonomialSetPosition();
      double val = coefs_[k](m);

      for (int i = 0; i < d_; ++i) {
        for (int j = 0; j < d_; ++j) index[j] = idx[j];

        if (index[i] > 1) {
          index[i] -= 2;
          m = MonomialSetPosition(index);
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
  os << "polynomial: order=" << p.order() << " d=" << p.dimension() 
     << " size=" << p.size() << std::endl;
  for (auto it = p.begin(); it.end() <= p.end(); ++it) {
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

    if (m == p.MonomialSet(k).NumRows() - 1) os << std::endl;
  } 
  os << "origin: " << p.origin() << std::endl;
  return os;
}

}  // namespace WhetStone
}  // namespace Amanzi


