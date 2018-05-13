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
    coefs_[i] = Monomial(d_, i);
    size_ += coefs_[i].size();
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
    coefs_[i] = Monomial(d_, i);
    size_ += coefs_[i].size();
  }

  int l = MonomialPosition(multi_index);
  coefs_[order_](l) = factor;
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
      coefs_[i] = Monomial(d_, i);
    }
  } else if (order_ != order) {
    coefs_.resize(order + 1);
    for (int i = order_ + 1; i <= order; ++i) {
      coefs_[i] = Monomial(d_, i);
    }
    order_ = order;
  }

  size_ = 0;
  for (int i = 0; i <= order_; ++i) {
    size_ += coefs_[i].size();
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
    auto& mono = coefs_[i].coefs();
    for (int k = 0; k < mono.size(); ++k) {
      mono[k] += poly.monomials(i).coefs()[k];
    }
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
    auto& mono = coefs_[i].coefs();
    for (int k = 0; k < mono.size(); ++k) {
      mono[k] -= poly.monomials(i).coefs()[k];
    }
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
    int k1 = it1.MonomialOrder();
    int m1 = it1.MonomialPosition();
    double val1 = arg1.monomials(k1).coefs()[m1];
    if (val1 == 0.0) continue;

    for (auto it2 = arg2->begin(); it2.end() <= arg2->end(); ++it2) {
      const int* idx2 = it2.multi_index();
      int k2 = it2.MonomialOrder();
      int m2 = it2.MonomialPosition();
      double val2 = arg2->monomials(k2).coefs()[m2];

      int n(0);
      for (int i = 0; i < d_; ++i) {
        index[i] = idx1[i] + idx2[i];
        n += index[i];
      }
      int l = MonomialPosition(index);
      coefs_[n](l) += val1 * val2;
    }
  }

  return *this;
}


Polynomial& Polynomial::operator*=(double val)
{
  for (int i = 0; i <= order_; ++i) {
    std::vector<double>& tmp = coefs_[i].coefs();
    for (auto it = tmp.begin(); it != tmp.end(); ++it) *it *= val;
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
    std::vector<std::vector<Polynomial> > powers(d_);

    for (int i = 0; i < d_; ++i) {
      powers[i].resize(order_ + 1);

      for (int k = 0; k <= order_; ++k) {
        int index[3] = {0, 0, 0};
        index[i] = k;

        Polynomial& p = powers[i][k];
        p.Reshape(d_, k, true);

        int cnk(1);
        double val(1.0), a(shift[i]);
        for (int n = 0; n <= k; ++n) {
          int l = p.MonomialPosition(index);
          index[i]--;

          p.monomials(k - n)(l) = val * cnk;
          cnk *= (k - n);
          cnk /= (n + 1);
          val *= a;
        }
      }
    }

    // iterate over polynomial and sum up products
    Polynomial rebased(d_, order_);
    for (auto it = begin(); it.end() <= end(); ++it) {
      int k = it.MonomialOrder();
      int m = it.MonomialPosition();
      double coef = monomials(k).coefs()[m];
      if (coef == 0.0) continue;

      const int* index = it.multi_index();
      Polynomial tmp(powers[0][index[0]]);
      for (int i = 1; i < d_; ++i) tmp *= powers[i][index[i]];

      tmp *= coef; 
      rebased += tmp;
    }
    
    *this = rebased;
  }
  origin_ = origin;   
}


/* ******************************************************************
* Reset all coefficients to the same number
****************************************************************** */
void Polynomial::PutScalar(double val)
{
  for (int i = 0; i <= order_; ++i) {
    std::vector<double>& tmp = coefs_[i].coefs();
    for (auto it = tmp.begin(); it != tmp.end(); ++it) *it = val;
  }
}


/* ******************************************************************
* Set polynomial coefficients to entries of the given vector.
****************************************************************** */
void Polynomial::SetPolynomialCoefficients(const DenseVector& coefs)
{
  ASSERT(size_ == coefs.NumRows());

  const double* data = coefs.Values();
  for (int i = 0; i <= order_; ++i) {
    std::vector<double>& tmp = coefs_[i].coefs();
    for (auto it = tmp.begin(); it != tmp.end(); ++it) {
      *it = *data; 
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
  for (int i = 0; i <= order_; ++i) {
    const std::vector<double>& tmp = coefs_[i].coefs();
    for (auto it = tmp.begin(); it != tmp.end(); ++it) {
      *data = *it; 
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
    int k = it.MonomialOrder();
    int l = it.MonomialPosition();
    const int* index = it.multi_index();

    double tmp(coefs_[k](l));
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

  for (int n = 0; n <= order_; ++n) {
    const std::vector<double>& coefs = coefs_[n].coefs();
    for (int i = 0; i < coefs.size(); ++i) {
      tmp = std::max(tmp, fabs(coefs[i]));
    }
  }

  return tmp;
}


/* ******************************************************************
* Position in monomial defined by multi_index. 2D algorithm
****************************************************************** */
int Polynomial::MonomialPosition(const int* multi_index) const
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

  return nk + MonomialPosition(multi_index);
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
    int m = it.MonomialOrder();
    int k = it.MonomialPosition();
    if (dnew == 1) {
      double coef(coefs_[m](k));
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
      int m = it.MonomialOrder();
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
    int k = it.MonomialOrder();
    if (k > 1) {
      const int* idx = it.multi_index();
      int m = it.MonomialPosition();
      double val = coefs_[k](m);

      for (int i = 0; i < d_; ++i) {
        for (int j = 0; j < d_; ++j) index[j] = idx[j];

        if (index[i] > 1) {
          index[i] -= 2;
          m = MonomialPosition(index);
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
    int k = it.MonomialOrder();
    int m = it.MonomialPosition();
    double val = p.monomials(k).coefs()[m];

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

    if (m == p.monomials(k).size() - 1) os << std::endl;
  } 
  os << "origin: " << p.origin() << std::endl;
  return os;
}

}  // namespace WhetStone
}  // namespace Amanzi


