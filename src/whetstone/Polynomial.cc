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
* Constructor
****************************************************************** */
Polynomial::Polynomial(int d, int order) : d_(d), order_(order), it_(d)
{
  size_ = 0;
  coefs_.resize(order_ + 1);
  for (int i = 0; i <= order_; ++i) {
    coefs_[i] = Monomial(d_, i);
    size_ += coefs_[i].size();
  }
}


/* ******************************************************************
* Resets all coefficients to zero
****************************************************************** */
void Polynomial::Reset()
{
  for (int i = 0; i <= order_; ++i) {
    std::vector<double>& tmp = coefs_[i].coefs();
    for (auto it = tmp.begin(); it != tmp.end(); ++it) *it = 0.0;
  }
}


/* ******************************************************************
* Calculate polynomial value
****************************************************************** */
double Polynomial::Value(const AmanziGeometry::Point& xp)
{
  double sum(0.0);

  for (auto it = begin(); it.end() <= end(); ++it) {
    const int* index = it.multi_index();
    int l = MonomialPosition(index);

    int k(0);
    double tmp(1.0);
    for (int i = 0; i < d_; ++i) {
      tmp *= std::pow(xp[i], index[i]);
      k += index[i];
    }
    sum += tmp * coefs_[k](l);
  }

  return sum;
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

}  // namespace WhetStone
}  // namespace Amanzi


