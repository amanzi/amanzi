/*
  WhetStone, version 2.1
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Operations with polynomials in 2D and 3D.
*/

#ifndef AMANZI_WHETSTONE_POLYNOMIAL_HH_
#define AMANZI_WHETSTONE_POLYNOMIAL_HH_

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <vector>

#include "Point.hh"

namespace Amanzi {
namespace WhetStone {

class Iterator {
 public:
  Iterator() : d_(0) {};
  Iterator(int d) : d_(d) {};
  ~Iterator() {};

  Iterator& begin() {
    k_ = 0;
    m_ = 0;
    multi_index_[0] = 0;
    multi_index_[1] = 0;
    multi_index_[2] = 0;

    return *this;
  }

  Iterator& operator++() {  // prefix only
    if (d_ == 2) {
      if (multi_index_[0] == 0) {
        k_++;  // next group of monomials
        m_ = 0;
        multi_index_[0] = k_;
        multi_index_[1] = 0;
        multi_index_[2] = 0;
      } else {
        m_++;
        multi_index_[0]--;
        multi_index_[1] = k_ - multi_index_[0];
      }
    } else if (d_ == 3) {
      if (multi_index_[0] == 0 && multi_index_[1] == 0) {
        k_++;  // next group of monomials
        m_ = 0;
        multi_index_[0] = k_;
        multi_index_[1] = 0;
        multi_index_[2] = 0;
      } else if (multi_index_[1] == 0) {
        m_++;
        multi_index_[0]--;
        multi_index_[1] = k_ - multi_index_[0];
        multi_index_[2] = 0;
      } else {
        m_++;
        multi_index_[1]--;
        multi_index_[2] = k_ - multi_index_[0] - multi_index_[1];
      }
    }
    return *this;
  }

  int end() { return k_; }

  // access
  int MonomialOrder() const { return k_; }
  int MonomialPosition() const { return m_; }
  const int* multi_index() const { return multi_index_; }

 private:
  int k_;  // current monomials order
  int m_;  // current position in the list of monomials
  int multi_index_[3];
  int d_;
};


class Monomial {
 public:
  Monomial() : d_(0), order_(-1) {};
  Monomial(int d, int order) : d_(d), order_(order) {
    int nk = (order_ == 0) ? 1 : d_;
    for (int i = 1; i < order_; ++i) {
      nk *= d_ + i;
      nk /= i + 1;
    }
    coefs_.resize(nk, 0.0);
  }
  ~Monomial() {};

  // access
  int order() const { return order_; }
  int size() const { return coefs_.size(); }

  std::vector<double>& coefs() { return coefs_; } 
  const std::vector<double>& coefs() const { return coefs_; } 

  // output 
  friend std::ostream& operator << (std::ostream& os, const Monomial& m) {
    const auto& v = m.coefs();
    for (int i = 0; i < v.size(); i++) {
      os << std::setw(10) << std::setprecision(10) << v[i] << " ";
    }
    return os;
  }

 private:
  // direct memory access operations
  double& operator()(int i) { return coefs_[i]; }
  const double& operator()(int i) const { return coefs_[i]; }

  friend class Polynomial;

 private:
  int d_, order_;
  std::vector<double> coefs_;
};
 

class Polynomial {
 public:
  Polynomial() : d_(0), order_(-1), size_(0) {};
  Polynomial(int d, int order);

  // reshape polynomial with erase (optionally) memory
  void Reshape(int d, int order, bool reset = false);

  // resets all coefficients to zero
  void Reset();

  // typical operations with polynomials
  // -- polynomial values
  double Value(const AmanziGeometry::Point& xp);

  // -- operators (ring algebra)
  Polynomial& operator+=(const Polynomial& poly);
  Polynomial& operator-=(const Polynomial& poly);
  Polynomial& operator*=(const Polynomial& poly);

  // -- multi index defines both monomial order and current monomial
  int MonomialPosition(const int* multi_index) const;
  int PolynomialPosition(const int* multi_index) const;

  // iterators
  Iterator& begin() const { return it_.begin(); }
  int end() const { return order_; }

  // access
  int dimension() const { return d_; }
  int order() const { return order_; }
  int size() const { return size_; }

  Monomial& monomials(int i) { return coefs_[i]; }
  const Monomial& monomials(int i) const { return coefs_[i]; }

  // output 
  friend std::ostream& operator << (std::ostream& os, const Polynomial& p) {
    for (int i = 0; i <= p.order(); i++) {
      const Monomial& m = p.monomials(i);
      os << "k=" << i << "  coefs: " << m << "\n";
    } 
    return os;
  }

 protected:
  mutable Iterator it_; 

 private:
  int d_, order_, size_;
  std::vector<Monomial> coefs_;
};

}  // namespace WhetStone
}  // namespace Amanzi

#endif

