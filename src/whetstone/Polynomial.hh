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

#ifndef AMANZI_WHETSTONE_POLYNOMIAL_HH_
#define AMANZI_WHETSTONE_POLYNOMIAL_HH_

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <vector>

namespace Amanzi {
namespace WhetStone {

class Iterator {
 public:
  Iterator() {};
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
    if (multi_index_[0] == 0) {
      k_++;
      m_ = 0;
      multi_index_[0] = k_;
      multi_index_[1] = 0;
      multi_index_[2] = 0;
    } else {
      m_++;
      multi_index_[0]--;
      multi_index_[1] = k_ - multi_index_[0];
    }
    return *this;
  }

  int end() { return k_; }

  // access
  const int* multi_index() const { return multi_index_; }

 private:
  int k_;  // current monomials order
  int m_;  // current position in the list of monomials
  int multi_index_[3];
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
      os << std::setw(12) << std::setprecision(12) << v[i] << " ";
    }
    return os;
  }

 private:
  int d_, order_;
  std::vector<double> coefs_;
};
 

class Polynomial {
 public:
  Polynomial() : d_(0), order_(-1), size_(0) {};
  Polynomial(int d, int order);

  // elemental operations with polynomials
  int MonomialPosition(const int* multi_index) const;
  int PolynomialPosition(const int* multi_index) const;

  // iterators
  Iterator& begin() { return it_.begin(); }
  int end() { return order_; }

  // access
  int order() const { return order_; }
  int size() const { return size_; }

  Monomial& monomials(int i) { return coefs_[i]; }
  const Monomial& monomials(int i) const { return coefs_[i]; }

  // output 
  friend std::ostream& operator << (std::ostream& os, const Polynomial& p) {
    for (int i = 0; i <= p.order(); i++) {
      const Monomial& m = p.monomials(i);
      os << m << "\n";
    } 
    return os;
  }

 protected:
  Iterator it_; 

 private:
  int d_, order_, size_;
  std::vector<Monomial> coefs_;
};

}  // namespace WhetStone
}  // namespace Amanzi

#endif

