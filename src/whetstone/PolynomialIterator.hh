/*
  WhetStone, version 2.1
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Iterator starts with a given monomial order (e.g. x^2) runs
  through all monomials of the same order (resp., xy, y^2), and
  jumps to the next order (resp., x^3). 
*/

#ifndef AMANZI_WHETSTONE_POLYNOMIAL_ITERATOR_HH_
#define AMANZI_WHETSTONE_POLYNOMIAL_ITERATOR_HH_

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <vector>

#include "Point.hh"

#include "DenseVector.hh"

namespace Amanzi {
namespace WhetStone {

class PolynomialIterator {
 public:
  PolynomialIterator() : d_(0) {};
  PolynomialIterator(int d) : d_(d) {};
  ~PolynomialIterator() {};

  // set iterator to the monomials group of order k0
  PolynomialIterator& begin(int k0 = 0) {
    k_ = k0;
    m_ = 0;
    multi_index_[0] = k0;
    multi_index_[1] = 0;
    multi_index_[2] = 0;
    count_ = 0;

    return *this;
  }

  // Move iterator either to the next monomial in the group or
  // to the begining of the next group of monomials.
  // Only prefix version of this operator is used in the code.
  PolynomialIterator& operator++() {
    if (d_ == 1) {
      k_++;
      m_ = 0;
      multi_index_[0] = k_;
    } else if (d_ == 2) {
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
    count_++;

    return *this;
  }

  // One way to terminate a for-loop is to capture the moment when
  // the iterator moved to the next group of monomials. Returning
  // the current monomial order is not ideal solution, but robust.
  int end() { return k_; }

  // access
  int MonomialOrder() const { return k_; }
  int MonomialPosition() const { return m_; }
  int PolynomialPosition() const { return count_; }
  const int* multi_index() const { return multi_index_; }

 private:
  int k_;  // current monomials order
  int m_;  // current position in the list of monomials
  int multi_index_[3];
  int d_;
  int count_;  // iterator count
};

}  // namespace WhetStone
}  // namespace Amanzi

#endif

