/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Simple iterator on top of the polynomial iterator design to simplify
  implementation of discretization scheme with vector degrees of freedom
  Iterator runs faster through polynomial dimension and then jumps to the
  next vector component. The polynomial iterator stops at the maximum
  polynomial degree which is the same for all components in the vector.
*/

#ifndef AMANZI_WHETSTONE_VECTOR_OBJECTS_ITERATOR_HH_
#define AMANZI_WHETSTONE_VECTOR_OBJECTS_ITERATOR_HH_

#include <cstdlib>

#include "PolynomialIterator.hh"
#include "VectorObjects.hh"

namespace Amanzi {
namespace WhetStone {

template<class T>
class VectorObjectsIterator {
 public:
  VectorObjectsIterator();
  VectorObjectsIterator(int d, int size, int order);
  ~VectorObjectsIterator() {};

  VectorObjectsIterator begin();
  VectorObjectsIterator end();
};


template<>
class VectorObjectsIterator<Polynomial> {
 public:
  VectorObjectsIterator<Polynomial>() 
    : size_(0), order_(0), it_(0), it1_(0) {};
  VectorObjectsIterator<Polynomial>(int d, int size, int order)
    : size_(size), order_(order), it_(d), it1_(d) {
    it1_.begin(order_ + 1);
  }
  ~VectorObjectsIterator<Polynomial>() {};

  // sometimes we need to iterate without creating an actual object
  // -- set iterator to the initial position
  VectorObjectsIterator<Polynomial> begin() {
    m_ = 0;
    count_ = 0;
    it_.begin(0);
    return *this;
  }

  // -- set iterator to the final position
  VectorObjectsIterator<Polynomial> end() {
    m_ = size_;
    it_ = it1_;
    count_ = size_ * it1_.PolynomialPosition();
    return *this;
  }

  // Move iterator either to the next term in the polynomial or 
  // to the vector component.
  VectorObjectsIterator<Polynomial>& operator++() {
    ++it_;
    if (!(it_ < it1_)) {
      m_++;
      if (m_ < size_) it_.begin(0);
    }
    count_++;

    return *this;
  }

  // Comparison of iterators
  friend bool operator<(const VectorObjectsIterator<Polynomial>& it1, const VectorObjectsIterator<Polynomial>& it2) {
    return it1.VectorPolynomialPosition() < it2.VectorPolynomialPosition();
  }

  // access
  int MonomialSetOrder() const { return it_.MonomialSetOrder(); }
  int MonomialSetPosition() const { return it_.MonomialSetPosition(); }
  int PolynomialPosition() const { return it_.PolynomialPosition(); }
  const int* multi_index() const { return it_.multi_index(); }

  int VectorPolynomialPosition() const { return count_; }
  int VectorComponent() const { return m_; }

 private:
  int m_;  // current component in the vector
  int size_;  // vector size
  int order_;  // polynomial degree
  int count_;  // global iterator counter

  PolynomialIterator it_;  // polynomial iterator
  PolynomialIterator it1_; // pointer after the end of polynomial
};

typedef VectorObjectsIterator<Polynomial> VectorPolynomialIterator;

}  // namespace WhetStone
}  // namespace Amanzi

#endif

