/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*
  Utils

  Helper classes which define const iterators that point to const.
*/

#pragma once

#include "boost/iterator_adaptors.hpp"
#include "Teuchos_RCP.hpp"
#include "Kokkos_Core.hpp"

//
// Iterators on Kokkos views
//
// This simply allows ranged-based for loops on Kokkos Views that are on host.
// We use this a lot...
//
namespace Kokkos {
namespace Impl {

template<typename View_type>
struct View_iter {
  using iterator_category = std::forward_iterator_tag;
  using value_type = typename View_type::value_type;
  using difference_type = int;
  using pointer = value_type*;
  using reference = value_type&;

  KOKKOS_INLINE_FUNCTION View_iter(const View_type& v) : View_iter(v,0) {}
  KOKKOS_INLINE_FUNCTION View_iter(const View_type& v, int i) : v_(v), i_(i) {}

  KOKKOS_INLINE_FUNCTION reference operator*() const { return v_(i_); }
  KOKKOS_INLINE_FUNCTION pointer operator->() { return &v_(i_); }

  // prefix
  KOKKOS_INLINE_FUNCTION View_iter& operator++() { i_++; return *this; }
  KOKKOS_INLINE_FUNCTION View_iter& operator--() { i_--; return *this; }
  // postfix
  KOKKOS_INLINE_FUNCTION View_iter operator++(int) { View_iter tmp(*this); ++(*this); return tmp; }

  KOKKOS_INLINE_FUNCTION friend View_iter operator+(const View_iter& v, const int& d) {
    View_iter tmp(v); tmp+=d; return tmp;
  }
  KOKKOS_INLINE_FUNCTION friend View_iter operator+(const int& d, const View_iter& v) {
    return v + d;
  }
  KOKKOS_INLINE_FUNCTION friend int operator-(const View_iter& l, const View_iter& r){
    return l.i_-r.i_;
  }
  KOKKOS_INLINE_FUNCTION friend bool operator==(const View_iter& a, const View_iter& b) {
    return a.v_ == b.v_ && a.i_ == b.i_;
  }
  KOKKOS_INLINE_FUNCTION friend bool operator!=(const View_iter& a, const View_iter& b) {
    return !(a == b);
  }
  KOKKOS_INLINE_FUNCTION friend bool operator<(const View_iter& l, const View_iter& r) {
    return l.v_ == r.v_ && l.i_ < r.i_;
  }
  KOKKOS_INLINE_FUNCTION friend bool operator<=(const View_iter& l, const View_iter& r) {
    return l.v_ == r.v_ && l.i_ <= r.i_;
  }
  KOKKOS_INLINE_FUNCTION friend bool operator>(const View_iter& l, const View_iter& r) {
    return l.v_ == r.v_ && l.i_ > r.i_;
  }
  KOKKOS_INLINE_FUNCTION friend bool operator>=(const View_iter& l, const View_iter& r) {
    return l.v_ == r.v_ && l.i_ >= r.i_;
  }
  KOKKOS_INLINE_FUNCTION View_iter& operator+=(const int& incr){
    this->i_ += incr;
    return *this;
  }
  KOKKOS_INLINE_FUNCTION View_iter& operator-=(const int& decr){
    this->i_ -= decr;
    return *this;
  }
  KOKKOS_INLINE_FUNCTION View_iter operator-(const int& decr){
    this->i_ -= decr;
    return *this;
  }

private:
  int i_;
  View_type v_;
};

} // namespace Impl

template<typename View_type>
Impl::View_iter<View_type>
KOKKOS_INLINE_FUNCTION begin(const View_type& view) {
  return Impl::View_iter<View_type>(view);
}

template<typename View_type>
Impl::View_iter<View_type>
KOKKOS_INLINE_FUNCTION end(const View_type& view) {
  return Impl::View_iter<View_type>(view, view.size());
}

} // namespace Kokkos


//
// Iterators of const from const containers
//
namespace Amanzi {
namespace Utils {

template <class Container_t, class Contained_t>
struct iterator
  : boost::iterator_adaptor<iterator<Container_t, Contained_t>, typename Container_t::iterator> {
  iterator() {}
  iterator(typename Container_t::iterator i)
    : boost::iterator_adaptor<iterator<Container_t, Contained_t>, typename Container_t::iterator>(i)
  {}
};

template <typename Container_t, typename Contained_t>
struct const_iterator : boost::iterator_adaptor<const_iterator<Container_t, Contained_t>,
                                                typename Container_t::const_iterator,
                                                Teuchos::RCP<const Contained_t>,
                                                boost::use_default,
                                                Teuchos::RCP<const Contained_t>> {
  const_iterator() {}
  const_iterator(iterator<Container_t, Contained_t> i)
    : boost::iterator_adaptor<const_iterator<Container_t, Contained_t>,
                              typename Container_t::const_iterator,
                              Teuchos::RCP<const Contained_t>,
                              boost::use_default,
                              Teuchos::RCP<const Contained_t>>(i.base())
  {}
  const_iterator(typename Container_t::const_iterator i)
    : boost::iterator_adaptor<const_iterator<Container_t, Contained_t>,
                              typename Container_t::const_iterator,
                              Teuchos::RCP<const Contained_t>,
                              boost::use_default,
                              Teuchos::RCP<const Contained_t>>(i)
  {}
};

} // namespace Utils
} // namespace Amanzi


