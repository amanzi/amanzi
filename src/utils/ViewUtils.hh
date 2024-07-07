/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Julien Loiseau (jloiseau@lanl.gov)
*/

#pragma once
#include <set>
#include <vector>

#include "AmanziTypes.hh"
#include "MeshView.hh"

//
// Utility functions for Kokkos View operations
//
namespace Amanzi {

//
// Get the right view from a dual view
//
template <MemSpace_kind M, typename DualView>
KOKKOS_INLINE_FUNCTION auto& // Kokkos::MeshView of the same type as DV, on M
view(DualView& dv)
{
  if constexpr (M == MemSpace_kind::HOST) {
    return dv.h_view;
  } else {
    return dv.d_view;
  }
}


//
// Conversion from view on host to vector
//
template <typename T, typename... Args>
std::vector<std::remove_cv_t<T>>
asVector(const Kokkos::MeshView<T*, Args...> view)
{
  static_assert(Kokkos::SpaceAccessibility<typename Kokkos::MeshView<T*, Args...>::execution_space,
                                           typename Kokkos::HostSpace>::accessible);
  // currently no fix for this -- C++20 will use span
  std::vector<std::remove_cv_t<T>> vec;
  vec.reserve(view.size());
  for (int i = 0; i != view.size(); ++i) vec.emplace_back(view[i]);
  return vec;
}


//
// Conversion from vector to non-owning view on host.
//
template <typename T>
Kokkos::MeshView<const T*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
toNonOwningView(const std::vector<T>& vec)
{
  using View_type =
    Kokkos::MeshView<const T*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  View_type my_view(vec.data(), vec.size());
  return my_view;
}

template <typename V, typename T>
void
vectorToConstView(V& view, const std::vector<T> vec)
{
  static_assert(Kokkos::SpaceAccessibility<typename V::execution_space,
                                           typename Kokkos::HostSpace>::accessible);
  Kokkos::MeshView<T*, typename V::execution_space> lview;
  Kokkos::resize(lview, vec.size());
  for (int i = 0; i < lview.size(); ++i) lview[i] = vec[i];
  view = std::move(lview);
}

template <typename V, typename T>
void
vectorToView(V& view, const std::vector<T> vec)
{
  static_assert(Kokkos::SpaceAccessibility<typename V::execution_space,
                                           typename Kokkos::HostSpace>::accessible);
  Kokkos::resize(view, vec.size());
  for (int i = 0; i < view.size(); ++i) view[i] = vec[i];
}

template <typename V, typename T>
void
setToView(V& view, const std::set<T> vec)
{
  static_assert(Kokkos::SpaceAccessibility<typename V::execution_space,
                                           typename Kokkos::HostSpace>::accessible);

  Kokkos::resize(view, vec.size());
  int i = 0;
  for (const auto& v : vec) view[i++] = v;
}

//
// Conversion between list and dual view through deep copy and sync.
//
template <typename T, typename MemSpace = DefaultMemorySpace>
Kokkos::MeshDualView<typename std::remove_const<T>::type*, MemSpace>
asDualView(std::string name, const std::vector<T>& in)
{
  using DV_type = Kokkos::MeshDualView<typename std::remove_const<T>::type*, MemSpace>;
  DV_type dv(name, in.size());
  Kokkos::deep_copy(dv.h_view, toNonOwningView(in));
  Kokkos::deep_copy(dv.d_view, dv.h_view);
  return dv;
}


//
// Conversion between view and dual view through deep copy and sync.
//
template <typename T, typename MemSpace = DefaultMemorySpace, typename... Args>
Kokkos::MeshDualView<typename std::remove_const<T>::type*, MemSpace>
asDualView(const Kokkos::MeshView<T*, Args...>& in)
{
  Kokkos::MeshDualView<typename std::remove_const<T>::type*, MemSpace> dv(in.label(), in.size());
  Kokkos::deep_copy(dv.h_view, in);
  Kokkos::deep_copy(dv.d_view, dv.h_view);
  return dv;
}

//
// Conversion between view and dual view through deep copy and sync.
//
template <typename T, typename MemSpace = DefaultMemorySpace, typename... Args>
Kokkos::DualView<typename std::remove_const<T>::type*, MemSpace>
asDualView(const Kokkos::View<T*, Args...>& in)
{
  Kokkos::DualView<typename std::remove_const<T>::type*, MemSpace> dv(in.label(), in.size());
  Kokkos::deep_copy(dv.h_view, in);
  Kokkos::deep_copy(dv.d_view, dv.h_view);
  return dv;
}


template <typename T, typename MemSpace = DefaultMemorySpace, typename... Args>
Kokkos::DualView<typename std::remove_const<T>::type**, MemSpace>
asDualView(const Kokkos::View<T**, Args...>& in)
{
  Kokkos::DualView<typename std::remove_const<T>::type**, MemSpace> dv(in.label(), in.extent(0), in.extent(1));
  Kokkos::deep_copy(dv.h_view, in);
  Kokkos::deep_copy(dv.d_view, dv.h_view);
  return dv;
}


namespace Impl {

template <class T>
struct function_traits : function_traits<decltype(&T::operator())> {};

template <class T, class C>
struct function_traits<T C::*> : function_traits<T> {
  using class_type = C;
};

template <class R, class... PP>
struct function_traits<R(PP...)> {
  using return_type = R;
  using parameter_types = std::tuple<PP...>;
};

template <class R, class... PP>
struct function_traits<R(PP...) const> : function_traits<R(PP...)> {};

// In C++23, these can be merged into the above:
template <class R, class... PP>
struct function_traits<R(PP...) noexcept> : function_traits<R(PP...)> {};

template <class R, class... PP>
struct function_traits<R(PP...) const noexcept> : function_traits<R(PP...)> {};

template <class V>
using dual_view_t = Kokkos::MeshDualView<typename V::data_type, Kokkos::DefaultHostExecutionSpace>;

template <typename Func, typename F, typename... PP>
auto
asDualViews(Func& mesh_func, int count, std::tuple<F, PP...>*)
{
  std::tuple<dual_view_t<PP>...> dvs;
  std::apply([count](auto&... dv) { (dv.resize(count), ...); }, dvs);
  for (int i = 0; i < count; ++i) std::apply([&](auto&... dv) { mesh_func(i, dv.h_view...); }, dvs);
  std::apply([](auto&... dv) { (Kokkos::deep_copy(dv.d_view, dv.h_view), ...); }, dvs);
  return dvs;
}

} // namespace Impl


//
// Accepts a function that returns a tuple, returns a tuple of two dual views
// where each entry in the tuple goes in a different dual view.
//
template <typename Func>
auto
asDualView(Func&& mesh_func, int count)
{
  return Impl::asDualViews(
    mesh_func,
    count,
    static_cast<typename Impl::function_traits<std::remove_reference_t<Func>>::parameter_types*>(
      nullptr));
}


// contains
template <typename T, typename List>
KOKKOS_INLINE_FUNCTION bool
is_present(const T& v, const List& l)
{
  for (int i = 0; i < l.size(); ++i) {
    if (v == l[i]) return true;
  }
  return false;
}


// Create a non-const view from a const view
// In the old format, using vector, a copy was always created
// and thus the cache was never directly modified.
// Using the view interface, a const view is returned from the cache
// and the behavior of using this view directly is made impossible.
template <typename View_type>
auto
alloc_and_deep_copy(const View_type& const_view)
{
  auto non_const_view =
    Kokkos::MeshView<typename View_type::traits::non_const_data_type>("", const_view.size());
  Kokkos::deep_copy(non_const_view, const_view);
  return non_const_view;
}


template<class View_type>
KOKKOS_INLINE_FUNCTION
int
find(const View_type& view, typename View_type::const_value_type& val) {
  for (int i = 0; i != view.size(); ++i) {
    if (view(i) == val) return i;
  }
  return -1;
}




} // namespace Amanzi
