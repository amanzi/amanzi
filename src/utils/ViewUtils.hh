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

template <typename View, typename T>
KOKKOS_INLINE_FUNCTION void
initView(View& v, const T& t)
{
  for (auto& vv : v) vv = t;
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
  Kokkos::MeshView<T*> lview;
  Kokkos::resize(lview, vec.size());
  for (int i = 0; i < lview.size(); ++i) lview[i] = vec[i];
  view = lview;
}

template <typename V, typename T>
void
vectorToView(V& view, const std::vector<T> vec)
{
  Kokkos::resize(view, vec.size());
  for (int i = 0; i < view.size(); ++i) view[i] = vec[i];
}

template <typename V, typename T>
void
setToView(V& view, const std::set<T> vec)
{
  Kokkos::resize(view, vec.size());
  int i = 0;
  for (const auto& v : vec) view[i++] = v;
}

//
// Conversion between list and dual view through deep copy and sync.
//
// NOTE: change this to DefaultDevice!
template <typename T, typename MemSpace = Kokkos::HostSpace>
Kokkos::MeshDualView<typename std::remove_const<T>::type*, MemSpace>
asDualView(const std::vector<T>& in)
{
  using DV_type = Kokkos::MeshDualView<typename std::remove_const<T>::type*, MemSpace>;
  DV_type dv;
  dv.resize(in.size());
  Kokkos::deep_copy(dv.h_view, toNonOwningView(in));
  Kokkos::deep_copy(dv.d_view, dv.h_view);
  return dv;
}

template <typename T, typename MemSpace = Kokkos::HostSpace, typename... Args>
Kokkos::MeshDualView<typename std::remove_const<T>::type*, MemSpace>
asDualView(const Kokkos::MeshView<T*, Args...>& in)
{
  Kokkos::MeshDualView<typename std::remove_const<T>::type*, MemSpace> dv;
  dv.resize(in.size());
  Kokkos::deep_copy(dv.h_view, in);
  Kokkos::deep_copy(dv.d_view, dv.h_view);
  return dv;
}

template <class...>
struct types {};
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

namespace detail {
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
} // namespace detail

template <typename Func>
auto
asDualView(Func&& mesh_func, int count)
{
  return detail::asDualViews(
    mesh_func,
    count,
    static_cast<typename function_traits<std::remove_reference_t<Func>>::parameter_types*>(
      nullptr));
}


// note, this template is left here despite not being used in case of future
// refactoring for a more general struct.
template <typename T>
struct RaggedArray_DualView {
  using View_type = Kokkos::MeshView<T*, Kokkos::DefaultHostExecutionSpace>;
  using Entity_ID_View = View_type;

  using type_t = T;

  template <MemSpace_kind MEM>
  using constview = Kokkos::MeshView<const T*, MemoryLocation<MEM>>;

  Kokkos::MeshDualView<int*> rows;
  Kokkos::MeshDualView<T*> entries;

  using host_mirror_space = typename Kokkos::MeshDualView<T*>::host_mirror_space;
  using execution_space = typename Kokkos::MeshDualView<T*>::execution_space;

  RaggedArray_DualView() {}


  RaggedArray_DualView(const std::vector<std::vector<T>>& vect)
  {
    rows.resize(vect.size() + 1);
    int count = 0;
    view<MemSpace_kind::HOST>(rows)[0] = count;
    for (int i = 1; i < vect.size() + 1; ++i) {
      count += vect[i - 1].size();
      view<MemSpace_kind::HOST>(rows)[i] = count;
    }

    entries.resize(count);

    int cur = 0;
    for (int i = 0; i < vect.size(); ++i) {
      for (int j = 0; j < vect[i].size(); ++j) {
        view<MemSpace_kind::HOST>(entries)[cur++] = vect[i][j];
      }
    }
    update<MemSpace_kind::DEVICE>();
  }

  template <MemSpace_kind MEM>
  KOKKOS_INLINE_FUNCTION auto getRowUnmanaged(int row)
  {
    const std::size_t size = view<MEM>(rows)[row + 1] - view<MEM>(rows)[row];
    const auto ptr = view<MEM>(entries).data() + view<MEM>(rows)[row];
    return Kokkos::MeshView<T*, MemoryLocation<MEM>>(ptr, size);
  }

  template <MemSpace_kind MEM>
  KOKKOS_INLINE_FUNCTION auto getRowUnmanaged(int row) const
  {
    const std::size_t size = view<MEM>(rows)[row + 1] - view<MEM>(rows)[row];
    const auto ptr = view<MEM>(entries).data() + view<MEM>(rows)[row];
    return Kokkos::MeshView<T*, MemoryLocation<MEM>>(ptr, size);
  }

  template <MemSpace_kind MEM>
  KOKKOS_INLINE_FUNCTION auto getRow(int row)
  {
    return Kokkos::subview(view<MEM>(entries),
                           Kokkos::make_pair(view<MEM>(rows)[row], view<MEM>(rows)[row + 1]));
  }

  template <MemSpace_kind MEM>
  KOKKOS_INLINE_FUNCTION auto getRow(int row) const
  {
    return Kokkos::subview(view<MEM>(entries),
                           Kokkos::make_pair(view<MEM>(rows)[row], view<MEM>(rows)[row + 1]));
  }

  template <MemSpace_kind MEM>
  KOKKOS_INLINE_FUNCTION T& get(int row, int i)
  {
    return view<MEM>(entries)[view<MEM>(rows)[row] + i];
  }

  template <MemSpace_kind MEM>
  KOKKOS_INLINE_FUNCTION const T& get(int row, int i) const
  {
    return view<MEM>(entries)[view<MEM>(rows)[row] + i];
  }

  template <MemSpace_kind MEM>
  KOKKOS_INLINE_FUNCTION int size() const
  {
    return view<MEM>(rows).size() - 1;
  }

  template <MemSpace_kind MEM>
  KOKKOS_INLINE_FUNCTION int size(int row) const
  {
    return view<MEM>(rows)[row + 1] - view<MEM>(rows)[row];
  }

  template <MemSpace_kind MEM>
  void update()
  {
    if constexpr (MEM == MemSpace_kind::HOST) {
      Kokkos::deep_copy(rows.view_host(), rows.view_device());
      Kokkos::deep_copy(entries.view_host(), entries.view_device());
    } else {
      Kokkos::deep_copy(rows.view_device(), rows.view_host());
      Kokkos::deep_copy(entries.view_device(), entries.view_host());
    }
  }

  bool operator!=(const RaggedArray_DualView& oth)
  {
    return oth.rows.view_host().size() != rows.view_host().size() &&
           oth.entries.view_host().size() != entries.view_host().size();
  }

  void resize(int s1, int s2)
  {
    rows.resize(s1 + 1);
    entries.resize(s1 * s2);
    rows.h_view[0] = 0;
    for (int i = 1; i < s1 + 1; ++i) { rows.h_view[i] = rows.h_view[i - 1] + s2; }
    update<MemSpace_kind::DEVICE>();
  }
};

//
// Cache a RaggedArray from a callable, e.g. getCellFaces()
//
template <typename T, typename Func>
auto
asRaggedArray_DualView(Func mesh_func, int count)
{
  RaggedArray_DualView<T> adj;
  adj.rows.resize(count + 1);

  // do a count first, setting rows
  std::vector<Kokkos::MeshView<const T*, Kokkos::HostSpace>> ents(count);
  int total = 0;
  for (int i = 0; i != count; ++i) {
    view<MemSpace_kind::HOST>(adj.rows)[i] = total;

    mesh_func(i, ents[i]);
    total += ents[i].size();
  }
  view<MemSpace_kind::HOST>(adj.rows)[count] = total;
  adj.entries.resize(total);

  for (int i = 0; i != count; ++i) {
    auto row_view = adj.template getRowUnmanaged<MemSpace_kind::HOST>(i);
    assert(row_view.extent(0) == ents[i].size());
    Kokkos::deep_copy(row_view, ents[i]);
  }
  Kokkos::deep_copy(adj.rows.view_device(), adj.rows.view_host());
  Kokkos::deep_copy(adj.entries.view_device(), adj.entries.view_host());
  return adj;
}

template <typename T, typename List>
KOKKOS_INLINE_FUNCTION bool
is_present(const T& v, const List& l)
{
  for (int i = 0; i < l.size(); ++i) {
    if (v == l[i]) return true;
  }
  return false;
}

// Find the right number of threads
template <typename T>
constexpr int
ThreadsPerTeams()
{
#ifdef KOKKOS_ENABLE_CUDA
  if constexpr (std::is_same_v<T, Kokkos::Cuda>) {
    return 32;
  } else // (std::is_same_v<T,Kokkos::Serial>){
#endif
  {
    return 1;
  }
}


} // namespace Amanzi
