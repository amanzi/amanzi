/*
  Copyright 2010-201x held jointly by LANL, ORNL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
           Julien Loiseau (jloiseau@lanl.gov)
           Rao Garimella (rao@lanl.gov)
*/

//! Implement detail utility functions for doing mesh cache work

#pragma once
#include "MeshDefs.hh"


//
// NOTE: begin/end must live in Kokkos namespace to work!
//
// This simply allows ranged-based for loops on Kokkos Views that are on host.
// We use this a lot...
//
namespace Kokkos {
namespace Impl {

template<typename T, typename ...Args>
class View_iter {
  using iterator_category = std::forward_iterator_tag;
  using difference_type = int;
  using value_type = T;
  using pointer = value_type*;
  using reference = value_type&;

  static_assert(Kokkos::SpaceAccessibility<typename Kokkos::View<T*, Args...>::execution_space,
                typename Kokkos::HostSpace>::accessible);
  using View_type = Kokkos::View<T*, Args...>;

public:
  View_iter(const View_type& v) : View_iter(v,0) {}
  View_iter(const View_type& v, int i) : v_(v), i_(i) {}

  reference operator*() const { return v_(i_); }
  pointer operator->() { return &v_(i_); }

  // prefix
  View_iter& operator++() { i_++; return *this; }
  // postfix
  View_iter operator++(int) { View_iter tmp(*this); ++(*this); return tmp; }

  friend bool operator==(const View_iter& a, const View_iter& b) {
    return a.v_ == b.v_ && a.i_ == b.i_;
  }
  friend bool operator!=(const View_iter& a, const View_iter& b) {
    return !(a == b);
  }

private:
  int i_;
  const View_type& v_;
};

} // namespace Impl

template<typename T, typename ...Args>
Impl::View_iter<T, Args...>
begin(const View<T*, Args...>& view) {
  return Impl::View_iter<T, Args...>(view);
}

template<typename T, typename ...Args>
Impl::View_iter<T, Args...>
end(const View<T*, Args...>& view) {
  return Impl::View_iter<T, Args...>(view, view.size());
}

} // namespace Kokkos



namespace Amanzi {
namespace AmanziMesh {

//
// Utility functions for Kokkos operations
//

//
// Get the right view from a dual view
//
template<MemSpace_type M, typename DualView>
auto& // Kokkos::View of the same type as DV, on M
view(DualView& dv)
{
  if constexpr(M == MemSpace_type::HOST) {
    return dv.h_view;
  } else {
    return dv.d_view;
  }
}


//
// Conversion from view on host to vector
//
template<typename T, typename ...Args>
std::vector<T>
asVector(const Kokkos::View<T*, Args...> view)
{
  static_assert(Kokkos::SpaceAccessibility<typename Kokkos::View<T*, Args...>::execution_space,
                typename Kokkos::HostSpace>::accessible);
  // currently no fix for this -- C++20 will use span
  std::vector<T> vec(view.size());
  for (int i=0; i!=view.size(); ++i) vec.emplace_back(view[i]);
  return vec;
}


//
// Conversion from vector to view on host
//
template<typename T>
Kokkos::View<const T*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
asView(const std::vector<T>& vec)
{
  using View_type = Kokkos::View<const T*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  View_type my_view(vec.data(), vec.size());
  return my_view;
}


//
// Conversion between list and view through deep_copy.
//
template<typename T, typename ...Args>
void my_deep_copy(Kokkos::View<T*, Args...>& out,
                  const std::vector<T>& in)
{
  auto in_view = asView(in);
  Kokkos::deep_copy(out, in_view);
}



//
// Conversion between list and dual view through deep copy and sync.
//
// NOTE: change this to DefaultDevice!
template<typename T, typename MemSpace = Kokkos::HostSpace>
Kokkos::DualView<typename std::remove_const<T>::type*, MemSpace>
asDualView(const std::vector<T>& in)
{
  using DV_type = Kokkos::DualView<typename std::remove_const<T>::type*, MemSpace>;
  DV_type dv;
  my_deep_copy(dv.h_view, in);
  dv.template modify<typename DV_type::host_mirror_space>();
  dv.template sync<typename DV_type::execution_space>();
  return dv;
}


// note, this template is left here despite not being used in case of future
// refactoring for a more general struct.
template<typename T = Entity_ID>
struct RaggedArray_DualView {
  DualView_type<int> rows;
  DualView_type<T> entries;

  using host_mirror_space = typename DualView_type<T>::host_mirror_space;
  using execution_space = typename DualView_type<T>::execution_space;

  template<MemSpace_type MEM>
  KOKKOS_INLINE_FUNCTION
  decltype(auto)
  getRow(int row) {
    return Kokkos::subview(view<MEM>(entries), std::make_pair(view<MEM>(rows)[row], view<MEM>(rows)[row+1]));
  }

  template<MemSpace_type MEM>
  KOKKOS_INLINE_FUNCTION
  T& get(int row, int i) {
    return view<MEM>(entries)[view<MEM>(rows)[row]+i];
  }
};


//
// Cache a RaggedArray from a callable, e.g. getCellFaces()
//
template<typename Func>
RaggedArray_DualView<Entity_ID>
asRaggedArray_DualView(Func mesh_func, Entity_ID count) {
  RaggedArray_DualView<Entity_ID> adj;
  adj.rows.resize(count+1);

  // do a count first, setting rows
  Entity_ID_List ents;
  int total = 0;
  for (Entity_ID i=0; i!=count; ++i) {
    view<MemSpace_type::HOST>(adj.rows)[i] = total;

    mesh_func(i, ents);
    total += ents.size();
  }
  view<MemSpace_type::HOST>(adj.rows)[count] = total;
  adj.entries.resize(total);

  for (Entity_ID i=0; i!=count; ++i) {
    mesh_func(i, ents);
    Kokkos::View<Entity_ID*, Kokkos::DefaultHostExecutionSpace> row_view = adj.getRow<MemSpace_type::HOST>(i);
    my_deep_copy(row_view, ents);
  }

  adj.rows.template modify<typename RaggedArray_DualView<Entity_ID>::host_mirror_space>();
  adj.rows.template sync<typename RaggedArray_DualView<Entity_ID>::execution_space>();
  adj.entries.template modify<typename RaggedArray_DualView<Entity_ID>::host_mirror_space>();
  adj.entries.template sync<typename RaggedArray_DualView<Entity_ID>::execution_space>();
  return adj;
}

} // namespace AmanziMesh
} // namespace Amanzi


