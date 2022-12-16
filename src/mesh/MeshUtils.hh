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

#include "Kokkos_Core.hpp"
#include "Kokkos_DualView.hpp"

enum class MemSpace_type {
  HOST,
  DEVICE
};

using size_type = Kokkos::View<int*, Kokkos::DefaultHostExecutionSpace>::size_type;

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

  using View_type = Kokkos::View<T*, Args...>;

public:
  KOKKOS_INLINE_FUNCTION View_iter(const View_type& v) : View_iter(v,0) {}
  KOKKOS_INLINE_FUNCTION View_iter(const View_type& v, int i) : v_(v), i_(i) {}

  KOKKOS_INLINE_FUNCTION reference operator*() const { return v_(i_); }
  KOKKOS_INLINE_FUNCTION pointer operator->() { return &v_(i_); }

  // prefix
  KOKKOS_INLINE_FUNCTION View_iter& operator++() { i_++; return *this; }
  // postfix
  KOKKOS_INLINE_FUNCTION View_iter operator++(int) { View_iter tmp(*this); ++(*this); return tmp; }

  KOKKOS_INLINE_FUNCTION friend bool operator==(const View_iter& a, const View_iter& b) {
    return a.v_ == b.v_ && a.i_ == b.i_;
  }
  KOKKOS_INLINE_FUNCTION friend bool operator!=(const View_iter& a, const View_iter& b) {
    return !(a == b);
  }

private:
  int i_;
  const View_type& v_;
};

} // namespace Impl

template<typename T, typename ...Args>
Impl::View_iter<T, Args...>
KOKKOS_INLINE_FUNCTION begin(const View<T*, Args...>& view) {
  return Impl::View_iter<T, Args...>(view);
}

template<typename T, typename ...Args>
Impl::View_iter<T, Args...>
KOKKOS_INLINE_FUNCTION end(const View<T*, Args...>& view) {
  return Impl::View_iter<T, Args...>(view, view.size());
}

} // namespace Kokkos



namespace Amanzi {

//
// Utility functions for Kokkos operations
//

//
// Get the right view from a dual view
//
template<MemSpace_type M, typename DualView>
KOKKOS_INLINE_FUNCTION
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
  std::vector<T> vec;
  vec.reserve(view.size());
  for (int i=0; i!=view.size(); ++i) vec.emplace_back(view[i]);
  return vec;
}


//
// Conversion from vector to non-owning view on host.
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
  dv.resize(in.size()); 
  my_deep_copy(dv.h_view, in);
  Kokkos::deep_copy(dv.d_view,dv.h_view); 
  return dv;
}


// note, this template is left here despite not being used in case of future
// refactoring for a more general struct.
template<typename T>
struct RaggedArray_DualView {
  using type_t = T; 

  template<MemSpace_type MEM>
  using constview = 
    Kokkos::View<const T*,
    std::conditional<
      MEM==MemSpace_type::DEVICE,
      Kokkos::DefaultExecutionSpace,
      Kokkos::HostSpace>>; 

  Kokkos::DualView<int*> rows;
  Kokkos::DualView<T*> entries;

  using host_mirror_space = typename Kokkos::DualView<T*>::host_mirror_space;
  using execution_space = typename Kokkos::DualView<T*>::execution_space;

  RaggedArray_DualView() {}


  RaggedArray_DualView(const std::vector<std::vector<T>>& vect){ 

    rows.resize(vect.size()+1); 
    int count = 0; 
    view<MemSpace_type::HOST>(rows)[0] = count; 
    for(int i = 1 ; i < vect.size()+1; ++i){
      count += vect[i-1].size(); 
      view<MemSpace_type::HOST>(rows)[i] = count; 
    }

    entries.resize(count);

    int cur = 0; 
    for(int i = 0; i < vect.size(); ++i){
      for(int j = 0 ; j < vect[i].size(); ++j){
        view<MemSpace_type::HOST>(entries)[cur++] = vect[i][j]; 
      }
    } 
    update<MemSpace_type::DEVICE>(); 
  }

  template<MemSpace_type MEM>
  KOKKOS_INLINE_FUNCTION
  decltype(auto)
  getRow(int row) {
    return Kokkos::subview(view<MEM>(entries), Kokkos::make_pair(view<MEM>(rows)[row], view<MEM>(rows)[row+1]));
  }

  template<MemSpace_type MEM>
  KOKKOS_INLINE_FUNCTION
  decltype(auto)
  getRow(int row) const {
    return Kokkos::subview(view<MEM>(entries), Kokkos::make_pair(view<MEM>(rows)[row], view<MEM>(rows)[row+1]));
  }

  template<MemSpace_type MEM>
  KOKKOS_INLINE_FUNCTION
  T& get(int row, int i) {
    return view<MEM>(entries)[view<MEM>(rows)[row]+i];
  }

  template<MemSpace_type MEM>
  KOKKOS_INLINE_FUNCTION
  const T& get(int row, int i) const {
    return view<MEM>(entries)[view<MEM>(rows)[row]+i];
  }

  template<MemSpace_type MEM>
  KOKKOS_INLINE_FUNCTION
  int size() const {return view<MEM>(rows).size()-1; }

  template<MemSpace_type MEM>
  KOKKOS_INLINE_FUNCTION
  int size(int row) const { return view<MEM>(rows)[row+1] - view<MEM>(rows)[row]; }

  template<MemSpace_type MEM> 
  void update(){ 
    if constexpr (MEM == MemSpace_type::HOST){
      Kokkos::deep_copy(rows.view_host(),rows.view_device()); 
      Kokkos::deep_copy(entries.view_host(),entries.view_device()); 
    }else{
      Kokkos::deep_copy(rows.view_device(),rows.view_host()); 
      Kokkos::deep_copy(entries.view_device(),entries.view_host()); 
    }
  }

  bool operator!=(const RaggedArray_DualView& oth){
    return oth.rows.view_host().size() != rows.view_host().size() 
    && oth.entries.view_host().size() != entries.view_host().size(); 
  }

  void resize(int s1, int s2) {
    rows.resize(s1);
    entries.resize(s1*s2);  
    rows.h_view[0] = 0; 
    for(int i = 1 ; i < s1; ++i){
      rows.h_view[i] = rows.h_view[i-1]+s2; 
    }
    update<MemSpace_type::DEVICE>(); 
  } 

};


//
// Cache a RaggedArray from a callable, e.g. getCellFaces()
//
template<typename T, typename Func>
RaggedArray_DualView<T>
asRaggedArray_DualView(Func mesh_func, int count) {
  RaggedArray_DualView<T> adj;
  adj.rows.resize(count+1);

  // do a count first, setting rows
  std::vector<T> ents;
  int total = 0;
  for (int i=0; i!=count; ++i) {
    view<MemSpace_type::HOST>(adj.rows)[i] = total;

    mesh_func(i, ents);
    total += ents.size();
  }
  view<MemSpace_type::HOST>(adj.rows)[count] = total;
  adj.entries.resize(total);

  for (int i=0; i!=count; ++i) {
    mesh_func(i, ents);
    Kokkos::View<T*, Kokkos::DefaultHostExecutionSpace> row_view = adj.template getRow<MemSpace_type::HOST>(i);
    my_deep_copy(row_view, ents);
  }
  Kokkos::deep_copy(adj.rows.view_device(), adj.rows.view_host()); 
  Kokkos::deep_copy(adj.entries.view_device(),adj.entries.view_host()); 
  return adj;
}

template<typename T, typename List> 
KOKKOS_INLINE_FUNCTION
bool is_present(const T& v, const List& l){ 
  for(int i = 0 ; i < l.size(); ++i){ 
    if(v == l[i])
      return true; 
  }
  return false; 
}

// Find the right number of threads 
template<typename T>
constexpr int ThreadsPerTeams(){ 
  #ifdef KOKKOS_ENABLE_CUDA
  if constexpr (std::is_same_v<T,Kokkos::Cuda>){
    return 32;
  }else // (std::is_same_v<T,Kokkos::Serial>){
  #endif
  { 
    return 1;
  }
}
} // namespace Amanzi