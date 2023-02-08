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

#include <set>
#include "Kokkos_Core.hpp"
#include "Kokkos_DualView.hpp"

enum class MemSpace_kind {
  HOST,
  DEVICE
};

//
// NOTE: begin/end must live in Kokkos namespace to work!
//
// This simply allows ranged-based for loops on Kokkos Views that are on host.
// We use this a lot...
//
namespace Kokkos {

template<class DataType, class... Properties>
struct MeshView; 

namespace Impl {

template<typename T, typename ...Args>
struct View_iter {
  using iterator_category = std::forward_iterator_tag;
  using value_type = T;
  using difference_type = int;
  using pointer = value_type*;
  using reference = value_type&;
  using View_type = MeshView<T*, Args...>;

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

template<class DataType, class... Properties>
struct MeshView: public Kokkos::View<DataType, Properties...>{

  using baseView = Kokkos::View<DataType, Properties...>; 
  using baseView::baseView; 
  using iterator = Impl::View_iter<std::remove_pointer_t<DataType>, Properties...>; 
  using const_iterator = Impl::View_iter<const std::remove_pointer_t<DataType>, Properties...>; 
  using traits = ViewTraits<DataType, Properties...>;

  MeshView(const baseView& bv): baseView(bv) {}
  MeshView(const MeshView& bv): baseView(bv) {}
  
  KOKKOS_FUNCTION MeshView& operator=(const MeshView& other){
    baseView::operator=(other);
    return *this; 
  }

  KOKKOS_FUNCTION MeshView& operator=(const MeshView&& other){
    baseView::operator=(other);
    return *this; 
  }
 
  using HostMirror =
    MeshView<typename traits::non_const_data_type, typename traits::array_layout,
              Device<DefaultHostExecutionSpace,
              typename traits::host_mirror_space::memory_space>>;

  KOKKOS_INLINE_FUNCTION iterator begin() const {
    return iterator(*this);
  }

  KOKKOS_INLINE_FUNCTION iterator end() const {
    return iterator(*this, this->size());
  }

  KOKKOS_INLINE_FUNCTION const_iterator cbegin() const {
    return const_iterator(*this);
  }

  KOKKOS_INLINE_FUNCTION const_iterator cend() const {
    return const_iterator(*this, this->size());
  }

  void insert(iterator v0_e, iterator v1_b, iterator v1_e) {
    //assert(v0_e - *this->end() != 0 && "Only insert at end supported for MeshViews"); 
    std::size_t size = v1_e-v1_b; 
    std::size_t csize = this->size(); 
    Kokkos::resize(*this, size+csize); 
    for(int i = csize; i < this->size(); ++i, ++v1_b) this->operator[](i) = *(v1_b); 
  }

}; 

template<class SubDataType, class... SubProperties, class... Args>
MeshView<SubDataType, SubProperties...> subview(Kokkos::MeshView<SubDataType, SubProperties...> v, Args... args){
  MeshView ret(Kokkos::subview((typename Kokkos::MeshView<SubDataType, SubProperties...>::baseView)v, std::forward<Args>(args)...));
  return ret;
}

template<class DataType, class Arg1Type = void, class Arg2Type = void, class Arg3Type = void>
struct MeshDualView: public Kokkos::ViewTraits<DataType, Arg1Type, Arg2Type, Arg3Type>{

  using traits = Kokkos::ViewTraits<DataType, Arg1Type, Arg2Type, Arg3Type>; 
  using host_mirror_space = typename traits::host_mirror_space; 
  using t_dev = MeshView<typename traits::data_type, Arg2Type, Arg2Type, Arg3Type>; 
  using t_host = typename t_dev::HostMirror; 
  using t_modified_flags = MeshView<unsigned int[2], LayoutLeft, Kokkos::HostSpace>;
  t_modified_flags modified_flags;

  KOKKOS_INLINE_FUNCTION t_host view_host() const {return h_view;}
  KOKKOS_INLINE_FUNCTION t_dev view_device() const {return d_view;}

  t_dev d_view; 
  t_host h_view;

  MeshDualView() = default;
  MeshDualView(const std::string& label,
           const size_t n0 = KOKKOS_IMPL_CTOR_DEFAULT_ARG)
      : modified_flags(t_modified_flags("DualView::modified_flags")),
        d_view(label, n0),
        h_view(create_mirror_view(d_view))  // without UVM, host View mirrors
  {}
  MeshDualView(const MeshDualView& src): d_view(src.d_view), h_view(src.h_view), modified_flags(src.modified_flags){}
  template <class SS, class LS, class DS, class MS>
  MeshDualView(const MeshDualView<SS, LS, DS, MS>& src)
      : modified_flags(src.modified_flags),
        d_view(src.d_view),
        h_view(src.h_view) {}

  void resize(const size_t n0 = KOKKOS_IMPL_CTOR_DEFAULT_ARG) {
    if (modified_flags.data() == nullptr) {
      modified_flags = t_modified_flags("DualView::modified_flags");
    }
    if (modified_flags(1) >= modified_flags(0)) {
      /* Resize on Device */
      ::Kokkos::resize(d_view, n0);
      h_view = create_mirror_view(d_view);

      /* Mark Device copy as modified */
      modified_flags(1) = modified_flags(1) + 1;

    } else {
      /* Realloc on Device */

      ::Kokkos::realloc(d_view, n0);

      const bool sizeMismatch = (h_view.extent(0) != n0) ;
      if (sizeMismatch)
        ::Kokkos::resize(h_view, n0);

      t_host temp_view = create_mirror_view(d_view);

      /* Remap on Host */
      Kokkos::deep_copy(temp_view, h_view);

      h_view = temp_view;

      d_view = create_mirror_view(typename t_dev::execution_space(), h_view);

      /* Mark Host copy as modified */
      modified_flags(0) = modified_flags(0) + 1;
    }
  }

}; 


} // namespace Kokkos

using size_type = Kokkos::MeshView<int*, Kokkos::DefaultHostExecutionSpace>::size_type;

namespace Amanzi {

//
// Utility functions for Kokkos operations
//

//
// Get the right view from a dual view
//
template<MemSpace_kind M, typename DualView>
KOKKOS_INLINE_FUNCTION
auto& // Kokkos::MeshView of the same type as DV, on M
view(DualView& dv)
{
  if constexpr(M == MemSpace_kind::HOST) {
    return dv.h_view;
  } else {
    return dv.d_view;
  }
}

template<typename View, typename T>
KOKKOS_INLINE_FUNCTION
void initView(View& v,const T& t){ 
  for(auto& vv: v) vv = t; 
}

//
// Conversion from view on host to vector
//
template<typename T, typename ...Args>
std::vector<std::remove_cv_t<T>>
asVector(const Kokkos::MeshView<T*, Args...> view)
{
  static_assert(Kokkos::SpaceAccessibility<typename Kokkos::MeshView<T*, Args...>::execution_space,
                typename Kokkos::HostSpace>::accessible);
  // currently no fix for this -- C++20 will use span
  std::vector<std::remove_cv_t<T>> vec;
  vec.reserve(view.size());
  for (int i=0; i!=view.size(); ++i) vec.emplace_back(view[i]);
  return vec;
}


//
// Conversion from vector to non-owning view on host.
//
template<typename T>
Kokkos::MeshView<const T*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
asView(const std::vector<T>& vec)
{
  using View_type = Kokkos::MeshView<const T*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  View_type my_view(vec.data(), vec.size());
  return my_view;
}

template<typename V, typename T>
void vectorToView(V& view, const std::vector<T> vec){
  Kokkos::resize(view,vec.size());
  for(int i = 0; i < view.size(); ++i)
    view[i] = vec[i]; 
}

template<typename V, typename T>
void setToView(V& view, const std::set<T> vec){
  Kokkos::resize(view,vec.size());
  int i = 0; 
  for(const auto& v: vec)
    view[i++] = v; 
}

//
// Conversion between list and dual view through deep copy and sync.
//
// NOTE: change this to DefaultDevice!
template<typename T, typename MemSpace = Kokkos::HostSpace>
Kokkos::MeshDualView<typename std::remove_const<T>::type*, MemSpace>
asDualView(const std::vector<T>& in)
{
  using DV_type = Kokkos::MeshDualView<typename std::remove_const<T>::type*, MemSpace>;
  DV_type dv;
  dv.resize(in.size()); 
  Kokkos::deep_copy(dv.h_view, asView(in));
  Kokkos::deep_copy(dv.d_view,dv.h_view); 
  return dv;
}

template<typename T, typename MemSpace = Kokkos::HostSpace, typename ...Args>
Kokkos::MeshDualView<typename std::remove_const<T>::type*, MemSpace>
asDualView(const Kokkos::MeshView<T*, Args...>& in)
{
  Kokkos::MeshDualView<typename std::remove_const<T>::type*, MemSpace> dv;
  dv.resize(in.size()); 
  Kokkos::deep_copy(dv.h_view,in); 
  Kokkos::deep_copy(dv.d_view,dv.h_view); 
  return dv;
}

// note, this template is left here despite not being used in case of future
// refactoring for a more general struct.
template<typename T>
struct RaggedArray_DualView {

  using View_type = Kokkos::MeshView<T*, Kokkos::DefaultHostExecutionSpace>;
  using Entity_ID_View = View_type;

  using type_t = T; 

  template<MemSpace_kind MEM>
  using constview = 
    Kokkos::MeshView<const T*,
    std::conditional<
      MEM==MemSpace_kind::DEVICE,
      Kokkos::DefaultExecutionSpace,
      Kokkos::HostSpace>>; 

  Kokkos::MeshDualView<int*> rows;
  Kokkos::MeshDualView<T*> entries;

  using host_mirror_space = typename Kokkos::MeshDualView<T*>::host_mirror_space;
  using execution_space = typename Kokkos::MeshDualView<T*>::execution_space;

  RaggedArray_DualView() {}


  RaggedArray_DualView(const std::vector<std::vector<T>>& vect){ 

    rows.resize(vect.size()+1); 
    int count = 0; 
    view<MemSpace_kind::HOST>(rows)[0] = count; 
    for(int i = 1 ; i < vect.size()+1; ++i){
      count += vect[i-1].size(); 
      view<MemSpace_kind::HOST>(rows)[i] = count; 
    }

    entries.resize(count);

    int cur = 0; 
    for(int i = 0; i < vect.size(); ++i){
      for(int j = 0 ; j < vect[i].size(); ++j){
        view<MemSpace_kind::HOST>(entries)[cur++] = vect[i][j]; 
      }
    } 
    update<MemSpace_kind::DEVICE>(); 
  }

  template<MemSpace_kind MEM>
  KOKKOS_INLINE_FUNCTION
  decltype(auto)
  getRow(int row) {
    return Kokkos::subview(view<MEM>(entries), Kokkos::make_pair(view<MEM>(rows)[row], view<MEM>(rows)[row+1]));
  }

  template<MemSpace_kind MEM>
  KOKKOS_INLINE_FUNCTION
  decltype(auto)
  getRow(int row) const {
    return Kokkos::subview(view<MEM>(entries), Kokkos::make_pair(view<MEM>(rows)[row], view<MEM>(rows)[row+1]));
  }

  template<MemSpace_kind MEM>
  KOKKOS_INLINE_FUNCTION
  T& get(int row, int i) {
    return view<MEM>(entries)[view<MEM>(rows)[row]+i];
  }

  template<MemSpace_kind MEM>
  KOKKOS_INLINE_FUNCTION
  const T& get(int row, int i) const {
    return view<MEM>(entries)[view<MEM>(rows)[row]+i];
  }

  template<MemSpace_kind MEM>
  KOKKOS_INLINE_FUNCTION
  int size() const {return view<MEM>(rows).size()-1; }

  template<MemSpace_kind MEM>
  KOKKOS_INLINE_FUNCTION
  int size(int row) const { return view<MEM>(rows)[row+1] - view<MEM>(rows)[row]; }

  template<MemSpace_kind MEM> 
  void update(){ 
    if constexpr (MEM == MemSpace_kind::HOST){
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
    update<MemSpace_kind::DEVICE>(); 
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
  Kokkos::MeshView<T*, Kokkos::DefaultHostExecutionSpace> ents;
  int total = 0;
  for (int i=0; i!=count; ++i) {
    view<MemSpace_kind::HOST>(adj.rows)[i] = total;

    mesh_func(i, ents);
    total += ents.size();
  }
  view<MemSpace_kind::HOST>(adj.rows)[count] = total;
  adj.entries.resize(total);

  for (int i=0; i!=count; ++i) {
    mesh_func(i, ents);
    Kokkos::MeshView<T*, Kokkos::DefaultHostExecutionSpace> row_view = adj.template getRow<MemSpace_kind::HOST>(i);
    Kokkos::resize(row_view, ents.size()); 
    Kokkos::deep_copy(row_view, ents); 
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
