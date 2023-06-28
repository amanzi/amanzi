/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Julien Loiseau
*/

#pragma once

#include "AmanziTypes.hh"
#include "ViewUtils.hh"
#include "MeshDefs.hh"

namespace Amanzi {
namespace AmanziMesh {

//
// These Getter and RaggedGetter functions are used within MeshCache to create
// a standard way of choosing between access patterns.
//
// Move these to an Impl namespace?
//

template<MemSpace_kind MEM = MemSpace_kind::HOST, AccessPattern_kind AP = AccessPattern_kind::DEFAULT>
struct Getter {
  template<typename DATA, typename MF, typename FF, typename CF>
  static KOKKOS_INLINE_FUNCTION decltype(auto)
  get(bool cached, DATA& d, MF& mf, FF&& f, CF&& c, const Entity_ID i){
      using type_t = typename DATA::t_dev::traits::value_type;
      // To avoid the cast to non-reference
      if (cached) return static_cast<type_t>(view<MEM>(d)(i));
      if constexpr(MEM == MemSpace_kind::HOST){
        if constexpr (!std::is_same_v<FF,decltype(nullptr)>)
          if(mf.get())
            return f(i);
      }
      if constexpr (!std::is_same_v<CF,decltype(nullptr)>)
        return c(i);
      assert(false);
      return type_t{};
  }
}; // Getter

template<MemSpace_kind MEM>
struct Getter<MEM,AccessPattern_kind::CACHE> {
  template<typename DATA, typename MF, typename FF, typename CF>
  static KOKKOS_INLINE_FUNCTION decltype(auto)
  get(bool cached, DATA& d, MF&, FF&&, CF&&, const Entity_ID i){
      assert(cached);
      return view<MEM>(d)(i);
  }
}; // Getter

template<MemSpace_kind MEM>
struct Getter<MEM,AccessPattern_kind::FRAMEWORK> {
  template<typename DATA, typename MF, typename FF, typename CF>
  static KOKKOS_INLINE_FUNCTION decltype(auto)
  get(bool, DATA&, MF& mf, FF&& f, CF&&, const Entity_ID i){
      static_assert(!std::is_same<FF,decltype(nullptr)>::value);
      static_assert(MEM == MemSpace_kind::HOST);
      assert(mf.get());
      return f(i);
  }
}; // Getter


template<MemSpace_kind MEM>
struct Getter<MEM,AccessPattern_kind::COMPUTE> {
  template<typename DATA, typename MF, typename FF, typename CF>
  static KOKKOS_INLINE_FUNCTION decltype(auto)
  get(bool, DATA&, MF&, FF&&, CF&& c, const Entity_ID i){
    static_assert(!std::is_same<CF,decltype(nullptr)>::value);
    return c(i);
  }
}; // Getter


// Getters for raggedViews
template<MemSpace_kind MEM = MemSpace_kind::HOST, AccessPattern_kind AP = AccessPattern_kind::DEFAULT>
struct RaggedGetter{
  template<typename DATA, typename MF, typename FF, typename CFD, typename CF>
  static KOKKOS_INLINE_FUNCTION decltype(auto)
  get (bool cached, DATA& d, MF& mf, FF&& f, CFD&& cd, CF&& c, const Entity_ID n) {
    using view_t = typename decltype(d.template getRow<MEM>(n))::const_type;
    if (cached) {
      return static_cast<view_t>(d.template getRowUnmanaged<MEM>(n));
    }
    if constexpr(MEM == MemSpace_kind::HOST){

      if constexpr (!std::is_same<FF,decltype(nullptr)>::value){
        if(mf.get()){
          return static_cast<view_t>(f(n));
        }
      }
      if constexpr (!std::is_same<CF,decltype(nullptr)>::value){
        return static_cast<view_t>(c(c));
      }
    } else {
      if constexpr (!std::is_same<CFD,decltype(nullptr)>::value) {
        static_cast<view_t>(cd(c));
      }
    }
    assert(false && "No access to cache/framework/compute available");
    return view_t{};
  }
};

template<MemSpace_kind MEM>
struct RaggedGetter<MEM,AccessPattern_kind::CACHE>{
  template<typename DATA, typename MF, typename FF, typename CFD, typename CF>
  static KOKKOS_INLINE_FUNCTION decltype(auto)
  get (bool cached, DATA& d, MF&, FF&&, CFD&&, CF&&, const Entity_ID n) {
    assert(cached);
    return d.template getRowUnmanaged<MEM>(n);
  }
};

template<MemSpace_kind MEM>
struct RaggedGetter<MEM,AccessPattern_kind::FRAMEWORK>{
  template<typename DATA, typename MF, typename FF, typename CFD, typename CF>
  static KOKKOS_INLINE_FUNCTION decltype(auto)
  get (bool, DATA&, MF& mf, FF&& f, CFD&&, CF&&, const Entity_ID n) {
    static_assert(!std::is_same<FF,decltype(nullptr)>::value);
    static_assert(MEM == MemSpace_kind::HOST);
    assert(mf.get());
    return f(n);
  }
};

template<MemSpace_kind MEM>
struct RaggedGetter<MEM,AccessPattern_kind::COMPUTE>{
  template<typename DATA, typename MF, typename FF, typename CFD, typename CF>
  static KOKKOS_INLINE_FUNCTION decltype(auto)
  get (bool, DATA&, MF&, FF&&, CFD&& cd, CF&& c, const Entity_ID n) {
    if constexpr(MEM == MemSpace_kind::HOST){
      static_assert(!std::is_same<CF,decltype(nullptr)>::value);
      return c(n);
    }else{
      static_assert(MEM==MemSpace_kind::DEVICE);
      static_assert(!std::is_same<CFD,decltype(nullptr)>::value);
      return cd(n);
    }
  }
};

} // namespace AmanziMesh
} // namespace Amanzi
