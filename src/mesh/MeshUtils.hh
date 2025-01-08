/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Julien Loiseau
*/

/*
  Note: these utilities are used internally to the mesh library and are not
  needed by client code.
*/

#pragma once

#include "AmanziTypes.hh"
#include "ViewUtils.hh"
#include "MeshDefs.hh"

namespace Amanzi {
namespace AmanziMesh {
namespace Impl {

//
// These Getter and RaggedGetter functions are used within MeshCache to create
// a standard way of choosing between access patterns.
//

// -----------------------------------------------------------------------------
// Helper functions and structs for getting things out of MeshCache
// -----------------------------------------------------------------------------

//
// An empty struct for when a compute function is not valid
//
struct NullFunc {};

//
// Struct for selecting between two types, a function pointer and a nullptr,
// based on the memory space.
//
template <MemSpace_kind MEM>
struct ComputeFunction {
  //
  // Returns an initialized function pointer on HOST, NoFunc (e.g. nullptr) on DEVICE
  //
  template <typename CF>
  static decltype(auto) hostOnly()
  {
    if constexpr (MEM == MemSpace_kind::HOST) {
      return CF();
    } else {
      return NullFunc();
    }
  }

  //
  // Returns an initialized function pointer on DEVICE, NoFunc (e.g. nullptr) on HOST
  //
  template <typename CF>
  static decltype(auto) deviceOnly()
  {
    if constexpr (MEM == MemSpace_kind::DEVICE) {
      return CF();
    } else {
      return NullFunc();
    }
  }
};


//
// A generic "recipe" for getting or computing a MeshCache quantity.
//
// A Getter works to provide any single entity or computed/derived quantity
// that may be cached, provided by the framework, or computed, on device or on
// host.  It does so following a generic recipe that can be used by all such
// get{Thing}() methods in MeshCache.
//
// See MeshCache::getCellVolume() for an example.
//
template <MemSpace_kind MEM = MemSpace_kind::HOST,
          AccessPattern_kind AP = AccessPattern_kind::DEFAULT>
struct Getter {
  // DATA: the View-type storing the quantity to be got
  // MF: the MeshFramework-like object
  // FF: the Framework Function -- a lambda that uses the framework, if it is
  //     not-null, to get the quantity.
  // CF: the compute function -- a lambda or NullFunc
  template <typename DATA, typename MF, typename FF, typename CF>
  static KOKKOS_INLINE_FUNCTION decltype(auto)
  get(bool cached, DATA& d, MF& mf, FF&& f, CF&& c, const Entity_ID i)
  {
    using type_t = typename DATA::t_dev::traits::value_type;
    // To avoid the cast to non-reference
    if (cached) return static_cast<type_t>(view<MEM>(d)(i));
    if constexpr (MEM == MemSpace_kind::HOST && (!std::is_same_v<FF, decltype(nullptr)>)) {
      if (mf.get()) return f(i);
    }
    if constexpr (std::is_invocable_v<CF, const Entity_ID>) { return c(i); }
    assert(false && "No access to cache/framework/compute available in Getter");
    return type_t{};
  }
}; // Getter

// specialization for CACHE
template <MemSpace_kind MEM>
struct Getter<MEM, AccessPattern_kind::CACHE> {
  template <typename DATA, typename MF, typename FF, typename CF>
  static KOKKOS_INLINE_FUNCTION decltype(auto)
  get(bool cached, DATA& d, MF&, FF&&, CF&&, const Entity_ID i)
  {
    assert(cached);
    return view<MEM>(d)(i);
  }
}; // Getter

// specialization for FRAMEWORK
template <MemSpace_kind MEM>
struct Getter<MEM, AccessPattern_kind::FRAMEWORK> {
  template <typename DATA, typename MF, typename FF, typename CF>
  static KOKKOS_INLINE_FUNCTION decltype(auto)
  get(bool, DATA&, MF& mf, FF&& f, CF&&, const Entity_ID i)
  {
    static_assert(!std::is_same<FF, decltype(nullptr)>::value);
    static_assert(MEM == MemSpace_kind::HOST);
    assert(mf.get());
    return f(i);
  }
}; // Getter

// specialization for COMPUTE
template <MemSpace_kind MEM>
struct Getter<MEM, AccessPattern_kind::COMPUTE> {
  template <typename DATA, typename MF, typename FF, typename CF>
  static KOKKOS_INLINE_FUNCTION decltype(auto)
  get(bool, DATA&, MF&, FF&&, CF&& c, const Entity_ID i)
  {
    static_assert(std::is_invocable_v<CF, const Entity_ID>);
    return c(i);
  }
}; // Getter


//
// A generic "recipe" for getting or computing a MeshCache quantity.
//
// A RaggedGetter works to provide collection of entities or computed/derived
// quantities that may be cached, provided by the framework, or computed, on
// device or on host.  It does so following a generic recipe that can be used
// by all such get{Thing}() methods in MeshCache.
//
// See MeshCache::getCellFaces() for an example.
//
template <MemSpace_kind MEM = MemSpace_kind::HOST,
          AccessPattern_kind AP = AccessPattern_kind::DEFAULT>
struct RaggedGetter {
  // DATA: the View-type storing the quantity to be got
  // MF: the MeshFramework-like object
  // FF: the Framework Function -- a lambda that uses the framework, if it is
  //     not-null, to get the quantity.
  // CF: the compute function -- a lambda or NullFunc
  template <typename DATA, typename MF, typename FF, typename CF>
  static KOKKOS_INLINE_FUNCTION decltype(auto)
  get(bool cached, DATA& d, MF& mf, FF&& f, CF&& c, const Entity_ID n)
  {
    using view_t = typename decltype(d.template getRow<MEM>(n))::const_type;
    if (cached) { return static_cast<view_t>(d.template getRowUnmanaged<MEM>(n)); }

    if constexpr (MEM == MemSpace_kind::HOST && (!std::is_same_v<FF, decltype(nullptr)>)) {
      if (mf.get()) { return static_cast<view_t>(f(n)); }
    }

    if constexpr (std::is_invocable_v<CF, const Entity_ID>) {
      view_t v = c(n);
      return v;
    }
    assert(false && "No access to cache/framework/compute available in RaggedGetter");
    return view_t{};
  }
};

// specialization for CACHE
template <MemSpace_kind MEM>
struct RaggedGetter<MEM, AccessPattern_kind::CACHE> {
  template <typename DATA, typename MF, typename FF, typename CF>
  static KOKKOS_INLINE_FUNCTION decltype(auto)
  get(bool cached, DATA& d, MF&, FF&&, CF&&, const Entity_ID n)
  {
    assert(cached);
    return d.template getRowUnmanaged<MEM>(n);
  }
};

// specialization for FRAMEWORK
template <MemSpace_kind MEM>
struct RaggedGetter<MEM, AccessPattern_kind::FRAMEWORK> {
  template <typename DATA, typename MF, typename FF, typename CF>
  static KOKKOS_INLINE_FUNCTION decltype(auto)
  get(bool, DATA&, MF& mf, FF&& f, CF&&, const Entity_ID n)
  {
    static_assert(!std::is_same_v<FF, decltype(nullptr)>);
    static_assert(MEM == MemSpace_kind::HOST);
    assert(mf.get());
    return f(n);
  }
};

// specialization for COMPUTE
template <MemSpace_kind MEM>
struct RaggedGetter<MEM, AccessPattern_kind::COMPUTE> {
  template <typename DATA, typename MF, typename FF, typename CF>
  static KOKKOS_INLINE_FUNCTION decltype(auto)
  get(bool, DATA&, MF&, FF&&, CF&& c, const Entity_ID n)
  {
    static_assert(std::is_invocable_v<CF, const Entity_ID>);
    return c(n);
  }
};

} // namespace Impl
} // namespace AmanziMesh
} // namespace Amanzi


namespace Errors {
//
// Exception to be thrown when Framework does not or cannot implement a concept.
//
class FrameworkNotImplemented : public Message {
  using Message::Message;
};

//
// Exception to be thrown when MeshCache cannot find a concept.
//
class MeshNotImplemented : public Message {
  using Message::Message;
};


} // namespace Errors
