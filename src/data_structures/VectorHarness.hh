/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

//!

#ifndef VECTOR_HARNESS_HPP_
#define VECTOR_HARNESS_HPP_

#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_Vector.hpp>

namespace Amanzi {
namespace VectorHarness {

enum AccessMode { ReadOnly, WriteOnly, ReadWrite };

namespace Impl {
// Given a global object, get its default memory space (both the
// type and the default instance thereof).
template <class GlobalObjectType>
struct DeviceOnlyMemorySpace {
  using type = typename GlobalObjectType::device_type::memory_space;

  // Given a global object, get its (default) memory space instance.
  static type space(const GlobalObjectType& /* G */)
  {
    // This stub just assumes that 'type' is default constructible.
    // In Kokkos, default-constructing a memory space instance just
    // gives the default memory space.
    return type();
  }
};

// Struct that tells withLocalAccess how to access a global object's
// local data.  Do not use this directly; start with readOnly,
// writeOnly, or readWrite.
template <class GlobalObjectType, class MemorySpace,
          const AccessMode am>
class LocalAccess; // forward declaration

// Mapping from LocalAccess to the "master" local object type.  The
// latter gets the local data from a global object, and holds on to
// it until after the user's function (input to withLocalAccess)
// returns.
template <class LocalAccessType>
struct GetMasterLocalObject {};

// Given a LocalAccess instance (which has a reference to a global
// object), get an instance of its master local object.  This may be
// a heavyweight operation.
//
// If you need to specialize this, just specialize get() in
// GetMasterLocalObject above.
template <class LocalAccessType>
typename GetMasterLocalObject<LocalAccessType>::master_local_object_type
getMasterLocalObject(LocalAccessType LA)
{
  return GetMasterLocalObject<LocalAccessType>::get(LA);
}

// Mapping from "master" local object type to the nonowning "local
// view" type that users see (as arguments to the function that they
// give to withLocalAccess).  The master local object may encode the
// memory space and access mode, but the mapping to local view type
// may also need run-time information.
template <class MasterLocalObjectType>
struct GetNonowningLocalObject {};

// Given a master local object, get an instance of a nonowning local
// object.  Users only ever see the nonowning local object, and
// subviews (slices) thereof.  This is supposed to be a lightweight
// operation.
//
// If you need to specialize this, just specialize get() in
// GetNonowningLocalObject above.
template <class MasterLocalObjectType>
typename GetNonowningLocalObject<
  MasterLocalObjectType>::nonowning_local_object_type
getNonowningLocalObject(const MasterLocalObjectType& master)
{
  return GetNonowningLocalObject<MasterLocalObjectType>::get(master);
}

// Use the LocalAccess type as the template parameter to determine
// the type of the nonowning local view to the global object's data.
// This only works if GetMasterLocalObject has been specialized for
// these template parameters, and if GetNonowningLocalObject has
// been specialized for the resulting "master" local object type.
template <class LocalAccessType>
class LocalAccessFunctionArgument {
 private:
  using gmlo = GetMasterLocalObject<LocalAccessType>;
  using master_local_object_type = typename gmlo::master_local_object_type;
  using gnlo = GetNonowningLocalObject<master_local_object_type>;

 public:
  using type = typename gnlo::nonowning_local_object_type;
};
} // namespace Impl

//////////////////////////////////////////////////////////////////////
// Users call readOnly, writeOnly, and readWrite, in order to declare
// how they intend to access a global object's local data.
//////////////////////////////////////////////////////////////////////

// Declare that you want to access the given global object's local
// data in read-only mode.
template <class GlobalObjectType>
Impl::LocalAccess<GlobalObjectType,
                  typename Impl::DeviceOnlyMemorySpace<GlobalObjectType>::type,
                  ReadOnly>
readOnly(GlobalObjectType&);

// Declare that you want to access the given global object's local
// data in write-only mode.
template <class GlobalObjectType>
Impl::LocalAccess<GlobalObjectType,
                  typename Impl::DeviceOnlyMemorySpace<GlobalObjectType>::type,
                  WriteOnly>
writeOnly(GlobalObjectType&);

// Declare that you want to access the given global object's local
// data in read-and-write mode.
template <class GlobalObjectType>
Impl::LocalAccess<GlobalObjectType,
                  typename Impl::DeviceOnlyMemorySpace<GlobalObjectType>::type,
                  ReadWrite>
readWrite(GlobalObjectType&);

namespace Impl {
// Declaration of access intent for a global object.
//
// Users aren't supposed to make instances of this class.  They
// should use readOnly, writeOnly, or readWrite instead, then call
// methods like on() and valid() on the resulting LocalAccess
// instance.
template <class GlobalObjectType, class MemorySpace, const AccessMode am>
class LocalAccess {
 public:
  using global_object_type = GlobalObjectType;
  using memory_space = typename MemorySpace::memory_space;
  static constexpr AccessMode access_mode = am;

  // Users must NOT call the LocalAccess constructor directly.  They
  // should instead start by calling readOnly, writeOnly, or
  // readWrite above.  They may then use methods like on() or
  // valid() (see below).
  //
  // G is a reference, because we only access it in a delimited
  // scope.  G is nonconst, because even read-only local access may
  // modify G.  For example, G may give access to its local data via
  // lazy allocation of a data structure that differs from its
  // normal internal storage format.
  //
  // Memory spaces should behave like Kokkos memory spaces.  Default
  // construction should work and should get the default instance of
  // the space.  Otherwise, it may make sense to get the default
  // memory space from G.
  LocalAccess(global_object_type& G, memory_space space = memory_space(),
              const bool isValid = true)
    : G_(G), space_(space), valid_(isValid)
  {}

  // Type users see, that's an argument to the function that they give
  // to withLocalAccess.
  using function_argument_type = typename LocalAccessFunctionArgument<
    LocalAccess<global_object_type, memory_space, access_mode>>::type;

 public:
  // Give users run-time control over whether they actually want to
  // access the object.  If isValid is false, implementations should
  // not spend any effort getting the master local object.  This may
  // save time on allocating temporary space, copying from device to
  // host, etc.  This implies that implementations must be able to
  // construct "null" / empty master local objects.
  LocalAccess<GlobalObjectType, MemorySpace, am> valid(const bool isValid) const
  {
    // std::cout << "  .valid(" << (isValid ? "true" : "false") << ")"
    //           << std::endl;
    return { this->G_, this->space_, isValid };
  }

  // Let users access this object in a different memory space.
  //
  // NOTE: This also works for PGAS.  'space' in that case could be
  // something like an MPI_Win, or a UPC-style "shared pointer."
  //
  // NewMemorySpace should behave like a Kokkos memory space
  // instance.
  template <class NewMemorySpace>
  LocalAccess<GlobalObjectType, NewMemorySpace, am>
  on(NewMemorySpace space) const
  {
    return { this->G_, space, this->valid_ };
  }

  // Is access supposed to be valid?  (See valid() above.)
  bool isValid() const { return this->valid_; }

  memory_space getSpace() const { return space_; }

 public:
  // Keep by reference, because this struct is only valid in a
  // delimited scope.
  global_object_type& G_;
  memory_space space_; // assume shallow-copy semantics
  bool valid_;         // will I actually need to access this object?

 private:
  // Nonmember "constructors"; see above for declarations.  This are
  // friends, because they are the only ways that users are supposed
  // to construct LocalAccess instances.
  template <class GOT>
  friend LocalAccess<GOT, typename Impl::DeviceOnlyMemorySpace<GOT>::type,
                     ReadOnly>
  readOnly(GOT&);
  template <class GOT>
  friend LocalAccess<GOT, typename Impl::DeviceOnlyMemorySpace<GOT>::type,
                     WriteOnly>
  writeOnly(GOT&);
  template <class GOT>
  friend LocalAccess<GOT, typename Impl::DeviceOnlyMemorySpace<GOT>::type,
                     ReadWrite>
  readWrite(GOT&);
};
} // namespace Impl

template <class GOT>
Impl::LocalAccess<GOT, typename Impl::DeviceOnlyMemorySpace<GOT>::type, ReadOnly>
readOnly(GOT& G)
{
  // std::cout << "readOnly" << std::endl;
  return { G, Impl::DeviceOnlyMemorySpace<GOT>::space(G), true };
}

template <class GOT>
Impl::LocalAccess<GOT, typename Impl::DeviceOnlyMemorySpace<GOT>::type, WriteOnly>
writeOnly(GOT& G)
{
  // std::cout << "writeOnly" << std::endl;
  return { G, Impl::DeviceOnlyMemorySpace<GOT>::space(G), true };
}

template <class GOT>
Impl::LocalAccess<GOT, typename Impl::DeviceOnlyMemorySpace<GOT>::type, ReadWrite>
readWrite(GOT& G)
{
  // std::cout << "readWrite" << std::endl;
  return { G, Impl::DeviceOnlyMemorySpace<GOT>::space(G), true };
}


namespace Impl {
template <class SC, class LO, class GO, class NT, class MemorySpace>
using multivector_nonconst_nonowning_local_object_type =
  typename std::conditional<
    std::is_same<typename MemorySpace::memory_space, Kokkos::HostSpace>::value,
    typename Tpetra::MultiVector<SC, LO, GO, NT>::dual_view_type::t_host,
    typename Tpetra::MultiVector<SC, LO, GO, NT>::dual_view_type::t_dev>::type;

template <class SC, class LO, class GO, class NT, class MemorySpace>
using multivector_const_nonowning_local_object_type =
  typename multivector_nonconst_nonowning_local_object_type<
    SC, LO, GO, NT, MemorySpace>::const_type;

template <class SC, class LO, class GO, class NT, class MemorySpace,
          const AccessMode am>
using multivector_nonowning_local_object_type = typename std::conditional<
  am == ReadOnly,
  multivector_const_nonowning_local_object_type<SC, LO, GO, NT, MemorySpace>,
  multivector_nonconst_nonowning_local_object_type<SC, LO, GO, NT,
                                                   MemorySpace>>::type;
} // namespace Impl

template <class SC, class LO, class GO, class NT, class MemorySpace,
          const AccessMode access_mode>
Impl::multivector_nonowning_local_object_type<SC, LO, GO, NT, MemorySpace,
                                              access_mode>
getMultiVector(Impl::LocalAccess<Tpetra::MultiVector<SC, LO, GO, NT>,
                                 MemorySpace, access_mode>
                 LA)
{
  using memory_space = typename MemorySpace::memory_space;
  using ret_type = Impl::multivector_nonowning_local_object_type<SC,
                                                                 LO,
                                                                 GO,
                                                                 NT,
                                                                 memory_space,
                                                                 access_mode>;

  // if (access_mode == WriteOnly) { ...} // FIXME (mfh 22 Oct 2018)
  if (LA.G_.template need_sync<memory_space>()) {
    LA.G_.template sync<memory_space>();
  }
  if (access_mode != ReadWrite) { LA.G_.template modify<memory_space>(); }
  // FIXME (mfh 22 Oct 2018) This might break if we need copy-back
  // semantics, e.g., for a memory space for which the
  // Tpetra::MultiVector does not store data.  In that case, we
  // would need some kind of object whose destructor copies back,
  // and it would need to have the whole DualView, not just the
  // View on one side.  Watch out for copy elision.  The object
  // could just be std::shared_ptr and could handle copy-back via
  // custom deleter.
  if (LA.isValid()) {
    // this converts to const if applicable
    return ret_type(LA.G_.template getLocalView<memory_space>());
  } else {
    return ret_type(); // "null" Kokkos::View
  }
}

template <class SC, class LO, class GO, class NT, class MemorySpace,
          const AccessMode access_mode>
Impl::multivector_nonowning_local_object_type<SC, LO, GO, NT, MemorySpace,
                                              access_mode>
getMultiVector(Impl::LocalAccess<const Tpetra::MultiVector<SC, LO, GO, NT>,
                                 MemorySpace, access_mode>
                 LA)
{
  // Here is more minor evil.  See above.
  Tpetra::MultiVector<SC, LO, GO, NT>& G =
    const_cast<Tpetra::MultiVector<SC, LO, GO, NT>&>(LA.G_);
  Impl::
    LocalAccess<Tpetra::MultiVector<SC, LO, GO, NT>, MemorySpace, access_mode>
      LA2(G, LA.space_, LA.valid_);
  return getMultiVector(LA2);
}

template <class SC, class LO, class GO, class NT, class MemorySpace,
          const AccessMode access_mode>
auto
getVector(Impl::LocalAccess<Tpetra::MultiVector<SC, LO, GO, NT>, MemorySpace,
                            access_mode>
            LA,
          const int whichColumn = 0)
  -> decltype(Kokkos::subview(getMultiVector(LA), Kokkos::ALL(), whichColumn))
{
  return Kokkos::subview(getMultiVector(LA), Kokkos::ALL(), whichColumn);
}

template <class SC, class LO, class GO, class NT, class MemorySpace,
          const AccessMode access_mode>
auto
getVector(Impl::LocalAccess<const Tpetra::MultiVector<SC, LO, GO, NT>,
                            MemorySpace, access_mode>
            LA,
          const int whichColumn = 0)
  -> decltype(Kokkos::subview(getMultiVector(LA), Kokkos::ALL(), whichColumn))
{
  using nc_global_object_type = Tpetra::MultiVector<SC, LO, GO, NT>;
  nc_global_object_type& G = const_cast<nc_global_object_type&>(LA.G_);
  Impl::LocalAccess<nc_global_object_type, MemorySpace, access_mode> LA2(
    G, LA.space_, LA.valid_);
  return getVector(LA2, whichColumn);
}

template <class SC, class LO, class GO, class NT, class MemorySpace,
          const AccessMode access_mode>
auto
getVector(
  Impl::LocalAccess<Tpetra::Vector<SC, LO, GO, NT>, MemorySpace, access_mode>
    LA)
  -> decltype(Kokkos::subview(
    getMultiVector(Impl::LocalAccess<Tpetra::MultiVector<SC, LO, GO, NT>,
                                     MemorySpace, access_mode>(LA.G_, LA.space_,
                                                               LA.valid_)),
    Kokkos::ALL(), 0))
{
  return Kokkos::subview(
    getMultiVector(Impl::LocalAccess<Tpetra::MultiVector<SC, LO, GO, NT>,
                                     MemorySpace,
                                     access_mode>(LA.G_, LA.space_, LA.valid_)),
    Kokkos::ALL(),
    0);
}

template <class SC, class LO, class GO, class NT, class MemorySpace,
          const AccessMode access_mode>
auto
getVector(Impl::LocalAccess<const Tpetra::Vector<SC, LO, GO, NT>, MemorySpace,
                            access_mode>
            LA)
  -> decltype(getVector(
    Impl::LocalAccess<Tpetra::Vector<SC, LO, GO, NT>, MemorySpace, access_mode>(
      const_cast<Tpetra::Vector<SC, LO, GO, NT>&>(LA.G_), LA.space_,
      LA.valid_)))
{
  return getVector(
    Impl::LocalAccess<Tpetra::Vector<SC, LO, GO, NT>, MemorySpace, access_mode>(
      const_cast<Tpetra::Vector<SC, LO, GO, NT>&>(LA.G_),
      LA.space_,
      LA.valid_));
}


namespace Impl {
// Specialization of GetMasterLocalObject for Tpetra::MultiVector.
template <class SC, class LO, class GO, class NT, class MemorySpace,
          const AccessMode am>
struct GetMasterLocalObject<
  LocalAccess<Tpetra::MultiVector<SC, LO, GO, NT>, MemorySpace, am>> {
 public:
  using local_access_type =
    LocalAccess<Tpetra::MultiVector<SC, LO, GO, NT>, MemorySpace, am>;

 private:
  using global_object_type = typename local_access_type::global_object_type;
  using memory_space = typename local_access_type::memory_space;
  static constexpr AccessMode access_mode = local_access_type::access_mode;
  using non_const_value_type =
    typename Tpetra::MultiVector<SC, LO, GO, NT>::impl_scalar_type;
  using value_type = typename std::conditional<access_mode == ReadOnly,
                                               const non_const_value_type,
                                               non_const_value_type>::type;

 public:
  // FIXME (mfh 22 Oct 2018) See FIXME below.
  using master_local_object_type = Kokkos::View<
    value_type**,
    typename global_object_type::dual_view_type::t_dev::array_layout,
    MemorySpace>; // FIXME (mfh 22 Oct 2018) need to make sure execution_space
                  // matches

  static master_local_object_type get(local_access_type LA)
  {
    // std::cout << "Get master local object" << std::endl;

    // if (access_mode == WriteOnly) { ...} // FIXME (mfh 22 Oct 2018)
    if (LA.G_.template need_sync<memory_space>()) {
      LA.G_.template sync<memory_space>();
    }
    if (access_mode != ReadWrite) { LA.G_.template modify<memory_space>(); }
    // FIXME (mfh 22 Oct 2018) This might break if we need copy-back
    // semantics, e.g., for a memory space for which the
    // Tpetra::MultiVector does not store data.  In that case, we
    // would need some kind of object whose destructor copies back,
    // and it would need to have the whole DualView, not just the
    // View on one side.  Watch out for copy elision.  The object
    // could just be std::shared_ptr and could handle copy-back via
    // custom deleter.
    if (LA.isValid()) {
      // this converts to const if applicable
      return master_local_object_type(
        LA.G_.template getLocalView<memory_space>());
    } else {
      return master_local_object_type(); // "null" Kokkos::View
    }
  }
};

// Specialization of GetMasterLocalObject for const Tpetra::MultiVector.
template <class SC, class LO, class GO, class NT, class MemorySpace,
          const AccessMode am>
struct GetMasterLocalObject<
  LocalAccess<const Tpetra::MultiVector<SC, LO, GO, NT>, MemorySpace, am>> {
 public:
  using local_access_type =
    LocalAccess<const Tpetra::MultiVector<SC, LO, GO, NT>, MemorySpace, am>;

 private:
  using global_object_type = typename std::remove_const<
    typename local_access_type::global_object_type>::type;
  using memory_space = typename local_access_type::memory_space;
  static constexpr AccessMode access_mode = local_access_type::access_mode;
  using non_const_value_type =
    typename Tpetra::MultiVector<SC, LO, GO, NT>::impl_scalar_type;
  using value_type = typename std::conditional<access_mode == ReadOnly,
                                               const non_const_value_type,
                                               non_const_value_type>::type;

 public:
  // FIXME (mfh 22 Oct 2018) See FIXME below.
  using master_local_object_type = Kokkos::View<
    value_type**,
    typename global_object_type::dual_view_type::t_dev::array_layout,
    MemorySpace>; // FIXME (mfh 22 Oct 2018) need to make sure execution_space
                  // matches

  static master_local_object_type get(local_access_type LA)
  {
    // std::cout << "Get master local object" << std::endl;

    // This is a bit evil, but fixing this in a less evil way is hard.
    global_object_type& G = const_cast<global_object_type&>(LA.G_);

    // if (access_mode == WriteOnly) { ...} // FIXME (mfh 22 Oct 2018)
    if (G.template need_sync<memory_space>()) {
      G.template sync<memory_space>();
    }
    if (access_mode != ReadWrite) { G.template modify<memory_space>(); }
    // FIXME (mfh 22 Oct 2018) This might break if we need copy-back
    // semantics, e.g., for a memory space for which the
    // Tpetra::MultiVector does not store data.  In that case, we
    // would need some kind of object whose destructor copies back,
    // and it would need to have the whole DualView, not just the
    // View on one side.  Watch out for copy elision.  The object
    // could just be std::shared_ptr and could handle copy-back via
    // custom deleter.
    if (LA.isValid()) {
      // this converts to const if applicable
      return master_local_object_type(G.template getLocalView<memory_space>());
    } else {
      return master_local_object_type(); // "null" Kokkos::View
    }
  }
};

// Specialization of GetMasterLocalObject for Tpetra::Vector.
template <class SC, class LO, class GO, class NT, class MemorySpace,
          const AccessMode am>
struct GetMasterLocalObject<
  LocalAccess<Tpetra::Vector<SC, LO, GO, NT>, MemorySpace, am>> {
 public:
  using local_access_type =
    LocalAccess<Tpetra::Vector<SC, LO, GO, NT>, MemorySpace, am>;

 private:
  using global_object_type = typename local_access_type::global_object_type;
  using memory_space = typename local_access_type::memory_space;
  static constexpr AccessMode access_mode = local_access_type::access_mode;
  using non_const_value_type = typename global_object_type::impl_scalar_type;
  using value_type = typename std::conditional<access_mode == ReadOnly,
                                               const non_const_value_type,
                                               non_const_value_type>::type;

 public:
  // FIXME (mfh 22 Oct 2018) See FIXME below.
  using master_local_object_type = Kokkos::View<
    value_type**,
    typename global_object_type::dual_view_type::t_dev::array_layout,
    MemorySpace>; // FIXME (mfh 22 Oct 2018) need to make sure execution_space
                  // matches

  static master_local_object_type get(local_access_type LA)
  {
    std::cout << "Get master local object" << std::endl;

    // if (access_mode == WriteOnly) { ...} // FIXME (mfh 22 Oct 2018)
    if (LA.G_.template need_sync<memory_space>()) {
      LA.G_.template sync<memory_space>();
    }
    if (access_mode != ReadWrite) { LA.G_.template modify<memory_space>(); }
    // FIXME (mfh 22 Oct 2018) This might break if we need copy-back
    // semantics, e.g., for a memory space for which the
    // Tpetra::MultiVector does not store data.  In that case, we
    // would need some kind of object whose destructor copies back,
    // and it would need to have the whole DualView, not just the
    // View on one side.  Watch out for copy elision.  The object
    // could just be std::shared_ptr and could handle copy-back via
    // custom deleter.
    if (LA.isValid()) {
      // this converts to const if applicable
      return master_local_object_type(
        LA.G_.template getLocalView<memory_space>());
    } else {
      return master_local_object_type(); // "null" Kokkos::View
    }
  }
};

// Specialization of GetNonowningLocalObject for StubMasterLocalObject.
template <class ValueType, class LayoutType, class MemorySpace>
struct GetNonowningLocalObject<
  Kokkos::View<ValueType**, LayoutType, MemorySpace>> {
  using master_local_object_type =
    Kokkos::View<ValueType**, LayoutType, MemorySpace>;
  using nonowning_local_object_type =
    Kokkos::View<ValueType**, LayoutType, MemorySpace, Kokkos::MemoryUnmanaged>;
  static nonowning_local_object_type get(const master_local_object_type& M)
  {
    return nonowning_local_object_type(M); // standard Kokkos::View assignment
  }
};

} // namespace Impl

} // namespace VectorHarness
} // namespace Amanzi

#endif // VECTOR_HARNESS_HPP
