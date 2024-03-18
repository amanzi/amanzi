/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

// Reductions for min-loc and max-loc based on those from
// Teuchos_TimeMonitor.cpp

#pragma once

#ifdef KOKKOS_ENABLE_CUDA
#include <cuda/std/limits>
#endif

#include "Teuchos_ReductionOp.hpp"
#include "AmanziVector.hh"

namespace Amanzi {
namespace Reductions {

//
// Used to do min value and its location reductions
//
template <typename Scalar, typename GO>
struct MinLoc {
  Scalar val;
  GO loc;

  KOKKOS_INLINE_FUNCTION
  MinLoc() { init(); }

  KOKKOS_INLINE_FUNCTION
  MinLoc(Scalar s, GO g) : val(s), loc(g) {}

  KOKKOS_INLINE_FUNCTION
  MinLoc(const MinLoc& other) : val(other.val), loc(other.loc) {}

  KOKKOS_INLINE_FUNCTION
  void init()
  {
#ifdef __CUDA_ARCH__
    val = cuda::std::numeric_limits<Scalar>::max();
#else
    val = std::numeric_limits<Scalar>::max();
#endif
    loc = -1;
  }

  KOKKOS_INLINE_FUNCTION
  MinLoc& operator+=(const MinLoc& other)
  {
    if (other.val < val) {
      val = other.val;
      loc = other.loc;
    }
    return *this;
  }

  KOKKOS_INLINE_FUNCTION
  void operator+=(const volatile MinLoc& other) volatile
  {
    if (other.val < val) {
      val = other.val;
      loc = other.loc;
    }
  }
};


//
// Max value and GID location
//
template <typename Scalar, typename GO>
struct MaxLoc {
  Scalar val;
  GO loc;

  KOKKOS_INLINE_FUNCTION
  MaxLoc() { init(); }

  KOKKOS_INLINE_FUNCTION
  MaxLoc(Scalar s, GO g) : val(s), loc(g) {}

  KOKKOS_INLINE_FUNCTION
  MaxLoc(const MaxLoc& other) : val(other.val), loc(other.loc) {}

  KOKKOS_INLINE_FUNCTION
  void init()
  {
#ifdef __CUDA_ARCH__
    val = cuda::std::numeric_limits<Scalar>::min();
#else
    val = std::numeric_limits<Scalar>::min();
#endif
    loc = -1;
  }

  KOKKOS_INLINE_FUNCTION
  MaxLoc& operator+=(const MaxLoc& other)
  {
    if (other.val > val) {
      val = other.val;
      loc = other.loc;
    }
    return *this;
  }

  KOKKOS_INLINE_FUNCTION
  void operator+=(const volatile MaxLoc& other) volatile
  {
    if (other.val > val) {
      val = other.val;
      loc = other.loc;
    }
  }
};


//
// Reduction of array used by Teuchos, which knows how to serialize pairs so we
// convert Reduction_type objects to and from pairs.
//
template <typename Ordinal, typename Scalar, typename Reductor_type>
struct ReductionArray : public Teuchos::ValueTypeReductionOp<Ordinal, std::pair<Scalar, GO>> {
  void reduce(const Ordinal count,
              const std::pair<Scalar, GO> in_buffer[],
              std::pair<Scalar, GO> inout_buffer[]) const
  {
    for (Ordinal ind = 0; ind < count; ++ind) {
      Reductor_type inout(inout_buffer[ind].first, inout_buffer[ind].second);
      inout += Reductor_type(in_buffer[ind].first, in_buffer[ind].second);
      inout_buffer[ind].first = inout.val;
      inout_buffer[ind].second = inout.loc;
    }
  }
};


template <typename Scalar, typename Reductor_type>
Reductor_type
reduceAllLoc(const Comm_type& comm, const Reductor_type& local_val)
{
  using Reduction_type = ReductionArray<int, Scalar, Reductor_type>;
  std::pair<Scalar, GO> global_val_pair;
  std::pair<Scalar, GO> local_val_pair{ local_val.val, local_val.loc };
  Teuchos::reduceAll(comm, Reduction_type(), 1, &local_val_pair, &global_val_pair);
  return Reductor_type(global_val_pair.first, global_val_pair.second);
}


template <typename Scalar, typename Reductor_type>
Reductor_type
reduceAllLoc(const Vector_type_<Scalar>& vec,
             Kokkos::View<const LO*, DefaultMemorySpace>* indices = nullptr)
{
  Reductor_type local_val;
  if (indices) {
    auto vec_view = vec.getLocalViewDevice(Tpetra::Access::ReadOnly);
    // do the local reduction in Kokkos
    auto indices_loc = *indices;
    Kokkos::parallel_reduce(
      "Amanzi::reduceAll",
      indices_loc.extent(0),
      KOKKOS_LAMBDA(const int& i, Reductor_type& lval) {
        Reductor_type val(vec_view(indices_loc(i), 0), indices_loc(i));
        lval += val;
      },
      local_val);
  } else {
    auto vec_view = vec.getLocalViewDevice(Tpetra::Access::ReadOnly);
    // do the local reduction in Kokkos
    Kokkos::parallel_reduce(
      "Amanzi::reduceAll",
      vec_view.extent(0),
      KOKKOS_LAMBDA(const int& c, Reductor_type& lval) {
        Reductor_type val(vec_view(c, 0), c);
        lval += val;
      },
      local_val);
  }
  // now do the global reduction in Teuchos, first converting to GID from LID
  local_val.loc = vec.getMap()->getGlobalElement(local_val.loc);
  return reduceAllLoc<Scalar>(*vec.getMap()->getComm(), local_val);
}

template <typename Scalar>
MinLoc<Scalar, GO>
reduceAllMinLoc(const Vector_type_<Scalar>& vec,
                Kokkos::View<const LO*, DefaultMemorySpace>* indices = nullptr)
{
  return reduceAllLoc<Scalar, MinLoc<Scalar, GO>>(vec);
}

template <typename Scalar>
MaxLoc<Scalar, GO>
reduceAllMaxLoc(const Vector_type_<Scalar>& vec,
                Kokkos::View<const LO*, DefaultMemorySpace>* indices = nullptr)
{
  return reduceAllLoc<Scalar, MaxLoc<Scalar, GO>>(vec);
}

} // namespace Reductions
} // namespace Amanzi
