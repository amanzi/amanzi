/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

//! Helper data structure and factory for using WhetStone Tensors
#ifndef AMANZI_DENSEVECTOR_VECTOR_HH_
#define AMANZI_DENSEVECTOR_VECTOR_HH_

#include "AmanziTypes.hh"
#include "AmanziVector.hh"
#include "DenseVector.hh"
#include "CSR.hh"

namespace Amanzi {

template<class MemorySpace>
struct DenseVector_VectorData : CSR<double, 1, MemorySpace> {
  template <typename MEM>
  using Value_type = WhetStone::DenseVector<MEM>;

}


//
// A simple data structure that keeps a vector of WhetStone Tensors and a
// CompositeVectorSpace to describe the layout/ghosting/entities those tensors
// are associated with.
// -----------------------------------------------------------------------------
struct DenseVector_Vector {
  template <typename MEM>
  using Value_type = WhetStone::DenseVector<MEM>;

 public:
  using memory_space = CSR_Vector::memory_space;

  DenseVector_Vector() : inited(false) {}

  DenseVector_Vector(const int& count) : inited(false) { prealloc_(count); }

  KOKKOS_INLINE_FUNCTION size_t size() const { return data.size(); }

  void set_shape(const int& i, const int& size)
  {
    int loc[1] = { size };
    data.set_shape_host(i, loc);
  }

  void Init()
  {
    inited = true;
    data.prefix_sum();
  }

  void putScalar(double val) { data.putScalar(val); }

  // The operator[] return the value on device
  KOKKOS_INLINE_FUNCTION
  Value_type<DeviceOnlyMemorySpace> operator[](const int& i) const
  {
    assert(inited);
    return std::move(Value_type<DeviceOnlyMemorySpace>(data.at(i)));
  }

  KOKKOS_INLINE_FUNCTION
  Value_type<Kokkos::HostSpace> at_host(const int& i) const
  {
    // FIXME -- not const correct, but to do so needs a const-correct
    // WhetStone::Tensor, e.g. a WhetStone::Tensor that takes a
    // Kokkos::View<const double*> --etc
    return std::move(Value_type<Kokkos::HostSpace>(data.at_host(i)));
  }


  KOKKOS_INLINE_FUNCTION
  Value_type<DeviceOnlyMemorySpace> at(const int& i) const
  {
    // FIXME -- not const correct, but to do so needs a const-correct
    // WhetStone::Tensor, e.g. a WhetStone::Tensor that takes a
    // Kokkos::View<const double*> --etc
    return std::move(Value_type<DeviceOnlyMemorySpace>(data.at(i)));
  }

  CSR_Vector data;
  bool inited;

 private:
  void prealloc_(const int& count) { data = std::move(CSR_Vector(count)); }
};

} // namespace Amanzi


#endif
