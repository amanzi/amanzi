/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//! Helper data structure and factory for using WhetStone Tensors

#ifndef AMANZI_DENSEVECTOR_VECTOR_HH_
#define AMANZI_DENSEVECTOR_VECTOR_HH_

#include "AmanziTypes.hh"
#include "AmanziVector.hh"
#include "DenseVector.hh"
#include "CSR.hh"

namespace Amanzi {

//
// A simple data structure that keeps a vector of WhetStone Tensors and a
// CompositeVectorSpace to describe the layout/ghosting/entities those tensors
// are associated with.
// -----------------------------------------------------------------------------
struct DenseVector_Vector {

  template<typename MEM> 
  using type_t = WhetStone::DenseVector<MEM>; 

public: 

  using memory_space = CSR_Vector::memory_space; 

  DenseVector_Vector() : inited(false) {}

  DenseVector_Vector(const int& count)  : inited(false) {
    prealloc_(count);
  }

  KOKKOS_INLINE_FUNCTION size_t size() const { return data.size(); }

  void set_shape(const int& i, const int& size) {
    int loc[1] = {size};
    data.set_shape_host(i, loc);
  }

  void Init() {
    inited = true;
    data.prefix_sum();
  }
  
  // The operator[] return the value on device 
  KOKKOS_INLINE_FUNCTION
  type_t<DeviceOnlyMemorySpace> operator[](const int& i) const {
    assert(inited);
    return std::move(type_t<DeviceOnlyMemorySpace>(data.at(i), data.size(i,0)));
  }

  KOKKOS_INLINE_FUNCTION
  type_t<Kokkos::HostSpace> at_host(const int& i) const {
    // FIXME -- not const correct, but to do so needs a const-correct WhetStone::Tensor,
    // e.g. a WhetStone::Tensor that takes a Kokkos::View<const double*> --etc
    return std::move(type_t<Kokkos::HostSpace>(
        data.at_host(i), data.size_host(i,0)));
  }
  

  KOKKOS_INLINE_FUNCTION
  type_t<DeviceOnlyMemorySpace> at(const int& i) const {
    // FIXME -- not const correct, but to do so needs a const-correct WhetStone::Tensor,
    // e.g. a WhetStone::Tensor that takes a Kokkos::View<const double*> --etc
    return std::move(type_t<DeviceOnlyMemorySpace>(data.at(i), data.size(i,0)));
  }
  
  CSR_Vector data;
  bool inited;

 private:

  void prealloc_(const int& count) {
    data = std::move(CSR_Vector(count));
  }

}; 

} // namespace Amanzi


#endif
