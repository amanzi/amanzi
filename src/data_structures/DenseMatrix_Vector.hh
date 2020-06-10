/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//! Helper data structure and factory for using WhetStone Tensors

#ifndef AMANZI_DENSEMATRIX_VECTOR_HH_
#define AMANZI_DENSEMATRIX_VECTOR_HH_

#include "AmanziTypes.hh"
#include "AmanziVector.hh"
#include "DenseMatrix.hh"
#include "CSR.hh"

namespace Amanzi {

//
// A simple data structure that keeps a vector of WhetStone Tensors and a
// CompositeVectorSpace to describe the layout/ghosting/entities those tensors
// are associated with.
// -----------------------------------------------------------------------------
struct DenseMatrix_Vector {

  template<typename MEM> 
  using type_t = WhetStone::DenseMatrix<MEM>; 

public: 

  using memory_space = CSR_Matrix::memory_space; 

  DenseMatrix_Vector() : inited(false) {}

  DenseMatrix_Vector(const int& count) : inited(false) {
    prealloc_(count); 
  }

  KOKKOS_INLINE_FUNCTION size_t size() const { return data.size(); }

  void set_shape(const int& i, const int& nrows, const int& ncols) {
    int loc[2] = {nrows, ncols};
    data.set_shape_host(i, loc);
  }

  /**
    * Init the CSR based on a Host vector 
  */
  void Init(const std::vector<WhetStone::DenseMatrix<Kokkos::HostSpace>>& ht){
    // Extract max size 
    const int row_map_size = ht.size(); 
    // Set CSR 
    data.set_size(row_map_size); 
    // Fill row map  and size
    for(int i = 0 ; i < ht.size(); ++i){
      data.set_row_map_host(i,ht[i].Values().size());
      data.set_sizes_host(i,0,ht[i].NumRows());
      data.set_sizes_host(i,1,ht[i].NumCols()); 
    }
    data.prefix_sum(); 
    // Fill entries 
    for(int i = 0 ; i < ht.size(); ++i){
      int idx = data.row_map_host(i);
      for(int j = 0 ; j < ht[i].Values().extent(0); ++j){
        data.set_entries_host(idx+j,ht[i].Values()(j)); 
      }
    }
    data.update_entries_device(); 
    inited = true; 
  }

  void Init() {
    inited = true;
    data.prefix_sum();
  }
  
  // The operator[] return the value on device 
  KOKKOS_INLINE_FUNCTION
  type_t<DeviceOnlyMemorySpace> operator[](const int& i) const {
    assert(inited);
    return std::move(type_t<DeviceOnlyMemorySpace>(data.at(i), data.size(i,0), data.size(i,1)));
  }

  KOKKOS_INLINE_FUNCTION
  type_t<Kokkos::HostSpace> at_host(const int& i) const {
    // FIXME -- not const correct, but to do so needs a const-correct WhetStone::Tensor,
    // e.g. a WhetStone::Tensor that takes a Kokkos::View<const double*> --etc
    return std::move(type_t<Kokkos::HostSpace>(
        data.at_host(i), data.size_host(i,0), data.size_host(i,1)));
  }
  

  KOKKOS_INLINE_FUNCTION
  type_t<DeviceOnlyMemorySpace> at(const int& i) const {
    // FIXME -- not const correct, but to do so needs a const-correct WhetStone::Tensor,
    // e.g. a WhetStone::Tensor that takes a Kokkos::View<const double*> --etc
    return std::move(type_t<DeviceOnlyMemorySpace>(data.at(i), data.size(i,0), data.size(i,1)));
  }

  void update_entries_host(){
    data.update_entries_host(); 
  }

  void update_entries_device(){
    data.update_entries_device();
  }

  
  CSR_Matrix data;
  bool inited;

 private:

  void prealloc_(const int& count) {
    data = std::move(CSR_Matrix(count));
  }

}; 

} // namespace Amanzi


#endif
