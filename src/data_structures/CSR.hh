/*
  Copyright 2010-202x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Julien Loiseau (jloiseau@lanl.gov)
*/

#pragma once 

#include "Kokkos_Core.hpp"

#include "AmanziTypes.hh"


namespace Amanzi {

/**
 * @brief CSR data structure 
 * @param T Type used for the entries_ 
 * @param D Dimension = number of dimensions allocated in sizes
 * This is used to reconstruct the Matrix/Tensor on the device
 */
template<typename T, int D, class MEMSPACE = DefaultMemorySpace> 
class CSR{

static constexpr int dim = D;

public: 

  CSR() {}

  CSR(const int& row_map_size, const int& entries_size)
  {
    Kokkos::resize(row_map_,row_map_size+1); 
    Kokkos::resize(sizes_,row_map_size,dim);
    Kokkos::resize(entries_,entries_size);  
  }

  CSR(const int& row_map_size)
  {
    Kokkos::resize(row_map_,row_map_size+1); 
    Kokkos::resize(sizes_,row_map_size,dim);
  }

  /**
   * @brief Compute the prefix sum of the row_map_ 
   * The CSR is not usable before this operation
   */ 
  void prefix_sum(){
    int tmp1 = row_map_(0);
    row_map_(0) = 0;
    for(int i = 0 ; i < row_map_.extent(0); ++i){ 
      int tmp2 = row_map_(i+1);
      row_map_(i+1) = tmp1 + row_map_(i);
      tmp1 = tmp2;   
    }
    if (entries_.extent(0) < row_map_(row_map_.extent(0)-1))
      Kokkos::resize(entries_, row_map_(row_map_.extent(0)-1));
  }

  void prefix_sum_device(const int& total = -1){
    Kokkos::parallel_scan ("CSR:prefix_sum_device",row_map_.size(), KOKKOS_LAMBDA (
      const int& i, int& upd, const bool& final) {
        // Load old value in case we update it before accumulating
        const float val_i = row_map_(i); 
        if (final) {
          row_map_(i) = upd; // only update array on final pass
        }
        // For exclusive scan, change the update value after
        // updating array, like we do here. For inclusive scan,
        // change the update value before updating array.
        upd += val_i;
    });
    //if (entries_.extent(0) < row_map_(row_map_.extent(0)-1))
    if(total != -1){
      Kokkos::resize(entries_, total);
    }
  }

  // /**
  //  * @brief Compute the prefix sum of the row_map_ 
  //  * and resize the entries based on row_map_
  //  * The CSR is not usable before this operation
  //  */ 
  // void prefix_sum_resize(){
  //   int tmp1 = row_map_(0);
  //   row_map_(0) = 0;
  //   for(int i = 0 ; i < row_map_.extent(0); ++i){ 
  //     int tmp2 = row_map_(i+1);
  //     row_map_(i+1) = tmp1 + row_map_(i);
  //     tmp1 = tmp2;   
  //   }
  //   Kokkos::resize(entries_,row_map_(row_map_.extent(0)-1));
  // }

  /**
   * @brief Set the shape of the ith entry
   * @param i index of the element in row_map_
   * @param shape shape of the ith entry
   * @param size storage for the ith entry
   */ 
  KOKKOS_INLINE_FUNCTION
  void set_shape(const int& i, const int* shape, int size=-1) {
    int total = 1;
    for (int d=0; d!=D; ++d) {
      sizes_(i,d) = shape[d];
      total *= shape[d];
    }
    row_map_(i) = size < 0 ? total : size;
  }


  /**
    * @brief Set the size of the whole csr
  */
  void set_size(const int& row_map_size, const int& entries_size = 0){
    if(row_map_.size() != row_map_size+1)
      Kokkos::resize(row_map_,row_map_size+1);
    if(entries_size != entries_.size())
      Kokkos::resize(entries_,entries_size);
    if(sizes_.extent(0) != row_map_size)
      Kokkos::resize(sizes_,row_map_size,D);   
  }

  KOKKOS_INLINE_FUNCTION
  void set_row_map(const int& i, const int& v){
    row_map_[i] = v; 
  }

  KOKKOS_INLINE_FUNCTION
  void set_entries(const int& i, const double& v){
    entries_[i] = v; 
  }

  KOKKOS_INLINE_FUNCTION
  void set_sizes(const int& i, const int& j, const int& v){
    sizes_(i,j) = v; 
  }

  KOKKOS_INLINE_FUNCTION
  int row_map(const int& i ){
    return row_map_[i]; 
  }

  /** 
   * @brief Return subview corresponding to element i
   * @param i index of the element in row_map_
   */
  KOKKOS_INLINE_FUNCTION
  Kokkos::View<T*,MEMSPACE> at(const int& i) const {
    return Kokkos::subview(
      entries_,Kokkos::make_pair(row_map_(i),row_map_(i+1))); 
  }

  /** 
   * @brief Return size of element idx
   * @param idx Index of the element 
   * @param d Dimension
   */
  KOKKOS_INLINE_FUNCTION 
  int size(int idx, int d) const {
    return sizes_(idx,d);  
  }

  KOKKOS_INLINE_FUNCTION 
  int size(int i) const {
    return row_map_(i+1)-row_map_(i);
  }

  /**
   * @brief Number of elements stored in the CSR 
   */
  KOKKOS_INLINE_FUNCTION
  int size() const{
    return row_map_.extent(0)-1;
  }

  /**
   * @brief Copy constructor 
   */
  void assign(const CSR& csr){
    Kokkos::resize(row_map_,csr.row_map_.extent(0)); 
    Kokkos::deep_copy(row_map_,csr.row_map_);
    Kokkos::resize(entries_,csr.entries_.extent(0)); 
    Kokkos::deep_copy(entries_,csr.entries_);
    Kokkos::resize(sizes_,csr.sizes_.extent(0),csr.sizes_.extent(1)); 
    Kokkos::deep_copy(sizes_,csr.sizes_);
    return *this; 
  }

public: 
  Kokkos::View<int*,MEMSPACE> row_map_; // Indices: number of element +1 
  Kokkos::View<T*,MEMSPACE> entries_; // Values for all entries 
  Kokkos::View<int**,MEMSPACE> sizes_; // Represent the sizes for matrices/tensors

};

// Definitions 
using CSR_Vector = CSR<double,1>; 
using CSR_Matrix = CSR<double,2>; 
using CSR_Tensor = CSR<double,3>; 


} // namespace Amanzi
