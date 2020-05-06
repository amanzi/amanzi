/*
  Copyright 2010-202x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Julien Loiseau (jloiseau@lanl.gov)
*/

#pragma once 

#include "Kokkos_Core.hpp"
#include "Kokkos_DualView.hpp"
#include "AmanziTypes.hh"

#include "DenseMatrix.hh"
#include "DenseVector.hh"

namespace Amanzi {

/**
 * @brief CSR data structure 
 * @param T Type used for the entries_ 
 * @param D Dimension = number of dimensions allocated in sizes
 * This is used to reconstruct the Matrix/Tensor on the device
 */
template<typename T, int D, class MEMSPACE = DefaultHostMemorySpace> 
class CSR{

static constexpr int dim = D;

public: 
  using memory_space = MEMSPACE;
  
  CSR() {}

  CSR(const int& row_map_size, const int& entries_size)
  {
    row_map_.realloc(row_map_size+1); 
    sizes_.realloc(row_map_size,dim);
    entries_.realloc(entries_size);  
  }

  CSR(const int& row_map_size)
  {
    row_map_.realloc(row_map_size+1); 
    sizes_.realloc(row_map_size,dim);
  }

  template<typename T2>
  CSR(const CSR<T2,D,MEMSPACE>& other) :
      row_map_(other.row_map_),
      sizes_(other.sizes_)
  {
    entries_.realloc(other.entries_.extent(0));
  }

  // Update functions 
  void update_row_map_device(){
    row_map_.modify_host(); row_map_.sync_device(); 
  }
  void update_sizes_device(){
    sizes_.modify_host(); sizes_.sync_device(); 
  }
  void update_entries_device(){
    entries_.modify_host(); entries_.sync_device(); 
  }

  void update_row_map_host(){
    row_map_.modify_device(); row_map_.sync_host(); 
  }
  void update_sizes_host(){
    sizes_.modify_device(); sizes_.sync_host(); 
  }
  void update_entries_host(){
    entries_.modify_device(); entries_.sync_host(); 
  }

  /**
   * @brief Compute the prefix sum of the row_map_ 
   * The CSR is not usable before this operation
   */ 
  void prefix_sum(){
    int tmp1 = row_map_.view_host()(0);
    row_map_.view_host()(0) = 0;
    for(int i = 0 ; i < row_map_.view_host().extent(0); ++i){ 
      int tmp2 = row_map_.view_host()(i+1);
      row_map_.view_host()(i+1) = tmp1 + row_map_.view_host()(i);
      tmp1 = tmp2;   
    }
    if (entries_.view_host().extent(0) < 
          row_map_.view_host()(row_map_.view_host().extent(0)-1))
      entries_.realloc(row_map_.view_host()(row_map_.extent(0)-1));
    // Transfer to device 
    update_row_map_device(); 
  }

  void prefix_sum_device(const int& total = -1){
    Kokkos::parallel_scan (
      "CSR:prefix_sum_device",row_map_.view_device().size(), KOKKOS_LAMBDA (
      const int& i, int& upd, const bool& final) {
        const float val_i = row_map_.view_device()(i); 
        if (final) {
          row_map_.view_device()(i) = upd;
        }
        upd += val_i;
    });
    if(total != -1){
      entries_.realloc(total);
    }
  }



  // Host side functions ------------------------------------------------------

  /**
   * @brief Set the shape of the ith entry
   * @param i index of the element in row_map_
   * @param shape shape of the ith entry
   * @param size storage for the ith entry
   */ 
  void set_shape_host(const int& i, const int* shape, int size=-1) {
    int total = 1;
    for (int d=0; d!=D; ++d) {
      sizes_.view_host()(i,d) = shape[d];
      total *= shape[d];
    }
    row_map_.view_host()(i) = size < 0 ? total : size;
  }

  /**
    * @brief Set the size of the whole csr
  */
  void set_size(const int& row_map_size, const int& entries_size = 0){
    if(row_map_.view_host().size() != row_map_size+1)
      row_map_.realloc(row_map_size+1);
    if(entries_.view_host().size() != entries_size)
      entries_.realloc(entries_size);
    if(sizes_.view_host().extent(0) != row_map_size)
      sizes_.realloc(row_map_size,D);   
  }

  int size_host() const{
    return row_map_.view_host().extent(0)-1;
  }
  int size_host(int idx, int d) const {
    return sizes_.view_host()(idx,d);  
  }
  int size_host(int i) const {
    return row_map_.view_host()(i+1)-row_map_.view_host()(i);
  }
  /** 
   * @brief Return subview corresponding to element i
   * @param i index of the element in row_map_
   */
  Kokkos::View<T*,Kokkos::HostSpace> at_host(const int& i) const {
    return Kokkos::subview(
      entries_.view_host(),
        Kokkos::make_pair(
          row_map_.view_host()(i),row_map_.view_host()(i+1))); 
  }
  void set_row_map_host(const int& i, const int& v){
    row_map_.view_host()[i] = v; 
  }
  void set_entries_host(const int& i, const double& v){
    entries_.view_host()[i] = v; 
  }
  void set_sizes_host(const int& i, const int& j, const int& v){
    sizes_.view_host()(i,j) = v; 
  }
  int row_map_host(const int& i ){
    return row_map_.view_host()[i]; 
  }


  // Device side functions = default name -------------------------------------


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
      sizes_.view_device()(i,d) = shape[d];
      total *= shape[d];
    }
    row_map_.view_device()(i) = size < 0 ? total : size;
  }

  KOKKOS_INLINE_FUNCTION
  int size() const{
    return row_map_.view_device().extent(0)-1;
  }
  KOKKOS_INLINE_FUNCTION 
  int size(int idx, int d) const {
    return sizes_.view_device()(idx,d);  
  }
  KOKKOS_INLINE_FUNCTION 
  int size(int i) const {
    return row_map_.view_device()(i+1)-row_map_.view_device()(i);
  }
  /** 
   * @brief Return subview corresponding to element i
   * @param i index of the element in row_map_
   */
  KOKKOS_INLINE_FUNCTION
  Kokkos::View<T*,MEMSPACE> at(const int& i) const {
    return Kokkos::subview(
      entries_.view_device(),
        Kokkos::make_pair(
          row_map_.view_device()(i),row_map_.view_device()(i+1))); 
  }

  KOKKOS_INLINE_FUNCTION
  void set_row_map(const int& i, const int& v){
    row_map_.view_device()[i] = v; 
  }

  KOKKOS_INLINE_FUNCTION
  void set_entries(const int& i, const double& v){
    entries_.view_device()[i] = v; 
  }

  KOKKOS_INLINE_FUNCTION
  void set_sizes(const int& i, const int& j, const int& v){
    sizes_.view_device()(i,j) = v; 
  }

  KOKKOS_INLINE_FUNCTION
  int row_map(const int& i ){
    return row_map_.view_device()[i]; 
  }


  /**
   * @brief Copy constructor 
   */
  void assign(const CSR& csr){
    row_map_.realloc(csr.row_map_.view_host().extent(0)); 
    Kokkos::deep_copy(row_map_,csr.row_map_);
    entries_.realloc(csr.entries_.view_host().extent(0)); 
    Kokkos::deep_copy(entries_,csr.entries_);
    sizes_.realloc(csr.sizes_.view_host().extent(0),
                  csr.sizes_.view_host().extent(1)); 
    Kokkos::deep_copy(sizes_,csr.sizes_);
    return *this; 
  }

public: 
  Kokkos::DualView<int*,MEMSPACE> row_map_; // Indices: number of element +1 
  Kokkos::DualView<T*,MEMSPACE> entries_; // Values for all entries 
  Kokkos::DualView<int**,MEMSPACE> sizes_; // Represent the sizes for matrices/tensors

};


template<class MEMSPACE=DefaultHostMemorySpace> 
class CSR_DenseVector : public CSR<double,1,MEMSPACE> {
 public:
  KOKKOS_INLINE_FUNCTION
  WhetStone::DenseVector<MEMSPACE>
  at(const int& i) const {
    return WhetStone::DenseVector<MEMSPACE>(CSR<double,1,MEMSPACE>::at(i), this->size(i,0));
  }
};

template<class MEMSPACE=DefaultHostMemorySpace> 
class CSR_DenseMatrix : public CSR<double,2,MEMSPACE> {
 public:
  KOKKOS_INLINE_FUNCTION
  WhetStone::DenseMatrix<MEMSPACE>
  at(const int& i) const {
    return WhetStone::DenseVector<MEMSPACE>(CSR<double,1,MEMSPACE>::at(i), this->size(i,0), this->size(i,1));
  }
};

using CSR_Vector = CSR<double,1,DeviceOnlyMemorySpace>;
using CSR_IntVector = CSR<int,1,DeviceOnlyMemorySpace>;
using CSR_Matrix = CSR<double,2,DeviceOnlyMemorySpace>;
using CSR_Tensor = CSR<double,3,DeviceOnlyMemorySpace>; 




} // namespace Amanzi
