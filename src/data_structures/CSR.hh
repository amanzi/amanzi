
/*
  Copyright 2010-202x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Julien Loiseau (jloiseau@lanl.gov)
*/

#pragma once 

#include <tuple>

template<typename T, int D> 
class CSR{
  
static constexpr int dim = D;
using size_type_ = typename std::tuple_element<dim, 
    std::tuple<int, int*, int**>>::type; 

public: 

  CSR(){}

  CSR(const int& row_map_size, const int& entries_size){
    row_map_size_ = row_map_size; 
    entries_size_ = entries_size; 
    Kokkos::resize(row_map_,row_map_size_+1); 
    Kokkos::resize(sizes_,row_map_size_,dim);
    Kokkos::resize(entries_,entries_size_);  
  }

  // Build given the precomputed matrices? 
  //CSR(Kokkos::View<> row_map, Kokkos::View<> entries, ){}

  // Compute prefix sum 
  void prefix_sum(){
    int tmp1 = row_map_(0);
    row_map_(0) = 0;
    for(int i = 0 ; i < row_map_.extent(0); ++i){ 
      int tmp2 = row_map_(i+1);
      row_map_(i+1) = tmp1 + row_map_(i);
      tmp1 = tmp2;   
    }
  }

  void push_back(){

  }

  KOKKOS_INLINE_FUNCTION 
  int size(int idx, int d) const {
    //if(d >= sizes_.extent(1)){
    //  printf("size: %d < %d\n",d,sizes_.extent(1));
    //}
    assert(idx <= sizes_.extent(0)); 
    assert(d <= sizes_.extent(1));
    return sizes_(idx,d);  
  }

  KOKKOS_INLINE_FUNCTION
  int size() const{
    return row_map_.extent(0)-1; 
  }

  Kokkos::View<T*> operator()(const int i){
    Kokkos::View<T*> sv = Kokkos::subview(
      entries_,Kokkos::make_pair(row_map_(i),row_map_(i+1))); 
    return sv; 
  }

  CSR& operator=(const CSR& csr){
    Kokkos::resize(row_map_,csr.row_map_.extent(0)); 
    Kokkos::deep_copy(row_map_,csr.row_map_);
    Kokkos::resize(entries_,csr.entries_.extent(0)); 
    Kokkos::deep_copy(entries_,csr.entries_);
    Kokkos::resize(sizes_,csr.sizes_.extent(0),csr.sizes_.extent(1)); 
    Kokkos::deep_copy(sizes_,csr.sizes_);
    row_map_size_ = csr.row_map_size_; 
    entries_size_ = csr.entries_size_;     
    return *this; 
  }

  Kokkos::View<int*> row_map_;
  Kokkos::View<T*> entries_;
  Kokkos::View<size_type_> sizes_; 
  int row_map_size_, entries_size_; 
};

// Definitions 
using CSR_Matrix = CSR<double,2>; 
using CSR_Vector = CSR<double,1>; 
using CSR_Tensor = CSR<double,2>; 