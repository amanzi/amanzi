/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//! Helper data structure and factory for using WhetStone Tensors

#ifndef AMANZI_TENSOR_VECTOR_HH_
#define AMANZI_TENSOR_VECTOR_HH_

#include "AmanziTypes.hh"
#include "AmanziVector.hh"
#include "MeshDefs.hh"
#include "Tensor.hh"
#include "CSR.hh"
#include "CompositeVectorSpace.hh"
#include "Mesh.hh"

namespace Amanzi {

//
// A simple data structure that keeps a vector of WhetStone Tensors and a
// CompositeVectorSpace to describe the layout/ghosting/entities those tensors
// are associated with.
// -----------------------------------------------------------------------------
struct TensorVector {
  TensorVector() {}
  
  TensorVector(const CompositeVectorSpace& map_, bool ghosted_ = false)
      : map(map_), ghosted(ghosted_), inited(false)
  {
    prealloc_();
  }


  TensorVector(const CompositeVectorSpace& map_, int dim_, int rank_,
               bool ghosted_ = false)
      : map(map_), ghosted(ghosted_), inited(false)
  {
    prealloc_();
    for (int i=0; i!=size(); ++i) set_shape(i, dim_, rank_);
    Init();
  }

  KOKKOS_INLINE_FUNCTION size_t size() const { return data.size(); }

  void set_shape(const int& i, const int& d, const int& rank) {
    int tsize = WhetStone::WHETSTONE_TENSOR_SIZE[d-1][rank-1];
    int loc[3] = {d, rank, tsize};
    data.set_shape_host(i, loc, tsize*tsize);
  }

  /**
    * Init the CSR based on a Host vector 
  */
  void Init(const std::vector<WhetStone::Tensor<Kokkos::HostSpace>>& ht){
    // Extract max size 
    const int row_map_size = ht.size(); 
    // Set CSR 
    data.set_size(row_map_size); 
    // Fill row map  and size
    for(int i = 0 ; i < ht.size(); ++i){
      data.set_row_map_host(i,ht[i].mem());
      data.set_sizes_host(i,0,ht[i].dimension());
      data.set_sizes_host(i,1,ht[i].rank()); 
      data.set_sizes_host(i,2,ht[i].size());
    }
    data.prefix_sum(); 
    // Fill entries 
    for(int i = 0 ; i < ht.size(); ++i){
      int idx = data.row_map_host(i);
      for(int j = 0 ; j < ht[i].data().extent(0); ++j){
        data.set_entries_host(idx+j,ht[i].data()(j)); 
      }
    }
    data.update_entries_device(); 
    data.update_sizes_device(); 
    data.update_row_map_device(); 
  }

  void Init() {
    inited = true;
    data.prefix_sum();
  }
  
  // The operator[] return the value on device 
  KOKKOS_INLINE_FUNCTION
  WhetStone::Tensor<DeviceOnlyMemorySpace> operator[](const int& i) {
    assert(inited);
    return std::move(WhetStone::Tensor<DeviceOnlyMemorySpace>(data.at(i), data.size(i,0), data.size(i,1), data.size(i,2)));
  }

  KOKKOS_INLINE_FUNCTION
  WhetStone::Tensor<Kokkos::HostSpace> at_host(const int& i) const {
    // FIXME -- not const correct, but to do so needs a const-correct WhetStone::Tensor,
    // e.g. a WhetStone::Tensor that takes a Kokkos::View<const double*> --etc
    return std::move(WhetStone::Tensor<Kokkos::HostSpace>(
        data.at_host(i), data.size_host(i,0), data.size_host(i,1), data.size_host(i,2)));
  }
  

  KOKKOS_INLINE_FUNCTION
  WhetStone::Tensor<DeviceOnlyMemorySpace> at(const int& i) const {
    // FIXME -- not const correct, but to do so needs a const-correct WhetStone::Tensor,
    // e.g. a WhetStone::Tensor that takes a Kokkos::View<const double*> --etc
    return std::move(WhetStone::Tensor<DeviceOnlyMemorySpace>(data.at(i), data.size(i,0), data.size(i,1), data.size(i,2)));
  }
  
  CSR_Tensor data;
  CompositeVectorSpace map;
  bool ghosted;
  bool inited;

 private:

  void prealloc_() {
    // how many entities?
    int count = 0;
    for (const auto& comp : map)
      count += map.getMap(comp, ghosted)->getNodeNumElements();

    data = std::move(CSR_Tensor(count));
  }

};


//
// A factory for putting these into state.
// -----------------------------------------------------------------------------
class TensorVector_Factory {
 public:
  TensorVector_Factory() : d_(0), rank_(0), ghosted_(false) {}

  int dimension() const { return d_; }
  int rank() const { return rank_; }
  const CompositeVectorSpace& map() const { return map_; }

  void set_rank(int rank) { rank_ = rank; }
  void set_map(CompositeVectorSpace map)
  {
    map_ = std::move(map);
    d_ = map_.Mesh()->space_dimension();
  }
  void set_ghosted(bool ghosted) { ghosted_ = ghosted; }

  Teuchos::RCP<TensorVector> Create() const
  {
    return Teuchos::rcp(new TensorVector(map_, d_, rank_, ghosted_));
  }

 private:
  CompositeVectorSpace map_;
  int d_, rank_;
  bool ghosted_;
};

} // namespace Amanzi


#endif
