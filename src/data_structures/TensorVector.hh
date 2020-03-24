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
    data.set_shape(i, {d, rank, tsize}, tsize*tsize);
  }

  void Init() {
    inited = true;
    data.prefix_sum();
  }
  
  KOKKOS_INLINE_FUNCTION
  WhetStone::Tensor operator[](const int& i) {
    assert(inited);
    return std::move(WhetStone::Tensor(data.at(i), data.size(i,0), data.size(i,1), data.size(i,2)));
  }

  KOKKOS_INLINE_FUNCTION
  WhetStone::Tensor at(const int& i) const {
    // FIXME -- not const correct, but to do so needs a const-correct WhetStone::Tensor,
    // e.g. a WhetStone::Tensor that takes a Kokkos::View<const double*> --etc
    return std::move(WhetStone::Tensor(data.at(i), data.size(i,0), data.size(i,1), data.size(i,2)));
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
