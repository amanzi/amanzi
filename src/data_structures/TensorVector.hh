//! Helper factory for storing WhetStone Tensors in State

/*
  WhetStone, version X.Y
  Release name: dev

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)

  A simple wrapper for creating a vector of tensors.
*/

#ifndef AMANZI_TENSOR_VECTOR_HH_
#define AMANZI_TENSOR_VECTOR_HH_

#include <vector>
#include <algorithm>

#include "MeshDefs.hh"
#include "Tensor.hh"
#include "CompositeVectorSpace.hh"

namespace Amanzi {

//
// A simple data structure that keeps a vector of WhetStone Tensors and a
// CompositeVectorSpace to describe the layout/ghosting/entities those tensors
// are associated with.
// -----------------------------------------------------------------------------
struct TensorVector {
  TensorVector(CompositeVectorSpace map_, bool ghosted_=false) :
      map(std::move(map_)),
      ghosted(ghosted_) {
    data.resize(size_());
  }

  TensorVector(CompositeVectorSpace map_, int dim_, int rank_, bool ghosted_=false) :
      map(std::move(map_)),
      ghosted(ghosted_) {
    data.resize(size_(), WhetStone::Tensor(dim_, rank_));
  }

  std::vector<WhetStone::Tensor>::iterator begin() { return data.begin(); }
  std::vector<WhetStone::Tensor>::const_iterator begin() const { return data.begin(); }
  std::vector<WhetStone::Tensor>::iterator end() { return data.end(); }
  std::vector<WhetStone::Tensor>::const_iterator end() const { return data.end(); }
  std::size_t size() const { return data.size(); }

  const WhetStone::Tensor& operator[](std::size_t i) const { return data[i]; }
  WhetStone::Tensor& operator[](std::size_t i) { return data[i]; }

  TensorVector& operator=(const TensorVector& other) = default;

  std::vector<WhetStone::Tensor> data;
  CompositeVectorSpace map;
  bool ghosted;

 private:
  int size_() {
    int count = 0;
    for (auto& name : map) {
      count += map.Mesh()->num_entities(AmanziMesh::entity_kind(name), 
          ghosted ? AmanziMesh::Parallel_type::ALL : AmanziMesh::Parallel_type::OWNED);
    }
    return count;
  }
};


//
// A factory for putting these into state.
// -----------------------------------------------------------------------------
class TensorVector_Factory {
 public:
  TensorVector_Factory() :
      d_(0),
      rank_(0),
      ghosted_(false)      
  {};

  int dimension() const { return d_; }
  const CompositeVectorSpace& map() const { return map_; }

  int get_rank() const { return rank_; }
  void set_rank(int rank) { rank_ = rank; }

  void set_map(CompositeVectorSpace map) {
    map_ = std::move(map);
    d_ = map_.Mesh()->space_dimension();
  }
  void set_ghosted(bool ghosted) {
    ghosted_ = ghosted;
  }
  
  Teuchos::RCP<TensorVector> Create() const {
    return Teuchos::rcp(new TensorVector(map_, d_, rank_, ghosted_));
  }

 private:
  CompositeVectorSpace map_;
  int d_, rank_;
  bool ghosted_;
};

}  // namespace Amanzi

#endif

