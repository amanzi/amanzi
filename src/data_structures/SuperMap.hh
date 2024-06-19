/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

//! SuperMap class provides a convenient way of creating and using SuperMapLumped
/*
  Amanzi uses SuperMapLumped in a few different ways that make its natural
  interface not that convenient.  SuperMapLumped also has a few limitations that
  simplify its design and implementaiton a great deal.

  This class is a Helper class in that it wraps a SuperMapLumped, providing a
  related interface that is better designed for users.

  It also enforces and mitigates the design limitations of SuperMapLumped
  itself.
*/

#pragma once

#include <map>
#include "AmanziTypes.hh"
#include "SuperMapLumped.hh"

namespace Amanzi {

class CompositeSpace;
class CompositeVectorSpace;
class TreeVectorSpace;
class TreeVector;

// wrapper class
class SuperMap {
 public:
  using cView_type = typename BlockVector<LO>::cView_type;
  using cHostView_type = typename BlockVector<LO>::cHostView_type;

  SuperMap(const Comm_ptr_type& comm, const std::vector<Teuchos::Ptr<const BlockSpace>>& cvss);

  // map accessors
  // -- global map accessors
  Map_ptr_type getMap() const { return smap_->getMap(); }
  Map_ptr_type getGhostedMap() const { return smap_->getGhostedMap(); }

  // -- component map accessors
  BlockMap_ptr_type getComponentMap(std::size_t block_num, const std::string& compname) const
  {
    return smap_->getComponentMap(block_info_.at(std::make_tuple(block_num, compname, 0)).first);
  }

  BlockMap_ptr_type getComponentGhostedMap(std::size_t block_num, const std::string& compname) const
  {
    return smap_->getComponentGhostedMap(
      block_info_.at(std::make_tuple(block_num, compname, 0)).first);
  }

  // check if accessor is valid
  bool
  hasComponent(std::size_t block_num, const std::string& compname, std::size_t dof_num = 0) const
  {
    return block_info_.count(std::make_tuple(block_num, compname, dof_num)) != 0;
  }

  // index accessors
  template <MemSpace_kind MEM = MemSpace_kind::DEVICE>
  auto
  viewIndices(std::size_t block_num, const std::string& compname, std::size_t dof_num) const
  {
    auto bi = block_info_.find(std::make_tuple(block_num, compname, dof_num));
    if (bi == block_info_.end()) {
      Errors::Message msg;
      msg << "SuperMap does not have block component <" << (int)block_num << "," << compname << ","
          << (int)dof_num;
      Exceptions::amanzi_throw(msg);
    }
    return smap_->viewIndices<MEM>(bi->second.first, bi->second.second);
  }

  template <MemSpace_kind MEM = MemSpace_kind::DEVICE>
  auto
  viewGhostIndices(std::size_t block_num, const std::string& compname, std::size_t dof_num) const
  {
    auto bi = block_info_.find(std::make_tuple(block_num, compname, dof_num));
    if (bi == block_info_.end()) {
      Errors::Message msg;
      msg << "SuperMap does not have block component <" << (int)block_num << "," << compname << ","
          << (int)dof_num;
      Exceptions::amanzi_throw(msg);
    }
    return smap_->viewGhostIndices<MEM>(bi->second.first, bi->second.second);
  }

  // block indices.  This is an array of integers, length
  // Map().getLocalLength(), where each dof and component have a unique integer
  // value.  The returned int is the number of unique values, equal to
  // sum(NumDofs(comp) for comp in components), in this array.
  // std::pair<int, Teuchos::RCP<std::vector<int> > > BlockIndices() const {
  //   return smap_->BlockIndices();
  // }

 protected:
  std::unique_ptr<SuperMapLumped> smap_;
  std::map<std::tuple<int, std::string, int>, std::pair<std::string, int>> block_info_;
};


// Nonmember contructors/factories
Teuchos::RCP<SuperMap>
createSuperMap(const Teuchos::Ptr<const BlockSpace>& cv);
Teuchos::RCP<SuperMap>
createSuperMap(const CompositeVectorSpace& cv);
Teuchos::RCP<SuperMap>
createSuperMap(const TreeVectorSpace& cv);

// Copy in/out
int
copyToSuperVector(const SuperMap& map,
                  const BlockVector<Vector_type::scalar_type>& bv,
                  Vector_type& sv,
                  int block_num = 0);
int
copyFromSuperVector(const SuperMap& map,
                    const Vector_type& sv,
                    BlockVector<Vector_type::scalar_type>& bv,
                    int block_num = 0);
int
addFromSuperVector(const SuperMap& map,
                   const Vector_type& sv,
                   BlockVector<Vector_type::scalar_type>& bv,
                   int block_num = 0);

// Nonmember TreeVector to/from Super-vector
// -- simple schema version
int
copyToSuperVector(const SuperMap& map, const TreeVector& tv, Vector_type& sv);
int
copyFromSuperVector(const SuperMap& map, const Vector_type& sv, TreeVector& tv);
int
addFromSuperVector(const SuperMap& map, const Vector_type& sv, TreeVector& tv);


} // namespace Amanzi
