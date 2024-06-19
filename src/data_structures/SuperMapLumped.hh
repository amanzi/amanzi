/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

//! SuperMapLumped class for blocking same-type maps into a single, common map.
/*!
  Takes non-contiguous data structure spaces (CompositeVector, TreeVector)
  and converts them into a single map.

  DESIGN FLAW: Currently this assumes that component names are unique, and if
  two components share the same name, they share the same map.  This is
  obviously wrong when multple meshes are involved -- for instance a TV of
  surface + subsurface, both with "cell" components, would break miserably.

  That said, this is a simple and well-posed class.  It stays, but users
  should prefer to use the wrapper class, SuperMap, which deals with a
  few shortcomings of this class.
*/

#ifndef AMANZI_OPERATORS_SUPER_MAP_HH_
#define AMANZI_OPERATORS_SUPER_MAP_HH_

#include "Teuchos_RCP.hpp"

#include "AmanziTypes.hh"
#include "dbc.hh"
#include "Mesh.hh"
#include "CompositeSpace.hh"
#include "BlockVector.hh"

namespace Amanzi {

// template<typename Scalar> class BlockVector;

class SuperMapLumped {
 public:
  using cView_type = typename BlockVector<LO>::cView_type;
  using cHostView_type = typename BlockVector<LO>::cHostView_type;

  explicit SuperMapLumped(const Teuchos::RCP<const BlockSpace>& maps);

  SuperMapLumped(const SuperMapLumped& other) = delete;
  ~SuperMapLumped();

  // meta-data
  bool hasComponent(const std::string& compname) const;

  // map accessors
  Map_ptr_type getMap() const { return map_; }
  Map_ptr_type getGhostedMap() const { return ghosted_map_; }

  // -- component map accessors
  BlockMap_ptr_type getComponentMap(const std::string& compname) const
  {
    return comp_maps_->getComponentMap(compname, false);
  }

  BlockMap_ptr_type getComponentGhostedMap(const std::string& compname)
  {
    return comp_maps_->getComponentMap(compname, true);
  }

  // index accessors
  template <MemSpace_kind MEM = MemSpace_kind::DEVICE>
  auto
  viewIndices(const std::string& compname, int dofnum) const;

  template <MemSpace_kind MEM = MemSpace_kind::DEVICE>
  auto
  viewGhostIndices(const std::string& compname, int dofnum) const;

  // block indices.  This is an array of integers, length
  // Map().getLocalLength(), where each dof and component have a unique integer
  // value.  The returned int is the number of unique values, equal to
  // sum(getNumVectors(comp) for comp in components), in this array.
  //  std::pair<int, Teuchos::RCP<std::vector<int> > > BlockIndices() const;

#ifdef SUPERMAP_TESTING
 public:
#else
 protected:
#endif

  // meta-data accessors
  LO getOffset(const std::string& compname) const { return offsets_.at(compname); }
  LO getGhostedOffset(const std::string& compname) const { return ghosted_offsets_.at(compname); }
  LO getNumOwnedElements(const std::string& compname) const { return counts_.at(compname); }
  LO getNumUsedElements(const std::string& compname) const
  {
    return counts_.at(compname) + ghosted_counts_.at(compname);
  }
  int getNumVectors(const std::string& compname) const
  {
    return comp_maps_->getNumVectors(compname);
  }

 protected:
  // iterate over compnames
  using name_iterator = std::vector<std::string>::const_iterator;
  name_iterator begin() const { return comp_maps_->begin(); }
  name_iterator end() const { return comp_maps_->end(); }
  std::size_t size() const { return comp_maps_->size(); }

 public:
  void createIndexing_();

 protected:
  std::map<std::string, LO> offsets_;
  std::map<std::string, LO> counts_;
  std::map<std::string, LO> ghosted_offsets_;
  std::map<std::string, LO> ghosted_counts_;

  LO n_local_;
  LO n_local_ghosted_;

  std::unique_ptr<BlockVector<LO>> indices_;
  Teuchos::RCP<const BlockSpace> comp_maps_;

  Map_ptr_type map_;         // the supermap
  Map_ptr_type ghosted_map_; // the ghosted supermap
};


Teuchos::RCP<SuperMapLumped>
createSuperMapLumped(const BlockSpace& cv);

// implementation of templated member functions
template <MemSpace_kind MEM>
auto
SuperMapLumped::viewIndices(const std::string& compname, int dofnum) const
{
  // force returning const view
  const BlockVector<LO>& inds = *indices_;
  return inds.viewComponent<MEM>(compname, dofnum, false);
}


template <MemSpace_kind MEM>
auto
SuperMapLumped::viewGhostIndices(const std::string& compname, int dofnum) const
{
  // force returning const view
  const BlockVector<LO>& inds = *indices_;
  return inds.viewComponent<MEM>(compname, dofnum, true);
}


} // namespace Amanzi

#endif
