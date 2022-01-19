/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon (coonet@ornl.gov)
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

namespace Operators {

class SuperMapLumped {
 public:
  explicit SuperMapLumped(const Teuchos::RCP<const BlockSpace>& maps);

  SuperMapLumped(const SuperMapLumped& other) = delete;
  ~SuperMapLumped();

  // meta-data
  bool HasComponent(const std::string& compname) const;

  // map accessors
  Map_ptr_type getMap() const { return map_; }
  Map_ptr_type getGhostedMap() const { return ghosted_map_; }

  // -- component map accessors
  BlockMap_ptr_type ComponentMap(const std::string& compname) const
  {
    return comp_maps_->ComponentMap(compname, false);
  }

  BlockMap_ptr_type ComponentGhostedMap(const std::string& compname)
  {
    return comp_maps_->ComponentMap(compname, true);
  }

  // index accessors
  template <class DeviceType = DefaultDevice>
  cVectorView_type_<DeviceType, LO>
  Indices(const std::string& compname, int dofnum) const;
  template <class DeviceType = DefaultDevice>
  cVectorView_type_<DeviceType, LO>
  GhostIndices(const std::string& compname, int dofnum) const;

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
  LO Offset(const std::string& compname) const { return offsets_.at(compname); }
  LO GhostedOffset(const std::string& compname) const
  {
    return ghosted_offsets_.at(compname);
  }
  LO NumOwnedElements(const std::string& compname) const
  {
    return counts_.at(compname);
  }
  LO NumUsedElements(const std::string& compname) const
  {
    return counts_.at(compname) + ghosted_counts_.at(compname);
  }
  int getNumVectors(const std::string& compname) const
  {
    return comp_maps_->getNumVectors(compname);
  }

  // iterate over compnames
  using name_iterator = std::vector<std::string>::const_iterator;
  name_iterator begin() const { return comp_maps_->begin(); }
  name_iterator end() const { return comp_maps_->end(); }
  std::size_t size() const { return comp_maps_->size(); }

 public:
  void CreateIndexing_();

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
template <class DeviceType>
cVectorView_type_<DeviceType, LO>
SuperMapLumped::Indices(const std::string& compname, int dofnum) const
{
  return indices_->ViewComponent<DeviceType>(compname, dofnum, false);
}


template <class DeviceType>
cVectorView_type_<DeviceType, LO>
SuperMapLumped::GhostIndices(const std::string& compname, int dofnum) const
{
  return indices_->ViewComponent<DeviceType>(compname, dofnum, true);
}


} // namespace Operators
} // namespace Amanzi

#endif
