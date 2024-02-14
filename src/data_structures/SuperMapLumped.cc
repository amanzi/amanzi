/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  Data Structures

  Takes non-contiguous data structure spaces (CompositeVector, TreeVector)
  and converts them into a single map.
*/

//! <MISSING_ONELINE_DOCSTRING>
#include "AmanziComm.hh"
#include "AmanziMap.hh"
#include "AmanziVector.hh"

#include "BlockSpace.hh"
#include "BlockVector.hh"
#include "TreeVectorSpace.hh"
#include "SuperMapLumped.hh"

namespace Amanzi {


SuperMapLumped::SuperMapLumped(const Teuchos::RCP<const BlockSpace>& space)
  : comp_maps_(space), indices_(std::make_unique<BlockVector<LO>>(space))
{
  createIndexing_();
}

SuperMapLumped::~SuperMapLumped() = default;

//
// NOTE: After this is done, CreateIndices_() may be called.
void
SuperMapLumped::createIndexing_()
{
  // fill the offsets
  {
    int offset = 0;
    for (const auto& compname : *comp_maps_) {
      counts_[compname] = comp_maps_->getComponentMap(compname, false)->getLocalNumElements();
      ghosted_counts_[compname] =
        comp_maps_->getComponentMap(compname, true)->getLocalNumElements() - counts_[compname];
      offsets_[compname] = offset;
      offset += comp_maps_->getNumVectors(compname) * counts_[compname];
    }
    n_local_ = offset;

    // fill the ghosted offsets
    for (const auto& compname : *comp_maps_) {
      ghosted_offsets_[compname] = offset;
      offset += comp_maps_->getNumVectors(compname) * ghosted_counts_[compname];
    }
    n_local_ghosted_ = offset;
  }

  // create indices into supermap and the supermap
  // -- count the total owned points
  GO n_global = indices_->getGlobalLength(false);
  GO n_global_ghosted = indices_->getGlobalLength(true);

  // -- construct the flat map
  map_ = Teuchos::rcp(new Map_type(n_global, n_local_, 0, indices_->getComm()));

  //
  // We need to know the ghost GIDs of the SuperMap, but the connections must
  // be made to maintain the validity of the connection in the various
  // component maps.  Therefore we use the component maps to scatter out their
  // own GIDs in the supermap to the ghosted supermap.
  BlockVector<GO> gids_comp(comp_maps_);
  std::vector<GO> gids_flat(n_local_ghosted_, -1);

  // fill the LIDs and GIDs of owned entities
  for (const auto& compname : *comp_maps_) {
    LO nentities_owned = counts_.at(compname);
    LO nentities = nentities_owned + ghosted_counts_.at(compname);

    auto index_view = indices_->viewComponent<Amanzi::HostSpaceSpecial>(compname, false);
    auto global_index_view = gids_comp.viewComponent<Amanzi::HostSpaceSpecial>(compname, false);

    int offset = offsets_.at(compname);
    int n_dofs = comp_maps_->getNumVectors(compname);
    for (LO i = 0; i != nentities_owned; ++i) {
      for (int j = 0; j != n_dofs; ++j) {
        LO lid = offset + j + i * n_dofs;
        GO gid = map_->getGlobalElement(lid);
        index_view(i, j) = lid;
        global_index_view(i, j) = gid;
        gids_flat[lid] = gid;
      }
    }
  }

  // communicate the GIDs
  gids_comp.scatterMasterToGhosted();

  // fill the LIDs and GIDs of the ghost entities
  for (const auto& compname : *comp_maps_) {
    LO nentities_owned = counts_.at(compname);
    LO nentities = nentities_owned + ghosted_counts_.at(compname);

    auto index_view = indices_->viewComponent<MirrorHost>(compname, true);
    auto global_index_view = gids_comp.viewComponent<MirrorHost>(compname, true);

    int ghosted_offset = ghosted_offsets_.at(compname);
    int n_dofs = comp_maps_->getNumVectors(compname);
    for (LO i = nentities_owned; i != nentities; ++i) {
      for (int j = 0; j != n_dofs; ++j) {
        LO lid = ghosted_offset + j + (i - nentities_owned) * n_dofs;
        index_view(i, j) = lid;
        gids_flat[lid] = global_index_view(i, j);
      }
    }
  }

  // create supermap
  // -- hopefully we found them all!
  if (n_local_ghosted_ > 0) {
    AMANZI_ASSERT(*std::min_element(gids_flat.begin(), gids_flat.end()) >= 0);
  }

  // -- construct
  ghosted_map_ = Teuchos::rcp(
    new Map_type(n_global_ghosted, gids_flat.data(), n_local_ghosted_, 0, indices_->getComm()));


  for (const auto& compname : *comp_maps_) {
    LO nentities_owned = counts_.at(compname);
    LO nentities = nentities_owned + ghosted_counts_.at(compname);
    int n_dofs = comp_maps_->getNumVectors(compname);

    auto index_view = indices_->viewComponent(compname, true);
  }
}


bool
SuperMapLumped::hasComponent(const std::string& key) const
{
  return indices_->hasComponent(key);
}


} // namespace Amanzi
