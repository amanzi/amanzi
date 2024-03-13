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

#include "Epetra_IntVector.h"

#include "CompositeVectorSpace.hh"
#include "TreeVectorSpace.hh"
#include "TreeVector_Utils.hh"
#include "SuperMapLumped.hh"

namespace Amanzi {
namespace Operators {

SuperMapLumped::SuperMapLumped(const Comm_ptr_type& comm,
                               const std::vector<std::string>& compnames,
                               const std::vector<int>& dofnums,
                               const std::vector<Teuchos::RCP<const Epetra_BlockMap>>& maps,
                               const std::vector<Teuchos::RCP<const Epetra_BlockMap>>& ghosted_maps)
  : compnames_(compnames)
{
  AMANZI_ASSERT(compnames.size() == dofnums.size());
  AMANZI_ASSERT(compnames.size() == maps.size());
  AMANZI_ASSERT(compnames.size() == ghosted_maps.size());

  for (int i = 0; i != compnames.size(); ++i) {
    comp_maps_[compnames[i]] = maps[i];
    comp_ghosted_maps_[compnames[i]] = ghosted_maps[i];
    num_dofs_[compnames[i]] = dofnums[i];
  }

  CreateIndexing_();
  CreateMap_(comm);
}

//
// NOTE: After this is done, CreateIndices_() may be called.
void
SuperMapLumped::CreateIndexing_()
{
  // fill the offsets
  int offset = 0;
  for (const auto& compname : compnames_) {
    counts_[compname] = comp_maps_[compname]->NumMyPoints();
    ghosted_counts_[compname] = comp_ghosted_maps_[compname]->NumMyPoints() - counts_[compname];
    offsets_[compname] = offset;
    offset += num_dofs_[compname] * counts_[compname];
  }
  n_local_ = offset;

  // fill the ghosted offsets
  for (const auto& compname : compnames_) {
    ghosted_offsets_[compname] = offset;
    offset += num_dofs_[compname] * ghosted_counts_[compname];
  }
  n_local_ghosted_ = offset;
}


// create the owned and ghosted SuperMapLumpeds
//
// NOTE: previous versions of this respected the GID ordering of the
// incoming map.  There is no real reason to do this -- we're reordering
// everything on purpose, so we can do whatever we want with the GIDs.  But
// the previous version took advantage of the fact that we know, on process,
// what our GID is, and therefore can set up ghost IDs for the SuperMapLumped
// without communication.  The BlockMap version doesn't seem to have a way
// to get the Point location of a ghosted entity.  Since we seemingly can't
// do this without communication, we might as well make it easier on
// ourselves, by just making Epetra assign SuperMapLumped's GIDs as it wants.
void
SuperMapLumped::CreateMap_(const Comm_ptr_type& comm)
{
  // -- count the total owned points
  int n_global = 0;
  int n_global_ghosted = 0;
  for (const auto& compname : compnames_) {
    n_global += num_dofs_[compname] * comp_maps_[compname]->NumGlobalPoints();
    n_global_ghosted += num_dofs_[compname] * comp_ghosted_maps_[compname]->NumGlobalPoints();
  }

  // -- construct
  map_ = Teuchos::rcp(new Epetra_Map(n_global, n_local_, 0, *comm));

  // create the ghosted map via communication using the provided BlockMaps
  //
  // We need to know the ghost GIDs of the SuperMapLumped, but the connections must
  // be made to maintain the validity of the connection in the various
  // component maps.  Therefore we use the component maps to scatter out their
  // own GIDs in the supermap to the ghosted supermap.
  std::vector<int> ghosted_gids(n_local_ghosted_, -1);
  for (const auto& compname : compnames_) {
    int n_dofs = num_dofs_[compname];

    const auto& comp_map = *comp_maps_[compname];
    const auto& comp_ghosted_map = *comp_ghosted_maps_[compname];

    // Put the supermap GID into an IntVector on the (owned) component map
    //
    // this is a bit ugly, but it is guaranteed safe due to ordering of the
    // construction.  A better way to do this would have CreateIndices_ make
    // IntVectors directly.
    //
    // Also, newer versions of Epetra have Epetra_IntMultiVector, which would
    // be useful here.
    Epetra_Import importer(comp_ghosted_map, comp_map);
    Epetra_IntVector owned_gids_c(comp_map);
    Epetra_IntVector ghosted_gids_c(comp_ghosted_map);

    for (int j = 0; j != n_dofs; ++j) {
      const auto& owned_lids = CreateIndices_(compname, j, false);
      for (int k = 0; k != owned_gids_c.MyLength(); ++k) {
        owned_gids_c[k] = map_->GID(owned_lids[k]);
      }

      // Scatter supermap GIDs to a ghosted component vector
      ghosted_gids_c.Import(owned_gids_c, importer, Insert);

      // gids for the ghosted supermap
      const auto& ghosted_inds = CreateIndices_(compname, j, true);
      for (int k = 0; k != ghosted_gids_c.MyLength(); ++k) {
        ghosted_gids[ghosted_inds[k]] = ghosted_gids_c[k];
      }
    }
  }

  // hopefully we found them all!
  if (n_local_ghosted_ > 0) {
    AMANZI_ASSERT(*std::min_element(ghosted_gids.begin(), ghosted_gids.end()) >= 0);
  }

  // -- construct
  ghosted_map_ =
    Teuchos::rcp(new Epetra_Map(n_global_ghosted, n_local_ghosted_, ghosted_gids.data(), 0, *comm));
}


bool
SuperMapLumped::HasComponent(const std::string& key) const
{
  std::map<std::string, int>::const_iterator lb = offsets_.lower_bound(key);
  if (lb != offsets_.end() && !(offsets_.key_comp()(key, lb->first))) {
    return true;
  } else {
    return false;
  }
}


const std::vector<int>&
SuperMapLumped::Indices(const std::string& compname, int dofnum) const
{
  if (indices_.count(compname)) {
    if (indices_[compname].count(dofnum)) { return indices_[compname][dofnum]; }
  }
  return CreateIndices_(compname, dofnum, false);
}


const std::vector<int>&
SuperMapLumped::GhostIndices(const std::string& compname, int dofnum) const
{
  if (ghosted_indices_.count(compname)) {
    if (ghosted_indices_[compname].count(dofnum)) { return ghosted_indices_[compname][dofnum]; }
  }
  return CreateIndices_(compname, dofnum, true);
}


std::pair<int, Teuchos::RCP<std::vector<int>>>
SuperMapLumped::BlockIndices() const
{
  auto block_indices = Teuchos::rcp(new std::vector<int>(Map()->NumMyElements()));
  int block_id = 0;
  for (const auto& comp : *this) {
    int ndofs = NumDofs(comp);
    for (int d = 0; d != ndofs; ++d) {
      const auto& inds = Indices(comp, d);
      for (int i = 0; i != inds.size(); ++i) (*block_indices)[inds[i]] = block_id;
      block_id++;
    }
  }
  return std::make_pair(block_id, block_indices);
}


const std::vector<int>&
SuperMapLumped::CreateIndices_(const std::string& compname, int dofnum, bool ghosted) const
{
  if (ghosted) {
    if (ghosted_indices_.count(compname) == 0) { ghosted_indices_[compname]; }

    // create the vector
    int nentities_owned = counts_.at(compname);
    int nentities = nentities_owned + ghosted_counts_.at(compname);

    std::vector<int> indices(nentities, -1);
    int offset = offsets_.at(compname);
    int num_dof = num_dofs_.at(compname);
    for (int i = 0; i != nentities_owned; ++i) { indices[i] = offset + dofnum + i * num_dof; }

    int ghosted_offset = ghosted_offsets_.at(compname);
    for (int i = nentities_owned; i != nentities; ++i) {
      indices[i] = ghosted_offset + dofnum + (i - nentities_owned) * num_dof;
    }

    // move-assign, indices is no longer valid
    ghosted_indices_[compname][dofnum] = std::move(indices);
    return ghosted_indices_[compname][dofnum];

  } else {
    if (indices_.count(compname) == 0) { indices_[compname]; }

    // create the vector
    int nentities = counts_.at(compname);

    std::vector<int> indices(nentities, -1);
    int offset = offsets_.at(compname);
    int num_dof = num_dofs_.at(compname);
    for (int i = 0; i != nentities; ++i) { indices[i] = offset + dofnum + i * num_dof; }

    // move-assign, indices is no longer valid
    indices_[compname][dofnum] = std::move(indices);

    return indices_[compname][dofnum];
  }
}


// nonmember function
Teuchos::RCP<SuperMapLumped>
createSuperMapLumped(const CompositeVectorSpace& cv)
{
  std::vector<std::string> names;
  std::vector<int> dofnums;
  std::vector<Teuchos::RCP<const Epetra_BlockMap>> maps, ghostmaps;

  for (const auto& name : cv) {
    names.push_back(name);
    dofnums.push_back(cv.NumVectors(name));
    maps.push_back(cv.Map(name, false));
    ghostmaps.push_back(cv.Map(name, true));
  }
  return Teuchos::rcp(new SuperMapLumped(cv.Comm(), names, dofnums, maps, ghostmaps));
}


} // namespace Operators
} // namespace Amanzi
