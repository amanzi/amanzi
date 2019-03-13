/*
  Data Structures

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)

  Takes non-contiguous data structure spaces (CompositeVector, TreeVector) 
  and converts them into a single map.
*/


#include "CompositeVectorSpace.hh"
#include "TreeVectorSpace.hh"
#include "TreeVector_Utils.hh"
#include "SuperMap.hh"

namespace Amanzi {
namespace Operators {

SuperMap::SuperMap(const Comm_ptr_type& comm,
                   const std::vector<std::string>& compnames,
                   const std::vector<int>& dofnums,
                   const std::vector<Teuchos::RCP<const Epetra_BlockMap> >& maps,
                   const std::vector<Teuchos::RCP<const Epetra_BlockMap> >& ghosted_maps) :
    compnames_(compnames)
{
  AMANZI_ASSERT(compnames.size() == dofnums.size());
  AMANZI_ASSERT(compnames.size() == maps.size());
  AMANZI_ASSERT(compnames.size() == ghosted_maps.size());

  for (int i=0; i!=compnames.size(); ++i) {
    comp_maps_[compnames[i]] = maps[i];
    comp_ghosted_maps_[compnames[i]] = ghosted_maps[i];
  }
  
  int offset = 0;

  // fill the offsets
  for (int i=0; i!=compnames.size(); ++i) {
    std::string compname = compnames[i];
    num_dofs_[compname] = dofnums[i];
    counts_[compname] = maps[i]->NumMyPoints();
    ghosted_counts_[compname] = ghosted_maps[i]->NumMyPoints() - counts_[compname];
    offsets_[compname] = offset;

    offset += dofnums[i]*counts_[compname];
  }
  int n_local = offset;

  // fill the ghosted offsets
  for (int i=0; i!=compnames.size(); ++i) {
    std::string compname = compnames[i];
    ghosted_offsets_[compname] = offset;
    offset += dofnums[i]*ghosted_counts_[compname];
  }
  int n_local_ghosted = offset;

  // populate the map GIDs
  int global_offset = 0;
  int n(0);
  std::vector<int> gids(offset, -1);
  for (int i=0; i!=compnames.size(); ++i) {
    for (int j=0; j!=counts_[compnames[i]]; ++j) {
      int base = global_offset + ghosted_maps[i]->GID(j)*dofnums[i];
      for (int dof=0; dof!=dofnums[i]; ++dof) gids[n++] = base+dof;
    }
    global_offset += maps[i]->NumGlobalPoints() * dofnums[i];
  }
  int n_global = global_offset;

  int n_global_ghosted = 0;
  global_offset = 0;
  for (int i=0; i!=compnames.size(); ++i) {
    for (int j=counts_[compnames[i]]; j!=(counts_[compnames[i]]+ghosted_counts_[compnames[i]]); ++j) {
      int base = global_offset + ghosted_maps[i]->GID(j)*dofnums[i];
      for (int dof=0; dof!=dofnums[i]; ++dof) gids[n++] = base+dof;
    }
    global_offset += maps[i]->NumGlobalPoints() * dofnums[i];
    n_global_ghosted += ghosted_maps[i]->NumGlobalPoints() * dofnums[i];
  }
  
  // create the maps.  Note these, unlike the inputs, are Maps, not BlockMaps!
  map_ = Teuchos::rcp(new Epetra_Map(n_global, n_local, &gids[0], 0, *comm)); 
  ghosted_map_ = Teuchos::rcp(new Epetra_Map(n_global_ghosted, n_local_ghosted, &gids[0], 0, *comm));
}


// Copy constructor specifically does not copy indices and ghosted_indices.
// These will be lazily recreated if needed.
SuperMap::SuperMap(const SuperMap& other) :
    compnames_(other.compnames_),
    offsets_(other.offsets_),
    num_dofs_(other.num_dofs_),
    counts_(other.counts_),
    ghosted_offsets_(other.ghosted_offsets_),
    ghosted_counts_(other.ghosted_counts_),
    map_(other.map_),
    ghosted_map_(other.ghosted_map_) {};


bool SuperMap::HasComponent(const std::string& key) const {
  std::map<std::string,int>::const_iterator lb = offsets_.lower_bound(key);
  if (lb != offsets_.end() && !(offsets_.key_comp()(key, lb->first))) {
    return true;
  } else {
    return false;
  }
}
  

const std::vector<int>&
SuperMap::Indices(const std::string& compname, int dofnum) const {
  
  if (indices_.count(compname)) {
    if (indices_[compname].count(dofnum)) {
      return indices_[compname][dofnum];
    }
  }
  return CreateIndices_(compname, dofnum, false);
}


const std::vector<int>&
SuperMap::GhostIndices(const std::string& compname, int dofnum) const {
  if (ghosted_indices_.count(compname)) {
    if (ghosted_indices_[compname].count(dofnum)) {
      return ghosted_indices_[compname][dofnum];
    }
  }
  return CreateIndices_(compname, dofnum, true);
}


std::pair<int, Teuchos::RCP<std::vector<int> > >
SuperMap::BlockIndices() const {
  auto block_indices = Teuchos::rcp(new std::vector<int>(Map()->NumMyElements()));
  int block_id = 0;
  for (const auto& comp : *this) {
    int ndofs = NumDofs(comp);
    for (int d=0; d!=ndofs; ++d) {
      const auto& inds = Indices(comp, d);
      for (int i=0; i!=inds.size(); ++i) (*block_indices)[inds[i]] = block_id;
      block_id++;
    }
  }
  return std::make_pair(block_id, block_indices);
}


const std::vector<int>&
SuperMap::CreateIndices_(const std::string& compname, int dofnum, bool ghosted) const {
  if (ghosted) {
    if (ghosted_indices_.count(compname) == 0) {
      ghosted_indices_[compname];
    }

    // create the vector
    int nentities_owned = counts_.at(compname);
    int nentities = nentities_owned + ghosted_counts_.at(compname);

    std::vector<int> indices(nentities, -1);
    int offset = offsets_.at(compname);
    int num_dof = num_dofs_.at(compname);
    for (int i=0; i!=nentities_owned; ++i) {
      indices[i] = offset + dofnum + i*num_dof;
    }

    int ghosted_offset = ghosted_offsets_.at(compname);
    for (int i=nentities_owned; i!=nentities; ++i) {
      indices[i] = ghosted_offset + dofnum + (i-nentities_owned)*num_dof;
    }

    // assign
    ghosted_indices_[compname][dofnum] = indices;
    return ghosted_indices_[compname][dofnum];

  } else {
    if (indices_.count(compname) == 0) {
      indices_[compname];
    }

    // create the vector
    int nentities = counts_.at(compname);

    std::vector<int> indices(nentities, -1);
    int offset = offsets_.at(compname);
    int num_dof = num_dofs_.at(compname);
    for (int i=0; i!=nentities; ++i) {
      indices[i] = offset + dofnum + i*num_dof;
    }

    // assign
    indices_[compname][dofnum] = indices;

    return indices_[compname][dofnum];
  }
}


} // namespace Operators
} // namespace Amanzi
