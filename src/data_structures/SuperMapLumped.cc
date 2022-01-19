/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon (coonet@ornl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

#include "AmanziComm.hh"
#include "AmanziMap.hh"
#include "AmanziVector.hh"

#include "BlockSpace.hh"
#include "BlockVector.hh"
#include "TreeVectorSpace.hh"
#include "TreeVector_Utils.hh"
#include "SuperMapLumped.hh"

namespace Amanzi {
namespace Operators {

SuperMapLumped::SuperMapLumped(const Teuchos::RCP<const BlockSpace>& space)
  : comp_maps_(space), indices_(std::make_unique<BlockVector<LO>>(space))
{
  CreateIndexing_();
}

SuperMapLumped::~SuperMapLumped() = default;

//
// NOTE: After this is done, CreateIndices_() may be called.
void
SuperMapLumped::CreateIndexing_()
{
  std::cout<<"SuperMapLumped::CreateIndexing_"<<std::endl;
  // fill the offsets
  {
    int offset = 0;
    for (const auto& compname : *comp_maps_) {
      counts_[compname] =
        comp_maps_->ComponentMap(compname, false)->getNodeNumElements();
      ghosted_counts_[compname] =
        comp_maps_->ComponentMap(compname, true)->getNodeNumElements() -
        counts_[compname];
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
  map_ = Teuchos::rcp(new Map_type(n_global, n_local_, 0, indices_->Comm()));

  //
  // We need to know the ghost GIDs of the SuperMap, but the connections must
  // be made to maintain the validity of the connection in the various
  // component maps.  Therefore we use the component maps to scatter out their
  // own GIDs in the supermap to the ghosted supermap.
  BlockVector<GO> gids_comp(comp_maps_);
  std::vector<GO> gids_flat(n_local_ghosted_, -1);

  // fill the LIDs and GIDs of owned entities
  for (const auto& compname : *comp_maps_) {
    #ifdef OUTPUT_CUDA
    std::cout<<"Owned entities of "<<compname<<std::endl;
    #endif
    LO nentities_owned = counts_.at(compname);
    LO nentities = nentities_owned + ghosted_counts_.at(compname);

    indices_->clear_sync_state(compname); 
    gids_comp.clear_sync_state(compname); 

    indices_->modify_host(compname); 
    gids_comp.modify_host(compname);

    auto index_view =
      indices_->ViewComponent<Kokkos::HostSpace>(compname, false);
    auto global_index_view =
      gids_comp.ViewComponent<Kokkos::HostSpace>(compname, false);

    int offset = offsets_.at(compname);
    int n_dofs = comp_maps_->getNumVectors(compname);
    for (LO i = 0; i != nentities_owned; ++i) {
      for (int j = 0; j != n_dofs; ++j) {
        LO lid = offset + j + i * n_dofs;
        GO gid = map_->getGlobalElement(lid);
        index_view(i, j) = lid;
        global_index_view(i, j) = gid;
        gids_flat[lid] = gid;
        #ifdef OUTPUT_CUDA 
        std::cout<<"index_view("<<i<<","<<j<<") = "<<index_view(i,j)<<std::endl;
        std::cout<<"global_index_view("<<i<<","<<j<<") = "<<global_index_view(i,j)<<std::endl;
        #endif 
      }
    }
    //indices_->sync_device(compname); 
    //gids_comp.sync_device(compname);
  }
     
  // communicate the GIDs
  gids_comp.ScatterMasterToGhosted();

  // fill the LIDs and GIDs of the ghost entities
  for (const auto& compname : *comp_maps_) {
    #ifdef OUTPUT_CUDA 
    std::cout<<"Creation loop HOST: "<<compname<<std::endl;
    #endif 
    LO nentities_owned = counts_.at(compname);
    LO nentities = nentities_owned + ghosted_counts_.at(compname);

    auto index_view =
      indices_->ViewComponent<MirrorHost>(compname, true);
    auto global_index_view =
      gids_comp.ViewComponent<MirrorHost>(compname, true);

    int ghosted_offset = ghosted_offsets_.at(compname);
    int n_dofs = comp_maps_->getNumVectors(compname);
    #ifdef OUTPUT_CUDA 
    std::cout<<" nentities-owned: "<<nentities_owned<<" nentities: "<<nentities<<" n_dofs: "<<n_dofs<<std::endl;
    #endif 
    for (LO i = nentities_owned; i != nentities; ++i) {
      for (int j = 0; j != n_dofs; ++j) {
        LO lid = ghosted_offset + j + (i - nentities_owned) * n_dofs;
        index_view(i, j) = lid;
        gids_flat[lid] = global_index_view(i, j);
        #ifdef OUTPUT_CUDA 
        std::cout<<"index_view("<<i<<","<<j<<") = "<<index_view(i,j)<<std::endl;
        std::cout<<"global_index_view("<<i<<","<<j<<") = "<<global_index_view(i,j)<<std::endl;
        #endif 
      }
    }
  }
  #ifdef OUTPUT_CUDA 
  std::cout<<"DONE  Host"<<std::endl;
  #endif 

  // create supermap
  // -- hopefully we found them all!
  if (n_local_ghosted_ > 0) {
    AMANZI_ASSERT(*std::min_element(gids_flat.begin(), gids_flat.end()) >= 0);
  }

  // -- construct
  ghosted_map_ = Teuchos::rcp(new Map_type(
    n_global_ghosted, gids_flat.data(), n_local_ghosted_, 0, indices_->Comm()));


  for (const auto& compname : *comp_maps_) {
    LO nentities_owned = counts_.at(compname);
    LO nentities = nentities_owned + ghosted_counts_.at(compname);
    int n_dofs = comp_maps_->getNumVectors(compname);

    auto index_view =
      indices_->ViewComponent(compname, true);

    #ifdef OUTPUT_CUDA 
    std::cout<<"indices component "<<compname<<" = "; 
    Kokkos::parallel_for(
      "", 
      nentities_owned, 
      KOKKOS_LAMBDA(const int i ){
        for(int j = 0 ; j < n_dofs; ++j){
          printf("(%d/%d) %d - ",i,j,index_view(i,j)); 
        }
      });
    Kokkos::fence(); 
    std::cout<<std::endl; 
    #endif 
  }

}


bool
SuperMapLumped::HasComponent(const std::string& key) const
{
  return indices_->HasComponent(key);
}


} // namespace Operators
} // namespace Amanzi
