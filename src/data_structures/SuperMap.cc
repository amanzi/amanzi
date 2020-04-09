/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon (coonet@ornl.gov)
*/

//! SuperMap class provides a convenient way of creating and using
//! SuperMapLumped

/*
  Amanzi uses SuperMapLumped in a few different ways that make its natural
  interface not that convenient.  SuperMapLumped also has a few limitations that
  simplify its design and implementaiton a great deal.

  This class is a Helper class in that it wraps a SuperMapLumped, providing a
  related interface that is better designed for users.

  It also enforces and mitigates the design limitations of SuperMapLumped
  itself.
*/

#include "AmanziVector.hh"

#include "UniqueHelpers.hh"
#include "MeshDefs.hh"
#include "Mesh.hh"
#include "BlockSpace.hh"
#include "BlockVector.hh"
#include "TreeVectorSpace.hh"
#include "TreeVector_Utils.hh"
#include "TreeVector.hh"
#include "SuperMap.hh"

namespace Amanzi {
namespace Operators {

// Nonmember contructors/factories
Teuchos::RCP<SuperMap>
createSuperMap(const Teuchos::Ptr<const BlockSpace>& cv)
{
  return Teuchos::rcp(new SuperMap(cv->Comm(), { cv }));
}

Teuchos::RCP<SuperMap>
createSuperMap(const TreeVectorSpace& tvs)
{
  if (tvs.Data() != Teuchos::null) {
    // TVS with only a CVS inside
    return createSuperMap(tvs.Data().ptr());
  } else {
    // multiple children
    // grab the leaf nodes
    auto tvss = collectTreeVectorLeaves_const<TreeVectorSpace>(tvs);

    std::vector<Teuchos::Ptr<const BlockSpace>> cvss;
    for (auto a_tvs : tvss) { cvss.push_back(a_tvs->Data().ptr()); }
    return Teuchos::rcp(new SuperMap(tvs.Comm(), cvss));
  }
}


// Copy in/out
int
copyToSuperVector(const SuperMap& map, const BlockVector<Vector_type::scalar_type>& bv,
                  Vector_type& sv, int block_num)
{
  auto svv = sv.getLocalViewDevice();
  for (const auto& compname : bv) {
    if (map.HasComponent(block_num, compname)) {
      auto data = bv.ViewComponent(compname, false);
      for (int dofnum = 0; dofnum != bv.getNumVectors(compname); ++dofnum) {
        auto inds = map.Indices(block_num, compname, dofnum);
        Kokkos::parallel_for(
          "SuperMap::copyToSuperVector",
          data.extent(0),
                             KOKKOS_LAMBDA(const int i) {
                               svv(inds(i),0) = data(i,dofnum);
                             });
      }
    }
  }
  return 0;
}

int
copyFromSuperVector(const SuperMap& map, const Vector_type& sv,
                    BlockVector<Vector_type::scalar_type>& bv, int block_num)
{
  auto svv = sv.getLocalViewDevice();
  for (const auto& compname : bv) {
    if (map.HasComponent(block_num, compname)) {
      auto data = bv.ViewComponent(compname, false);
      for (int dofnum = 0; dofnum != bv.getNumVectors(compname); ++dofnum) {
        auto inds = map.Indices(block_num, compname, dofnum);
        Kokkos::parallel_for(
          "SuperMap::copyFromSuperVector",
          data.extent(0),
                             KOKKOS_LAMBDA(const int i) {
                               data(i,dofnum) = svv(inds(i),0);
                             });
      }
    }
  }  
  return 0;
}

int
addFromSuperVector(const SuperMap& map, const Vector_type& sv,
                   BlockVector<Vector_type::scalar_type>& bv, int block_num)
{
  auto svv = sv.getLocalViewDevice();
  for (const auto& compname : bv) {
    if (map.HasComponent(block_num, compname)) {
      auto data = bv.ViewComponent(compname, false);
      for (int dofnum = 0; dofnum != bv.getNumVectors(compname); ++dofnum) {
        auto inds = map.Indices(block_num, compname, dofnum);
        Kokkos::parallel_for(
          "SuperMap::addFromSuperVector",
          data.extent(0),
                             KOKKOS_LAMBDA(const int i) {
                               data(i,dofnum) += svv(inds(i),0);
                             });
      }
    }
  }  
  return 0;
}

// Nonmember TreeVector to/from Super-vector
// -- simple schema version
int
copyToSuperVector(const SuperMap& map, const TreeVector& tv,
                  Vector_type& sv)
{
  int ierr(0);

  if (tv.Data().get()) {
    return copyToSuperVector(map, *tv.Data(), sv, 0);
  } else {
    auto sub_tvs = collectTreeVectorLeaves_const(tv);
    int block_num = 0;
    for (const auto& sub_tv : sub_tvs) {
      ierr |= copyToSuperVector(map, *sub_tv->Data(), sv, block_num);
      block_num++;
    }
    return ierr;
  }
}

int
copyFromSuperVector(const SuperMap& map, const Vector_type& sv,
                    TreeVector& tv)
{
  int ierr(0);

  if (tv.Data().get()) {
    return copyFromSuperVector(map, sv, *tv.Data(), 0);
  } else {
    auto sub_tvs = collectTreeVectorLeaves(tv);
    int block_num = 0;
    for (auto& sub_tv : sub_tvs) {
      ierr |= copyFromSuperVector(map, sv, *sub_tv->Data(), block_num);
      block_num++;
    }
    return ierr;
  }
}

int
addFromSuperVector(const SuperMap& map, const Vector_type& sv,
                   TreeVector& tv)
{
  int ierr(0);

  if (tv.Data().get()) {
    return addFromSuperVector(map, sv, *tv.Data(), 0);
  } else {
    auto sub_tvs = collectTreeVectorLeaves(tv);
    int block_num = 0;
    for (auto& sub_tv : sub_tvs) {
      ierr |= addFromSuperVector(map, sv, *sub_tv->Data(), block_num);
      block_num++;
    }
    return ierr;
  }
}



SuperMap::SuperMap(const Comm_ptr_type& comm,
                   const std::vector<Teuchos::Ptr<const BlockSpace>>& cvss)
{
  std::vector<std::string> names;
  std::map<std::string, std::size_t> dofnums;
  std::map<std::string, BlockMap_ptr_type> master_maps;
  std::map<std::string, BlockMap_ptr_type> ghost_maps;

  // this groups maps by map equivalence, not component names.

  // loop over nodes, finding unique component names on unique meshes
  std::size_t block_num = 0;
  for (auto cvs : cvss) {
    for (const auto& compname : *cvs) {
      // check if this component's map matches any previously accepted map
      auto this_map = cvs->ComponentMap(compname, false);
      auto index_iter = std::find_if(
        master_maps.begin(),
        master_maps.end(),
        [=](const std::pair<const std::string, const BlockMap_ptr_type>& m) {
          return this_map->locallySameAs(*m.second);
        });
      if (index_iter == master_maps.end()) {
        // new map with no previous matches, keep the map
        // -- ensure this name is unique by appending the block_num
        std::string compname_unique =
          compname + "-" + std::to_string(block_num);
        names.push_back(compname_unique);
        master_maps[compname_unique] = this_map;
        ghost_maps[compname_unique] = cvs->ComponentMap(compname, true);

        std::size_t this_dofnum = cvs->getNumVectors(compname);
        dofnums[compname_unique] = this_dofnum;

        // map the block_num, compname, dof_num tuple to their corresponding
        // values in the SuperMapLumped
        for (std::size_t dnum = 0; dnum != this_dofnum; ++dnum) {
          block_info_[std::make_tuple(block_num, compname, dnum)] =
            std::make_pair(compname_unique, dnum);
        }

      } else {
        // previously found this map.
        //
        // If the master map is the same, the ghost map had better be the same
        // too?  If this assert ever throws, it can be added to the above
        // lambda, to ensure that both the owned map and the ghosted map are
        // the same.
        AMANZI_ASSERT(cvs->ComponentMap(compname, true)
                        ->locallySameAs(*ghost_maps[index_iter->first]));

        // map the block_num, compname, dof_num tuple to their corresponding
        // values in the SuperMapLumped
        std::size_t this_dofnum = cvs->getNumVectors(compname);
        for (std::size_t dnum = 0; dnum != this_dofnum; ++dnum) {
          block_info_[std::make_tuple(block_num, compname, dnum)] =
            std::make_pair(index_iter->first,
                           dnum + dofnums[index_iter->first]);
        }
        dofnums[index_iter->first] += this_dofnum;
      }
    }
    block_num++;
  }

  if (cvss.size() > 0) {
    auto block_map = Teuchos::rcp(
      new BlockSpace(comm, names, master_maps, ghost_maps, dofnums));
    smap_ = std::make_unique<SuperMapLumped>(block_map);
  }
}

} // namespace Operators
} // namespace Amanzi
