/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)
*/

//! SuperMap class provides a convenient way of creating and using SuperMapLumped

/*
  Amanzi uses SuperMapLumped in a few different ways that make its natural interface
  not that convenient.  SuperMapLumped also has a few limitations that simplify its
  design and implementaiton a great deal.

  This class is a Helper class in that it wraps a SuperMapLumped, providing a
  related interface that is better designed for users.

  It also enforces and mitigates the design limitations of SuperMapLumped itself.
*/

#include "UniqueHelpers.hh"
#include "MeshDefs.hh"
#include "Mesh.hh"
#include "CompositeVectorSpace.hh"
#include "TreeVectorSpace.hh"
#include "TreeVector_Utils.hh"
#include "SuperMap.hh"

namespace Amanzi {
namespace Operators {

// Nonmember contructors/factories
Teuchos::RCP<SuperMap> createSuperMap(const CompositeVectorSpace& cv) {
  return Teuchos::rcp(new SuperMap({cv}));
}


Teuchos::RCP<SuperMap> createSuperMap(const TreeVectorSpace& tvs) {
  if (tvs.Data() != Teuchos::null) {
    // TVS with only a CVS inside
    return createSuperMap(*tvs.Data());
  } else {
    // multiple children
    // grab the leaf nodes
    auto tvss = collectTreeVectorLeaves_const<TreeVectorSpace>(tvs);

    std::vector<CompositeVectorSpace> cvss;
    for (auto tvs_tmp : tvss) {
      cvss.push_back(*tvs_tmp->Data());
    }
    return Teuchos::rcp(new SuperMap(cvss));
  }
}


SuperMap::SuperMap(const std::vector<CompositeVectorSpace>& cvss)
{
  std::vector<std::string> names;
  std::vector<int> dofnums;
  std::vector<Teuchos::RCP<const Epetra_BlockMap> > maps;
  std::vector<Teuchos::RCP<const Epetra_BlockMap> > ghost_maps;

  // this groups maps by map equivalence, not component names.
  
  // loop over nodes, finding unique component names on unique meshes
  int block_num = 0;
  for (auto& cvs : cvss) {
    for (auto compname : cvs) {
      // check if this component's map matches any previously accepted map
      auto this_map = cvs.Map(compname, false);
      int index = std::find_if(maps.begin(), maps.end(),
              [=](const Teuchos::RCP<const Epetra_BlockMap>& m)
              { return this_map->SameBlockMapDataAs(*m); } ) - maps.begin();
      if (index == names.size()) {
        // new map with no previous matches, keep the map
        maps.push_back(this_map);
        ghost_maps.push_back(cvs.Map(compname, true));

        // ensure this name is unique by appending the block_num
        std::string compname_unique = compname+"-"+std::to_string(block_num);
        names.push_back(compname_unique);

        // grab dof count as well
        int this_dofnum = cvs.NumVectors(compname);
        dofnums.push_back(this_dofnum);

        // map the block_num, compname, dof_num tuple to their corresponding
        // values in the SuperMapLumped
        for (int dnum=0; dnum!=this_dofnum; ++dnum) {
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
        AMANZI_ASSERT(cvs.Map(compname, true)->SameBlockMapDataAs(*ghost_maps[index]));

        // map the block_num, compname, dof_num tuple to their corresponding
        // values in the SuperMapLumped
        std::string compname_unique = names[index];
        int this_dofnum = cvs.NumVectors(compname);
        for (int dnum=0; dnum!=this_dofnum; ++dnum) {
          block_info_[std::make_tuple(block_num, compname, dnum)] =
              std::make_pair(compname_unique, dnum+dofnums[index]);
        }

        dofnums[index] += this_dofnum;
      }
    }
    block_num++;
  }

  if (cvss.size() > 0) {
    smap_ = std::make_unique<SuperMapLumped>(cvss[0].Comm(), names, dofnums, maps, ghost_maps);
  }
}

} // namespace Operators
} // namespace Amanzi
