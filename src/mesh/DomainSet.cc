/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

//! A collection of meshes, indexed from a parent entity set.

#include "Mesh.hh"
#include "DomainSet.hh"

namespace Amanzi {
namespace AmanziMesh {


DomainSet::DomainSet(const Teuchos::RCP<const Mesh>& parent_,
                     const std::vector<std::string>& regions_,
                     Entity_kind kind_,
                     const std::string& name_,
                     bool by_region_,
                     std::vector<std::string>* region_aliases)
    : parent(parent_),
      regions(regions_),
      kind(kind_),
      name(name_),
      by_region(by_region_)
{
  // construct the collection of mesh names and, if by entity ID, the
  // corresponding ID.
  if (region_aliases == nullptr) region_aliases = &regions;
  AMANZI_ASSERT(region_aliases->size() == regions.size());

  int lcv = 0;
  for (const auto& region : regions) {
    Entity_ID_List ids;
    parent->get_set_entities(region, kind, Parallel_type::OWNED, &ids);
    const auto& map = parent->map(kind, false);

    if (by_region_) {
      // one domain for each region
      if (ids.size() > 0) {
        Key mesh_name = Keys::getDomainInSet(name, (*region_aliases)[lcv]);
        meshes[mesh_name] = lcv;
      }

    } else {
      for (Entity_ID id : ids) {
        Key mesh_name = Keys::getDomainInSet(name, map.GID(id));
        meshes[mesh_name] = id;
      }
    }
    lcv++;
  }
}

} // namespace AmanziMesh
} // namespace Amanzi
