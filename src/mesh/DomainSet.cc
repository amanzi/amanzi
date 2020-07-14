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
                 const std::vector<std::string>& sets_,
                 Entity_kind kind_,
                 const std::string& name_)
    : parent(parent_),
      sets(sets_),
      kind(kind_),
      name(name_)
{
  // construct a list of names by GID
  //
  // NOTE: we do not check, but sets really probably should be non-overlapping?
  for (const auto& set : sets) {
    Entity_ID_List ids;
    parent->get_set_entities(set, kind, Parallel_type::OWNED, &ids);
    const auto& map = parent->map(kind, false);

    for (Entity_ID id : ids) {
      Key mesh_name = name+"_" + std::to_string(map.GID(id));
      meshes[mesh_name] = id;
    }
  }
}  


} // namespace AmanziMesh
} // namespace Amanzi
