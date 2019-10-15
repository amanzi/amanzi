/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
      Daniil Svyatsky (dasvyat@lanl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

/*
  A CompositeSpace is a BlockSpace where all maps come from a specific mesh.

  This should be thought of as a vector-space: it lays out data components as a
  mesh along with entities on the mesh.  This meta-data can be used with the
  mesh's *_map() methods to create the data.

  This class is very light weight as it maintains only meta-data.
*/

#include "errors.hh"

#include "Mesh.hh"
#include "CompositeSpace.hh"

namespace Amanzi {

std::pair<Map_ptr_type, Map_ptr_type>
getMaps(const AmanziMesh::Mesh& mesh, AmanziMesh::Entity_kind location)
{
  return std::make_pair(mesh.map(location, false), mesh.map(location, true));
}


// constructor
CompositeSpace::CompositeSpace(
  const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
  const std::vector<std::string>& names,
  const std::map<std::string, AmanziMesh::Entity_kind> locations,
  const std::map<std::string, BlockMap_ptr_type>& master_maps,
  const std::map<std::string, BlockMap_ptr_type>& ghost_maps,
  const std::map<std::string, std::size_t>& num_vectors, bool ghosted)
  : BlockSpace(mesh->get_comm(), names, master_maps, ghost_maps, num_vectors),
    locations_(locations),
    mesh_(mesh),
    ghosted_(ghosted)
{
  for (const auto& name : names) {
    if (master_maps_.at(name) == Teuchos::null &&
        locations_.at(name) != AmanziMesh::Entity_kind::UNKNOWN) {
      auto maps = getMaps(*mesh, locations.at(name));
      master_maps_[name] = maps.first;
      if (ghosted) ghost_maps_[name] = maps.second;
    }
  }
}

AmanziMesh::Entity_kind
CompositeSpace::Location(const std::string& name) const
{
  if (!HasComponent(name)) {
    Errors::Message message("Map: Requested component (" + name +
                            ") does not exist.");
    throw(message);
  }
  return locations_.at(name);
}


} // namespace Amanzi
