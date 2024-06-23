/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

//! A collection of meshes, indexed from a parent entity set.

#include "Key.hh"
#include "AmanziVector.hh"
#include "Mesh.hh"
#include "MeshCache.hh"
#include "DomainSet.hh"

namespace Amanzi {
namespace AmanziMesh {


DomainSet::DomainSet(const std::string& name,
                     const Teuchos::RCP<const Mesh>& indexing_parent,
                     const std::vector<std::string>& indices)
  : name_(name), indexing_parent_(indexing_parent)
{
  // construct the collection of mesh names
  for (const auto& index : indices) { meshes_.emplace_back(Keys::getDomainInSet(name, index)); }
}


DomainSet::DomainSet(const std::string& name,
                     const Teuchos::RCP<const Mesh>& indexing_parent,
                     const std::vector<std::string>& indices,
                     const Teuchos::RCP<const Mesh>& referencing_parent,
                     const std::map<std::string, MeshCache::cEntity_ID_View> maps)
  : DomainSet(name, indexing_parent, indices)
{
  referencing_parent_ = referencing_parent;
  maps_ = maps;
}


MeshCache::cEntity_ID_View
DomainSet::getSubdomainMap(const std::string& subdomain) const
{
  if (maps_.size() == 0) {
    Errors::Message msg("DomainSet: subdomain map was requested, but no reference maps were "
                        "created on construction.");
    Exceptions::amanzi_throw(msg);
  } else if (!maps_.count(subdomain)) {
    Errors::Message msg("DomainSet: subdomain map \"");
    msg << subdomain << "\" requested, but no such subdomain map exists.";
    Exceptions::amanzi_throw(msg);
  }
  return maps_.at(subdomain);
}

void
DomainSet::setSubdomainMap(const std::string& subdomain, const MeshCache::cEntity_ID_View& map)
{
  maps_[subdomain] = map;
}


//
// helper functions to create importers
//
// Creates an importer from an extracted mesh entity to its parent entities.
MeshCache::Entity_ID_View
createMapToParent(const AmanziMesh::Mesh& subdomain_mesh, const AmanziMesh::Entity_kind& src_kind)
{
  const auto& src_map = subdomain_mesh.getMap(AmanziMesh::Entity_kind::CELL, false);

  // create the target map
  MeshCache::Entity_ID_View map("parent map", src_map->getLocalNumElements());
  const MeshCache& m = subdomain_mesh.getCache();

  Kokkos::parallel_for(
    "DomainSet::createMapToParent", src_map->getLocalNumElements(), KOKKOS_LAMBDA(const int c) {
      map(c) = m.getEntityParent(src_kind, c);
    });
  return map;
}


//
// Creates an importer from a surface mesh lifted from an extracted subdomain
// mesh to the global surface mesh (which itself was lifted from the extracted
// subdomain's parent mesh).
Mesh::Entity_ID_View
createMapSurfaceToSurface(const AmanziMesh::Mesh& subdomain_mesh,
                          const AmanziMesh::Mesh& parent_mesh)
{
  const auto& src_map = subdomain_mesh.getMap(AmanziMesh::Entity_kind::CELL, false);
  int parent_ncells =
    parent_mesh.getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

  // create the target map
  Mesh::Entity_ID_View map("parent map", src_map->getLocalNumElements());

  const MeshCache& subdomain_m = subdomain_mesh.getCache();
  const MeshCache& subdomain_parent_m = subdomain_mesh.getParentMesh()->getCache();
  const MeshCache& parent_m = parent_mesh.getCache();

  Kokkos::parallel_for(
    "DomainSet::createMapSurfaceToSurface",
    src_map->getLocalNumElements(),
    KOKKOS_LAMBDA(const int c) {
      Entity_ID f = subdomain_m.getEntityParent(AmanziMesh::Entity_kind::CELL, c);
      Entity_ID parent_f = subdomain_parent_m.getEntityParent(AmanziMesh::Entity_kind::FACE, f);
      map(c) = -1;
      for (int parent_c = 0; parent_c != parent_ncells; ++parent_c) {
        if (parent_m.getEntityParent(AmanziMesh::Entity_kind::CELL, parent_c) == parent_f) {
          map(c) = parent_c;
          break;
        }
      }
    });
  return map;
}


} // namespace AmanziMesh
} // namespace Amanzi
