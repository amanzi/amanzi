/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

//! A collection of meshes, indexed from a parent entity set.
/*!

A DomainSet is a collection of meshes, generated from a parent mesh and a set of
entities on that mesh.  For instance, if the parent mesh is a 3D, vertically
structured columnar mesh, one might create a set of meshes, each of which are
1D and correspond to cells of the subsurface mesh.  These column meshes might
be indexed from a set of surface faces (or from the surface mesh's cells).

Alternatively, for e.g. a multi-layer snow model, a DomainSet might be a bunch of
mini meshes constructed at each surface cell.

Note, this does not actually know how to construct theses meshes, it only deals
with their names and their relationship to the parent mesh.  To get/store the
actual mesh, this works with State.

*/
#pragma once

#include <string>
#include "Teuchos_RCP.hpp"
#include "Tpetra_Access.hpp"

#include "errors.hh"
#include "exceptions.hh"
#include "AmanziTypes.hh"
#include "MeshDefs.hh"
#include "Mesh.hh"

namespace Amanzi {
namespace AmanziMesh {

//
// A class for iterating over subdomains and storing importers
//
class DomainSet {
 public:
  // this constructor is for domain sets indexed by an entity + regions
  // (e.g. one for each entity in that collection of regions)
  DomainSet(const std::string& name,
            const Teuchos::RCP<const Mesh>& indexing_parent,
            const std::vector<std::string>& indices);

  DomainSet(const std::string& name,
            const Teuchos::RCP<const Mesh>& indexing_parent,
            const std::vector<std::string>& indices,
            const Teuchos::RCP<const Mesh>& referencing_parent,
            const std::map<std::string, Mesh::cEntity_ID_View> maps);

  // iterate over mesh names
  using const_iterator = std::vector<std::string>::const_iterator;
  const_iterator begin() const { return meshes_.begin(); }
  const_iterator end() const { return meshes_.end(); }
  std::size_t size() const { return meshes_.size(); }

  Comm_ptr_type getComm() const { return referencing_parent_->getComm(); }
  const std::string& name() const { return name_; }
  Teuchos::RCP<const Mesh> getIndexingParent() const { return indexing_parent_; }
  Teuchos::RCP<const Mesh> getReferencingParent() const { return referencing_parent_; }

  // exporters and maps to/from the parent
  const Mesh::cEntity_ID_View& getSubdomainMap(const std::string& subdomain) const
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
  void setSubdomainMap(const std::string& subdomain, const Mesh::cEntity_ID_View& map)
  {
    maps_[subdomain] = map;
  }
  const std::map<std::string, Mesh::cEntity_ID_View>& getSubdomainMaps() const { return maps_; }

  // import from subdomain to parent domain
  template <typename scalar_type>
  void doImport(const std::string& subdomain,
                const MultiVector_type_<scalar_type>& src,
                MultiVector_type_<scalar_type>& target) const;

  // import from parent domain to subdomain
  template <typename scalar_type>
  void doExport(const std::string& subdomain,
                const MultiVector_type_<scalar_type>& src,
                MultiVector_type_<scalar_type>& target) const;

 protected:
  std::string name_;
  Teuchos::RCP<const Mesh> indexing_parent_;
  Teuchos::RCP<const Mesh> referencing_parent_;
  std::vector<std::string> meshes_;

  std::map<std::string, Mesh::cEntity_ID_View> maps_;
};

//
// helper functions to create importers
//
// Creates an importer from an extracted mesh entity to its parent entities.
Mesh::Entity_ID_View
createMapToParent(const AmanziMesh::Mesh& subdomain_mesh,
                  const AmanziMesh::Entity_kind& src_kind = AmanziMesh::Entity_kind::CELL);

//
// Creates an importer from a surface mesh lifted from an extracted subdomain
// mesh to the global surface mesh (which itself was lifted from the extracted
// subdomain's parent mesh).
Mesh::Entity_ID_View
createMapSurfaceToSurface(const AmanziMesh::Mesh& subdomain_mesh,
                          const AmanziMesh::Mesh& parent_mesh);


// import from subdomain to parent domain
template <typename scalar_type>
void
DomainSet::doImport(const std::string& subdomain,
                    const MultiVector_type_<scalar_type>& src_vec,
                    MultiVector_type_<scalar_type>& target_vec) const
{
  const auto& map = getSubdomainMap(subdomain);
  auto src = src_vec.getLocalViewDevice(Tpetra::Access::ReadOnly);
  auto target = target_vec.getLocalViewDevice(Tpetra::Access::ReadWrite);
  Kokkos::parallel_for(
    "DomainSet::doImport",
    Kokkos::MDRangePolicy<Kokkos::Rank<2>>({ 0, 0 }, { src.extent(0), src.extent(1) }),
    KOKKOS_LAMBDA(const int i, const int j) {
      assert(map(i) >= 0 && map(i) < target.extent(0));
      target(map(i), j) = src(i, j);
    });
}

// import from parent domain to subdomain
template <typename scalar_type>
void
DomainSet::doExport(const std::string& subdomain,
                    const MultiVector_type_<scalar_type>& src_vec,
                    MultiVector_type_<scalar_type>& target_vec) const
{
  const auto& map = getSubdomainMap(subdomain);
  auto src = src_vec.getLocalViewDevice(Tpetra::Access::ReadOnly);
  auto target = target_vec.getLocalViewDevice(Tpetra::Access::ReadWrite);
  Kokkos::parallel_for(
    "DomainSet::doExport",
    Kokkos::MDRangePolicy<Kokkos::Rank<2>>({ 0, 0 }, { target.extent(0), target.extent(1) }),
    KOKKOS_LAMBDA(const int i, const int j) { target(i, j) = src(map(i), j); });
}


} // namespace AmanziMesh
} // namespace Amanzi
