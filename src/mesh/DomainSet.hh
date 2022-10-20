/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
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
#include <vector>
#include "Teuchos_RCP.hpp"

#include "errors.hh"
#include "exceptions.hh"
#include "MeshDefs.hh"

class Epetra_MultiVector;

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
            const std::map<std::string, Teuchos::RCP<const std::vector<int>>> maps);

  // iterate over mesh names
  using const_iterator = std::vector<std::string>::const_iterator;
  const_iterator begin() const { return meshes_.begin(); }
  const_iterator end() const { return meshes_.end(); }
  std::size_t size() const { return meshes_.size(); }

  const std::string& get_name() const { return name_; }
  Teuchos::RCP<const Mesh> get_indexing_parent() const { return indexing_parent_; }
  Teuchos::RCP<const Mesh> get_referencing_parent() const { return referencing_parent_; }

  // exporters and maps to/from the parent
  const std::vector<int>& get_subdomain_map(const std::string& subdomain) const {
    if (maps_.size() == 0) {
      Errors::Message msg("DomainSet: subdomain map was requested, but no reference maps were created on construction.");
      Exceptions::amanzi_throw(msg);
    } else if (!maps_.count(subdomain)) {
      Errors::Message msg("DomainSet: subdomain map \"");
      msg << subdomain << "\" requested, but no such subdomain map exists.";
      Exceptions::amanzi_throw(msg);
    }
    return *maps_.at(subdomain);
  }
  void set_subdomain_map(const std::string& subdomain,
                         const Teuchos::RCP<const std::vector<int>>& map) {
    maps_[subdomain] = map;
  }
  const std::map<std::string, Teuchos::RCP<const std::vector<int>>>&
  get_subdomain_maps() const { return maps_; }

  // import from subdomain to parent domain
  void DoImport(const std::string& subdomain,
                const Epetra_MultiVector& src, Epetra_MultiVector& target) const;

  // import from subdomains to parent domain
  void DoImport(const std::vector<const Epetra_MultiVector*>& subdomains,
                Epetra_MultiVector& target) const;

  // import from parent domain to subdomain
  void DoExport(const std::string& subdomain,
                const Epetra_MultiVector& src, Epetra_MultiVector& target) const;
  // import from parent domain to subdomains
  void DoExport(const Epetra_MultiVector& src,
                const std::vector<Epetra_MultiVector*>& targets) const;

 protected:
  std::string name_;
  Teuchos::RCP<const Mesh> indexing_parent_;
  Teuchos::RCP<const Mesh> referencing_parent_;
  std::vector<std::string> meshes_;

  std::map<std::string, Teuchos::RCP<const std::vector<int>>> maps_;
};

//
// helper functions to create importers
//
// Creates an importer from an extracted mesh entity to its parent entities.
Teuchos::RCP<const std::vector<int>>
createMapToParent(const AmanziMesh::Mesh& subdomain_mesh,
                  const AmanziMesh::Entity_kind& src_kind=AmanziMesh::Entity_kind::CELL);

//
// Creates an importer from a surface mesh lifted from an extracted subdomain
// mesh to the global surface mesh (which itself was lifted from the extracted
// subdomain's parent mesh).
Teuchos::RCP<const std::vector<int>>
createMapSurfaceToSurface(const AmanziMesh::Mesh& subdomain_mesh,
        const AmanziMesh::Mesh& parent_mesh);

} // namespace AmanziMesh
} // namespace Amanzi
