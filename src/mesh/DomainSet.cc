/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

//! A collection of meshes, indexed from a parent entity set.

#include "Epetra_MultiVector.h"

#include "Key.hh"
#include "Mesh.hh"
#include "DomainSet.hh"

namespace Amanzi {
namespace AmanziMesh {


DomainSet::DomainSet(const std::string& name,
                     const Teuchos::RCP<const Mesh>& indexing_parent,
                     const std::vector<std::string>& indices)
  : name_(name),
    indexing_parent_(indexing_parent)
{
  // construct the collection of mesh names
  for (const auto& index : indices) {
    meshes_.emplace_back(Keys::getDomainInSet(name, index));
  }
}


DomainSet::DomainSet(const std::string& name,
                     const Teuchos::RCP<const Mesh>& indexing_parent,
                     const std::vector<std::string>& indices,
                     const Teuchos::RCP<const Mesh>& referencing_parent,
                     const std::map<std::string, Teuchos::RCP<const std::vector<int>>> maps)
  : DomainSet(name, indexing_parent, indices)
{
  referencing_parent_ = referencing_parent;
  maps_ = maps;
}


void
DomainSet::doImport(const std::string& subdomain,
                    const Epetra_MultiVector& src, Epetra_MultiVector& target) const
{
  const auto& map = getSubdomainMap(subdomain);

  for (int j=0; j!=src.NumVectors(); ++j) {
    for (int c=0; c!=src.MyLength(); ++c) {
      target[j][map[c]] = src[j][c];
    }
  }
}

void
DomainSet::doExport(const std::string& subdomain,
                    const Epetra_MultiVector& src, Epetra_MultiVector& target) const
{
  const auto& map = getSubdomainMap(subdomain);

  for (int j=0; j!=src.NumVectors(); ++j) {
    for (int c=0; c!=target.MyLength(); ++c) {
      target[j][c] = src[j][map[c]];
    }
  }
}



//
// helper functions to create importers
//
// Creates an importer from an extracted mesh entity to its parent entities.
Teuchos::RCP<const std::vector<int>>
createMapToParent(const AmanziMesh::Mesh& subdomain_mesh,
                  const AmanziMesh::Entity_kind& src_kind)
{
  const Epetra_Map& src_map = subdomain_mesh.getMap(AmanziMesh::Entity_kind::CELL, false);

  // create the target map
  auto map = Teuchos::rcp(new std::vector<int>(src_map.NumMyElements(), -1));
  for (int c=0; c!=src_map.NumMyElements(); ++c) {
    (*map)[c] = subdomain_mesh.getEntityParent(src_kind, c);
  }
  return map;
}

//
// Creates an importer from a surface mesh lifted from an extracted subdomain
// mesh to the global surface mesh (which itself was lifted from the extracted
// subdomain's parent mesh).
Teuchos::RCP<const std::vector<int>>
createMapSurfaceToSurface(const AmanziMesh::Mesh& subdomain_mesh,
        const AmanziMesh::Mesh& parent_mesh)
{
  const Epetra_Map& src_map = subdomain_mesh.getMap(AmanziMesh::Entity_kind::CELL, false);

  int parent_ncells = parent_mesh.getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::OWNED);

  // create the target map
  auto map = Teuchos::rcp(new std::vector<int>(src_map.NumMyElements(), -1));
  for (int c=0; c!=src_map.NumMyElements(); ++c) {
    Entity_ID f = subdomain_mesh.getEntityParent(AmanziMesh::Entity_kind::CELL, c);
    Entity_ID parent_f = subdomain_mesh.getParentMesh()
      ->getEntityParent(AmanziMesh::Entity_kind::FACE, f);

    for (int parent_c=0; parent_c!=parent_ncells; ++parent_c) {
      if (parent_mesh.getEntityParent(AmanziMesh::Entity_kind::CELL, parent_c) == parent_f) {
        (*map)[c] = parent_c;
        break;
      }
    }
    AMANZI_ASSERT((*map)[c] >= 0);
  }
  return map;
}


} // namespace AmanziMesh
} // namespace Amanzi
