/*
  Copyright 2010-201x held jointly by LANL, ORNL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

//! Collects maps and importers for Mesh objects.

#pragma once

#include "AmanziTypes.hh"

namespace Amanzi {
namespace AmanziMesh {

using Map_type = Epetra_Map;
usign Map_ptr_type = Teuchos::RCP<Map_type>;

//
// Creates a pair of maps, <ALL, OWNED>, for a given entity_kind.
// Uses GIDs provided by the Mesh object.
//
template<class Mesh>
std::pair<Map_ptr_type, Map_ptr_type>
createMapsFromMeshGIDs(const Mesh& mesh, const Entity_kind kind);

//
// Creates a pair of maps, <ALL, OWNED>, for a given entity_kind.
// Uses a natural ordering of GIDs, proc 0 == 0...n, proc 1 = n..., etc.
//
template<class Mesh>
std::pair<Map_ptr_type, Map_ptr_type>
createMapsFromNaturalGIDs(const Mesh& mesh, const Entity_kind kind);





template<class Mesh>
class MeshMap {

 public:
  explicit MeshMap(const Mesh& mesh);

  const Epetra_Map& map(Entity_kind kind, bool include_ghost) const;
  const Epetra_Import& importer(Entity_kind kind) const;

  const Epetra_Map& cell_map(bool include_ghost) const {
    return map(Entity_kind::CELL, include_ghost);
  }
  const Epetra_Map& face_map(bool include_ghost) const {
    return map(Entity_kind::FACE, include_ghost);
  }
  const Epetra_Map& edge_map(bool include_ghost) const {
    return map(Entity_kind::EDGE, include_ghost);
  }
  const Epetra_Map& node_map(bool include_ghost) const {
    return map(Entity_kind::NODE, include_ghost);
  }
  const Epetra_Map& boundary_face_map(bool include_ghost) const {
    return map(Entity_kind::BOUNDARY_FACE, include_ghost);
  }
  const Epetra_Map& exterior_node_map(bool include_ghost) const {
    return map(Entity_kind::BOUNDARY_NODE, include_ghost);
  }

  const Epetra_Import& exterior_face_importer() const {
    return *ext_face_importer_;
  }

 private:
  std::map<Entity_kind, Teuchos::RCP<Epetra_Map>> owned_;
  std::map<Entity_kind, Teuchos::RCP<Epetra_Map>> all_;
  std::map<Entity_kind, Teuchos::RCP<Epetra_Import>> importer_;
  Teuchos::RCP<Epetra_Import> ext_face_importer_;
  Teuchos::RCP<Epetra_Import> ext_node_importer_;
}

} // namespace AmanziMesh
} // namespace Amanzi
