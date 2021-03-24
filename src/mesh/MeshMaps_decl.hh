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


class MeshMaps {

 public:
  MeshMaps() {}

  template<class Mesh>
  void initialize(const Mesh& mesh,
                  bool natural_ordering=false,
                  bool request_edges=false);

  const Entity_ID_View& get_boundary_faces() const {
    return boundary_faces_;
  }
  const Entity_ID_View& get_boundary_nodes() const {
    return boundary_nodes_;
  }
  const Map_type& get_map(Entity_kind kind, bool include_ghost) const;
  const Import_type& get_importer(Entity_kind kind) const;
  const Import_type& get_boundary_face_importer() const {
    return *boundary_face_importer_;
  }
  const Import_type& get_boundary_node_importer() const {
    return *boundary_node_importer_;
  }

 private:
  std::map<Entity_kind, Map_ptr_type> owned_;
  std::map<Entity_kind, Map_ptr_type> all_;
  std::map<Entity_kind, Import_ptr_type> importer_;
  Import_ptr_type boundary_face_importer_;
  Import_ptr_type boundary_node_importer_;
  Entity_ID_View boundary_faces_;
  Entity_ID_View boundary_nodes_;
};

} // namespace AmanziMesh
} // namespace Amanzi
