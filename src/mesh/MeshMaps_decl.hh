/*
  Copyright 2010-201x held jointly by LANL, ORNL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

//! Collects maps and importers for Mesh objects.

/*

  Maps are indexing objects, a way of going from the "local ID", or LID, to
  "global IDs" or GIDs.

  LIDs are a non-negative integer from [0,num_entities_owned_on_this_rank)

  GIDs may be any non-negative integer, though typically are in
  [0,num_entities_globally).

  A map, then is effectively an array of length num_entities_on_this_process
  indexed by LID, containing the corresponding GID.

  Boundary entity maps are somewhat different -- a boundary LID corresponds to
  the GID of the corresponding entity, so they are not in
  [0,num_global_boundary_entities), but in [0,num_global_entities).

  Similarly, in derived meshes, the GID is often the GID of the parent entity.

  Importers are used to provide scatter/gather halo exchange operations across
  processors.  Importers implement the scatter operation; so for instance, the
  following examples are valid code:

  ..code::
      // set up a ghosted vector
      auto& cell_map maps.getMap(Entity_kind::CELL, true);
      Vector_type vec(cell_map);
      // ... set vec data somehow...

      // scatter operation
      vec.Import(vec, maps.getImporter(Entity_kind::CELL), Insert);

      // gather operation
      vec.Export(vec, maps.getImporter(Entity_kind::CELL), Add);

  Note this is just as valid with BOUNDARY entities as with normal entities.

  Additionally, importers are provided to import from a FACE vector to a
  BOUNDARY_FACE vector (respectively NODE and BOUNDARY_NODE):


  ..code::
      // set up a face vector
      auto& face_map maps.getMap(Entity_kind::FACE, true);
      Vector_type face_vec(face_map);
      // ... set face_vec data somehow...

      // set up a boundary_face vector
      auto& bf_map maps.getMap(Entity_kind::BOUNDARY_FACE, true);
      Vector_type bf_vec(bf_map);

      // import operation
      bf_vec.Import(face_vec, maps.getBoundaryFaceImporter(), Insert);

      // export operation
      face_vec.Export(bf_vec, maps.getBoundaryFaceImporter(), Insert);

*/

#pragma once

#include "AmanziTypes.hh"
#include "ViewUtils.hh"

namespace Amanzi {
namespace AmanziMesh {

//
// Creates a pair of maps, <ALL, OWNED>, for a given entity_kind.
// Uses GIDs provided by the Mesh object.
//
template<class Mesh_type>
std::pair<Map_ptr_type, Map_ptr_type>
createMapsFromMeshGIDs(const Mesh_type& mesh, const Entity_kind kind);

//
// Creates a pair of maps, <ALL, OWNED>, for a given entity_kind.
// Uses a natural ordering of GIDs, proc 0 == 0...n, proc 1 = n..., etc.
//
template<class Mesh_type>
std::pair<Map_ptr_type, Map_ptr_type>
createMapsFromContiguousGIDs(const Mesh_type& mesh, const Entity_kind kind);


class MeshMaps {

 public:
  MeshMaps() {}

  // If renumber is True, map GIDs are renumbered contiguously, starting on
  // rank 0 with 0...num_owned_entities and continuing through the ranks.
  template<class Mesh_type>
  void initialize(const Mesh_type& mesh,
                  bool renumber=false);

  template<MemSpace_kind MEM>
  decltype(auto) // cEntity_ID_View
  getBoundaryFaces() const {
    return view<MEM>(boundary_faces_);
  }

  template<MemSpace_kind MEM>
  decltype(auto) // cEntity_ID_View
  getBoundaryNodes() const {
    return view<MEM>(boundary_nodes_);
  }

  std::size_t 
  getNBoundaryFaces(Parallel_kind ptype) {
    if(Parallel_kind::OWNED == ptype)
      return owned_[Entity_kind::BOUNDARY_FACE]->NumMyElements(); 
    else 
      return all_[Entity_kind::BOUNDARY_FACE]->NumMyElements(); 
  }

  std::size_t 
  getNBoundaryNodes(Parallel_kind ptype) {
    if(Parallel_kind::OWNED == ptype)
      return owned_[Entity_kind::BOUNDARY_NODE]->NumMyElements(); 
    else 
      return all_[Entity_kind::BOUNDARY_NODE]->NumMyElements(); 
  }

  inline const Map_type& getMap(Entity_kind kind, bool include_ghost) const;
  inline const Import_type& getImporter(Entity_kind kind) const;
  const Import_type& getBoundaryFaceImporter() const {
    return *boundary_face_importer_;
  }
  const Import_type& getBoundaryNodeImporter() const {
    return *boundary_node_importer_;
  }

 private:
  std::map<Entity_kind, Map_ptr_type> owned_;
  std::map<Entity_kind, Map_ptr_type> all_;
  std::map<Entity_kind, Import_ptr_type> importer_;
  Import_ptr_type boundary_face_importer_;
  Import_ptr_type boundary_node_importer_;
  Entity_ID_DualView boundary_faces_;
  Entity_ID_DualView boundary_nodes_;
};


} // namespace AmanziMesh
} // namespace Amanzi
