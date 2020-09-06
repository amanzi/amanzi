/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Rao Garimella, others

  Abstract class that adds sets and maps to the mini mesh.
*/

#ifndef AMANZI_MESH_STATIC_MESH_HH_
#define AMANZI_MESH_STATIC_MESH_HH_

#include "Epetra_Import.h"

// Amanzi::AmanziMesh
#include "GeometricModel.hh"
#include "MeshDefs.hh"
#include "StaticMeshMini.hh"

namespace Amanzi {
namespace AmanziMesh {

class StaticMesh : public StaticMeshMini {
 public:
  StaticMesh()
    : comm_(NULL),
      parent_(Teuchos::null),
      manifold_dim_(-1),
      logical_(false) {};

  Comm_ptr_type get_comm() const { return comm_; }

  //--------------------------------
  // Additional downward adjacencies
  //--------------------------------
  // Get the local index of a face edge in a cell edge list
  //
  // Example:
  // face_get_edges(face=5) --> {20, 21, 35, 9, 10}
  // cell_get_edges(cell=18) --> {1, 2, 3, 5, 8, 9, 10, 13, 21, 35, 20, 37, 40}
  // face_to_cell_edge_map(face=5,cell=18) --> {10, 8, 9, 5, 6}
  virtual void face_to_cell_edge_map(
          const Entity_ID f, const Entity_ID c, std::vector<int> *map) const = 0;


  //------------------------------
  // Additional upward adjacencies
  //------------------------------
  // Cells of type 'ptype' connected to a node
  // NOTE: The order of cells is not guaranteed to be the same for
  // corresponding nodes on different processors
  virtual void node_get_cells(
          const Entity_ID nodeid,
          const Parallel_type ptype,
          Entity_ID_List *cellids) const = 0;

  // Faces of type parallel 'ptype' connected to a node
  // NOTE: The order of faces is not guarnateed to be the same for
  // corresponding nodes on different processors
  virtual void node_get_faces(
          const Entity_ID nodeid,
          const Parallel_type ptype,
          Entity_ID_List *faceids) const = 0;

  // Faces of type 'ptype' connected to an edge
  // NOTE: The order of faces is not guaranteed to be the same for
  // corresponding edges on different processors
  virtual void edge_get_faces(
          const Entity_ID edgeid,
          const Parallel_type ptype,
          Entity_ID_List *faceids) const { AMANZI_ASSERT(false); }

  // Cells connected to a face
  // The cells are returned in no particular order. Also, the order of cells
  // is not guaranteed to be the same for corresponding faces on different
  // processors
  virtual void face_get_cells(
          const Entity_ID f,
          const Parallel_type ptype,
          Entity_ID_List *cells) const = 0;

  // ------------------------
  // Additional mesh topology
  // ------------------------
  // On a distributed mesh, this will count all entities, OWNED or GHOST.
  virtual unsigned int cell_get_num_faces(Entity_ID cellid) const = 0;


  // ------------------------
  // Additional mesh geometry
  // ------------------------
  // Node coordinates - 3 in 3D and 2 in 2D
  virtual void node_get_coordinates(
          const Entity_ID v, AmanziGeometry::Point *ncoord) const = 0;

  // Face coordinates - conventions same as face_to_nodes call
  // Number of nodes is the vector size.
  virtual void face_get_coordinates(
          const Entity_ID f, std::vector<AmanziGeometry::Point> *fcoords) const = 0;

  // Coordinates of cells in standard order (Exodus II convention)
  //
  // NOTE: Standard convention works only for standard cell types in 3D!
  // For a general polyhedron this will return the node coordinates in
  // arbitrary order.
  virtual void cell_get_coordinates(
          const Entity_ID c, std::vector<AmanziGeometry::Point> *ccoords) const = 0;


  //------------
  // Epetra maps
  //------------
  virtual const Epetra_Map& cell_map(bool include_ghost) const = 0;
  virtual const Epetra_Map& face_map(bool include_ghost) const = 0;
  virtual const Epetra_Map& node_map(bool include_ghost) const = 0;

  virtual const Epetra_Map& exterior_face_map(bool include_ghost) const = 0;

  // dummy implementation so that frameworks can skip or overwrite
  virtual const Epetra_Map& edge_map(bool include_ghost) const
  {
    Errors::Message mesg("Edges are not implemented in this framework.");
    Exceptions::amanzi_throw(mesg);
    throw(mesg);
  }

  const Epetra_Map& map(Entity_kind kind, bool include_ghost) const {
    if (kind == CELL) return cell_map(include_ghost);
    else if (kind == FACE) return face_map(include_ghost);
    else if (kind == EDGE) return edge_map(include_ghost);
    else if (kind == NODE) return node_map(include_ghost);
    else if (kind == BOUNDARY_FACE) return exterior_face_map(include_ghost);
    Errors::Message mesg("No such map type.");
    Exceptions::amanzi_throw(mesg);
    throw(mesg);
  }

  // Epetra importers that will allow apps to import values from a
  // Epetra vector defined on all owned faces into an Epetra vector
  // defined only on exterior faces
  virtual const Epetra_Import& exterior_face_importer() const = 0;
  virtual const Epetra_Map& exterior_node_map(bool include_ghost) const = 0;


  // ------------------------
  // General mesh information
  // ------------------------
  // Geometric model
  void set_geometric_model(const Teuchos::RCP<const AmanziGeometry::GeometricModel>& gm) { geometric_model_ = gm; }
  Teuchos::RCP<const AmanziGeometry::GeometricModel> geometric_model() const { return geometric_model_; }

  // Number of entities of any kind (cell, face, node) and in a
  // particular category (OWNED, GHOST, ALL)
  virtual unsigned int num_entities(
          const Entity_kind kind, const Parallel_type ptype) const = 0;

  // Were optional edges initialized?
  virtual bool valid_edges() const { return false; }

  // Usually this is null, but derived meshes may provide it.
  virtual Teuchos::RCP<const StaticMesh> parent() const { return parent_; }

  void set_manifold_dimension(const unsigned int dim) { manifold_dim_ = dim; }
  unsigned int manifold_dimension() const { return manifold_dim_; }

  bool is_logical() const { return logical_; }


  //--------------------------------------------------------------
  // Mesh Sets for ICs, BCs, Material Properties and whatever else
  //--------------------------------------------------------------
  // Is this is a valid ID of a set containing entities of 'kind'
  bool valid_set_id(const Set_ID setid, const Entity_kind kind) const;
  bool valid_set_name(const std::string setname, const Entity_kind kind) const;

  // Get number of entities of type 'category' in set
  virtual unsigned int get_set_size(
          const std::string& setname,
          const Entity_kind kind,
          const Parallel_type ptype) const;

  // Get list of entities of type 'category' in set by set name.
  // -- original interface. The returned volume fractions are ignored.
  virtual void get_set_entities(
          const std::string& setname,
          const Entity_kind kind,
          const Parallel_type ptype,
          Entity_ID_List *entids) const;

  // Since not all regions support volume fractions 
  // (vofs), this vector is optional and could be empty.
  virtual void get_set_entities_and_vofs(
          const std::string& setname,
          const Entity_kind kind,
          const Parallel_type ptype,
          Entity_ID_List *entids,
          std::vector<double> *vofs) const = 0;

 protected:
  void get_set_entities_box_vofs_(
      Teuchos::RCP<const AmanziGeometry::Region> region,
      const Entity_kind kind, 
      const Parallel_type ptype, 
      std::vector<Entity_ID>* setents,
      std::vector<double> *vofs) const;

 protected:
  Comm_ptr_type comm_;

  Teuchos::RCP<const AmanziGeometry::GeometricModel> geometric_model_;

  // miscallenous mesh information
  bool logical_;
  unsigned int manifold_dim_;
  Teuchos::RCP<const StaticMesh> parent_;

  // -- region data
  mutable std::map<std::string, std::vector<int> > region_ids;
  mutable std::map<std::string, std::vector<double> > region_vofs;
};

}  // namespace AmanziMesh
}  // namespace Amanzi

#endif

