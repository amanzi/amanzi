/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

//! This is a mesh for a single surface cell.

/*!

  This exists solely because we need "surface meshes" extracted from
  a MeshColumn.  This is really just 1 2D cell.  Really.

*/

#ifndef AMANZI_MESH_SURFACE_CELL_HH_
#define AMANZI_MESH_SURFACE_CELL_HH_

#include <vector>
#include <string>
#include <algorithm>

#include "Teuchos_ParameterList.hpp"
#include "Epetra_Map.h"
#include "AmanziComm.hh"
#include "Epetra_SerialComm.h"

#include "VerboseObject.hh"
#include "dbc.hh"
#include "errors.hh"

#include "Region.hh"
#include "Mesh.hh"

namespace Amanzi {
namespace AmanziMesh {

class MeshSurfaceCell : public Mesh {
 public:
  MeshSurfaceCell(const Teuchos::RCP<const Mesh>& parent_mesh,
                  const std::string& setname,
                  bool flatten=true);

  virtual ~MeshSurfaceCell() = default;

  // Get parallel type of entity - OWNED, GHOST, ALL (See MeshDefs.hh)
  virtual
  Parallel_type entity_get_ptype(const Entity_kind kind,
          const Entity_ID entid) const override {
    return Parallel_type::OWNED;
  }

  // Parent entity in the source mesh if mesh was derived from another mesh
  virtual
  Entity_ID entity_get_parent(const Entity_kind kind, const Entity_ID entid) const override {
    AMANZI_ASSERT(kind == CELL);
    AMANZI_ASSERT(entid == 0);
    return parent_face_;
  }

  // Get cell type - UNKNOWN, TRI, QUAD, ... See MeshDefs.hh
  virtual
  Cell_type cell_get_type(const Entity_ID cellid) const override {
    return cell_type_;
  }


  //
  // General mesh information
  // -------------------------
  //

  // Number of entities of any kind (cell, face, node) and in a
  // particular category (OWNED, GHOST, ALL)
  virtual
  unsigned int num_entities(const Entity_kind kind,
                            const Parallel_type ptype) const override;

  // Global ID of any entity
  virtual
  Entity_ID GID(const Entity_ID lid, const Entity_kind kind) const override {
    return lid;
  }


  //
  // Mesh Entity Adjacencies
  //-------------------------

  // Downward Adjacencies
  //---------------------

  // Get nodes of a cell
  virtual
  void cell_get_nodes(const Entity_ID cellid,
                      Entity_ID_List *nodeids) const override;

  // Get nodes of face
  // On a distributed mesh, all nodes (OWNED or GHOST) of the face
  // are returned
  // In 3D, the nodes of the face are returned in ccw order consistent
  // with the face normal
  // In 2D, nfnodes is 2
  virtual
  void face_get_nodes(const Entity_ID faceid,
                      Entity_ID_List *nodeids) const override;

  // Get nodes of edge
  virtual
  void edge_get_nodes(const Entity_ID edgeid,
                      Entity_ID *nodeid0, Entity_ID *nodeid1) const override;

  // Upward adjacencies
  //-------------------

  // Cells of type 'ptype' connected to a node - The order of cells is
  // not guaranteed to be the same for corresponding nodes on
  // different processors
  virtual
  void node_get_cells(const Entity_ID nodeid,
                      const Parallel_type ptype,
                      Entity_ID_List *cellids) const override;

  // Faces of type 'ptype' connected to a node - The order of faces is
  // not guarnateed to be the same for corresponding nodes on
  // different processors
  virtual
  void node_get_faces(const Entity_ID nodeid,
                      const Parallel_type ptype,
                      Entity_ID_List *faceids) const override;

  // Cells of type 'ptype' connected to an edges
  virtual
  void edge_get_cells(const Entity_ID edgeid,
                      const Parallel_type ptype,
                      Entity_ID_List *cellids) const override;

  // Same level adjacencies
  //-----------------------

  // Face connected neighboring cells of given cell of a particular ptype
  // (e.g. a hex has 6 face neighbors)

  // The order in which the cellids are returned cannot be
  // guaranteed in general except when ptype = ALL, in which case
  // the cellids will correcpond to cells across the respective
  // faces given by cell_get_faces
  virtual
  void cell_get_face_adj_cells(const Entity_ID cellid,
          const Parallel_type ptype,
          Entity_ID_List *fadj_cellids) const override ;


  //
  // Mesh entity geometry
  //--------------
  //

  // Node coordinates - 3 in 3D and 2 in 2D
  virtual
  void node_get_coordinates(const Entity_ID nodeid,
                            AmanziGeometry::Point *ncoord) const override {
    (*ncoord) = nodes_[nodeid];
  }


  // Face coordinates - conventions same as face_to_nodes call
  // Number of nodes is the vector size divided by number of spatial dimensions
  virtual
  void face_get_coordinates(const Entity_ID faceid,
                            std::vector<AmanziGeometry::Point> *fcoords) const override {
    fcoords->resize(2);
    (*fcoords)[0] = nodes_[faceid];
    (*fcoords)[1] = nodes_[(faceid + 1) % nodes_.size()];
  }

  // Coordinates of cells in standard order (Exodus II convention)
  // STANDARD CONVENTION WORKS ONLY FOR STANDARD CELL TYPES IN 3D
  // For a general polyhedron this will return the node coordinates in
  // arbitrary order
  // Number of nodes is vector size divided by number of spatial dimensions
  virtual
  void cell_get_coordinates(const Entity_ID cellid,
                            std::vector<AmanziGeometry::Point> *ccoords) const override {
    (*ccoords) = nodes_;
  }


  //
  // Mesh modification
  //-------------------

  // Set coordinates of node
  virtual
  void node_set_coordinates(const Entity_ID nodeid,
                            const AmanziGeometry::Point ncoord) override {
    nodes_[nodeid] = ncoord;
  }


  virtual
  void node_set_coordinates(const Entity_ID nodeid,
                            const double *ncoord) override {
    Errors::Message mesg("Not implemented");
    Exceptions::amanzi_throw(mesg);
  }


  // Deform a mesh so that cell volumes conform as closely as possible
  // to target volumes without dropping below the minimum volumes.  If
  // move_vertical = true, nodes will be allowed to move only in the
  // vertical direction (right now arbitrary node movement is not allowed)
  // Nodes in any set in the fixed_sets will not be permitted to move.
  virtual
  int deform(const std::vector<double>& target_cell_volumes_in,
             const std::vector<double>& min_cell_volumes_in,
             const Entity_ID_List& fixed_nodes,
             const bool move_vertical) override;
  //
  // Epetra maps
  //------------
  virtual
  const Epetra_Map& cell_map(bool include_ghost) const override {
    return *cell_map_;
  }

  virtual
  const Epetra_Map& face_map(bool include_ghost) const override {
    return *face_map_;
  }

  // dummy implementation so that frameworks can skip or overwrite
  const Epetra_Map& edge_map(bool include_ghost) const override {
    Errors::Message mesg("Edges not implemented in this framework");
    Exceptions::amanzi_throw(mesg);
    throw(mesg);
  };

  virtual
  const Epetra_Map& node_map(bool include_ghost) const override {
    return *face_map_;
  }

  virtual
  const Epetra_Map& exterior_face_map(bool include_ghost) const override {
    return *face_map_;
  }

  virtual
  const Epetra_Map& exterior_node_map(bool include_ghost) const override {
    Errors::Message mesg("Exterior node map is not implemented in this framework");
    Exceptions::amanzi_throw(mesg);
    throw(mesg);
  }


  // Epetra importer that will allow apps to import values from a
  // Epetra vector defined on all owned faces into an Epetra vector
  // defined only on exterior faces
  virtual
  const Epetra_Import& exterior_face_importer(void) const override {
    return *exterior_face_importer_;
  }


  //
  // Mesh Sets for ICs, BCs, Material Properties and whatever else
  //--------------------------------------------------------------
  //

  // Get number of entities of type 'category' in set
  virtual
  unsigned int get_set_size(const std::string& setname,
                            const Entity_kind kind,
                            const Parallel_type ptype) const override;

  virtual
  bool valid_set_type(const AmanziGeometry::RegionType rtype, const Entity_kind kind) const override;

  // Get list of entities of type 'category' in set
  virtual
  void get_set_entities_and_vofs(const std::string& setname,
          const Entity_kind kind,
          const Parallel_type ptype,
          Entity_ID_List *entids,
          std::vector<double> *vofs) const override;


  // Miscellaneous functions
  virtual
  void write_to_exodus_file(const std::string filename) const override;

 protected:

  // get faces and face dirs of a cell. This can be called by
  // cell_get_faces_and_dirs method of the base class and the data
  // cached or it can be called directly by the
  // cell_get_faces_and_dirs method of this class
  virtual
  void cell_get_faces_and_dirs_internal_(const Entity_ID cellid,
          Entity_ID_List *faceids,
          std::vector<int> *face_dirs,
          const bool ordered=false) const override;

  // Cells connected to a face - this function is implemented in each
  // mesh framework. The results are cached in the base class
  virtual
  void face_get_cells_internal_(const Entity_ID faceid,
          const Parallel_type ptype,
          Entity_ID_List *cellids) const override;

  // edges of a face - this function is implemented in each mesh
  // framework. The results are cached in the base class
  virtual
  void face_get_edges_and_dirs_internal_(const Entity_ID faceid,
          Entity_ID_List *edgeids,
          std::vector<int> *edge_dirs,
          const bool ordered=true) const override;

  // edges of a cell - this function is implemented in each mesh
  // framework. The results are cached in the base class.
  virtual
  void cell_get_edges_internal_(const Entity_ID cellid,
          Entity_ID_List *edgeids) const override;


 protected:
  std::vector<AmanziGeometry::Point> nodes_;
  std::map<Set_ID,bool> sets_;
  Entity_ID parent_face_;
  Cell_type cell_type_;

  Teuchos::RCP<Epetra_Map> cell_map_;
  Teuchos::RCP<Epetra_Map> face_map_;
  Teuchos::RCP<Epetra_Import> exterior_face_importer_;

};


} // close namespace AmanziMesh
} // close namespace Amanzi


#endif /* _MESH_MAPS_H_ */
