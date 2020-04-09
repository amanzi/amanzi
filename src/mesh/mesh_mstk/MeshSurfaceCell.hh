/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
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
#include "AmanziMap.hh"
#include "AmanziComm.hh"

#include "VerboseObject.hh"
#include "dbc.hh"
#include "errors.hh"

#include "Region.hh"
#include "Mesh_MSTK.hh"

namespace Amanzi {
namespace AmanziMesh {

class MeshSurfaceCell : public Mesh {
 public:
  MeshSurfaceCell(const Teuchos::RCP<const Mesh>& parent_mesh,
                  const std::string& setname, bool flatten = true);

  ~MeshSurfaceCell() = default;

  // Get parallel type of entity - OWNED, GHOST, ALL (See MeshDefs.hh)
  virtual Parallel_type
  entity_get_ptype(const Entity_kind kind, const Entity_ID entid) const
  {
    return Parallel_type::OWNED;
  }

  // Parent entity in the source mesh if mesh was derived from another mesh
  virtual Entity_ID
  entity_get_parent_type(const Entity_kind kind, const Entity_ID entid) const
  {
    AMANZI_ASSERT(kind == CELL);
    AMANZI_ASSERT(entid == 0);
    return parent_face_;
  }

  // Get cell type - UNKNOWN, TRI, QUAD, ... See MeshDefs.hh
  virtual Cell_type cell_get_type(const Entity_ID cellid) const
  {
    return cell_type_;
  }


  //
  // General mesh information
  // -------------------------
  //

  // Number of entities of any kind (cell, face, node) and in a
  // particular category (OWNED, GHOST, ALL)
  virtual unsigned int
  num_entities(const Entity_kind kind, const Parallel_type ptype) const;

  // Global ID of any entity
  virtual Entity_ID GID(const Entity_ID lid, const Entity_kind kind) const
  {
    return lid;
  }


  //
  // Mesh Entity Adjacencies
  //-------------------------

  // Downward Adjacencies
  //---------------------

  // Get nodes of a cell
  virtual void cell_get_nodes(const Entity_ID cellid,
                              Entity_ID_List& nodeids) const;

  // Get nodes of face
  // On a distributed mesh, all nodes (OWNED or GHOST) of the face
  // are returned
  // In 3D, the nodes of the face are returned in ccw order consistent
  // with the face normal
  // In 2D, nfnodes is 2
  virtual void face_get_nodes(const Entity_ID faceid,
                              Entity_ID_List& nodeids) const;

  // Get nodes of edge
  virtual void edge_get_nodes(const Entity_ID edgeid, Entity_ID* nodeid0,
                              Entity_ID* nodeid1) const;

  // Upward adjacencies
  //-------------------

  // Cells of type 'ptype' connected to a node - The order of cells is
  // not guaranteed to be the same for corresponding nodes on
  // different processors
  virtual void node_get_cells(const Entity_ID nodeid, const Parallel_type ptype,
                              Entity_ID_List& cellids) const;

  // Faces of type 'ptype' connected to a node - The order of faces is
  // not guarnateed to be the same for corresponding nodes on
  // different processors
  virtual void node_get_faces(const Entity_ID nodeid, const Parallel_type ptype,
                              Entity_ID_List& faceids) const;

  // Get faces of ptype of a particular cell that are connected to the
  // given node - The order of faces is not guarnateed to be the same
  // for corresponding nodes on different processors
  virtual void
  node_get_cell_faces(const Entity_ID nodeid, const Entity_ID cellid,
                      const Parallel_type ptype,
                      Entity_ID_List& faceids) const;
  // Cells of type 'ptype' connected to an edges
  virtual void edge_get_cells(const Entity_ID edgeid, const Parallel_type ptype,
                              Entity_ID_List& cellids) const;

  // Same level adjacencies
  //-----------------------

  // Face connected neighboring cells of given cell of a particular ptype
  // (e.g. a hex has 6 face neighbors)

  // The order in which the cellids are returned cannot be
  // guaranteed in general except when ptype = ALL, in which case
  // the cellids will correcpond to cells across the respective
  // faces given by cell_get_faces
  virtual void
  cell_get_face_adj_cells(const Entity_ID cellid, const Parallel_type ptype,
                          Entity_ID_List& fadj_cellids) const;

  // Node connected neighboring cells of given cell
  // (a hex in a structured mesh has 26 node connected neighbors)
  // The cells are returned in no particular order
  virtual void
  cell_get_node_adj_cells(const Entity_ID cellid, const Parallel_type ptype,
                          Entity_ID_List& nadj_cellids) const;


  //
  // Mesh entity geometry
  //--------------
  //

  // Node coordinates - 3 in 3D and 2 in 2D
  virtual void node_get_coordinates(const Entity_ID nodeid,
                                    AmanziGeometry::Point* ncoord) const
  {
    (*ncoord) = nodes_[nodeid];
  }


  // Face coordinates - conventions same as face_to_nodes call
  // Number of nodes is the vector size divided by number of spatial dimensions
  virtual void
  face_get_coordinates(const Entity_ID faceid,
                       std::vector<AmanziGeometry::Point>& fcoords) const
  {
    fcoords.resize(2); 
    fcoords[0] = nodes_[faceid];
    fcoords[1] = nodes_[(faceid + 1) % nodes_.size()];
  }

  // Coordinates of cells in standard order (Exodus II convention)
  // STANDARD CONVENTION WORKS ONLY FOR STANDARD CELL TYPES IN 3D
  // For a general polyhedron this will return the node coordinates in
  // arbitrary order
  // Number of nodes is vector size divided by number of spatial dimensions
  virtual void
  cell_get_coordinates(const Entity_ID cellid,
                       std::vector<AmanziGeometry::Point>& ccoords) const
  {
    ccoords.resize(nodes_.size());
    for (int i = 0; i < nodes_.size(); ++i) ccoords[i] = nodes_[i];
  }


  //
  // Mesh modification
  //-------------------

  // Set coordinates of node
  virtual void node_set_coordinates(const Entity_ID nodeid,
                                    const AmanziGeometry::Point ncoord)
  {
    nodes_[nodeid] = ncoord;
  }


  virtual void
  node_set_coordinates(const Entity_ID nodeid, const double* ncoord)
  {
    Errors::Message mesg("Not implemented");
    Exceptions::amanzi_throw(mesg);
  }


  // Deform a mesh so that cell volumes conform as closely as possible
  // to target volumes without dropping below the minimum volumes.  If
  // move_vertical = true, nodes will be allowed to move only in the
  // vertical direction (right now arbitrary node movement is not allowed)
  // Nodes in any set in the fixed_sets will not be permitted to move.
  virtual int
  deform(const std::vector<double>& target_cell_volumes_in,
         const std::vector<double>& min_cell_volumes_in,
         const Entity_ID_List& fixed_nodes, const bool move_vertical);
  //
  // Epetra maps
  //------------
  virtual Map_ptr_type cell_map(bool include_ghost) const { return cell_map_; }

  virtual Map_ptr_type face_map(bool include_ghost) const { return face_map_; }

  // dummy implementation so that frameworks can skip or overwrite
  Map_ptr_type edge_map(bool include_ghost) const
  {
    Errors::Message mesg("Edges not implemented in this framework");
    Exceptions::amanzi_throw(mesg);
    throw(mesg);
  };

  virtual Map_ptr_type node_map(bool include_ghost) const { return face_map_; }

  virtual Map_ptr_type exterior_face_map(bool include_ghost) const
  {
    return face_map_;
  }


  // Epetra importer that will allow apps to import values from a
  // Epetra vector defined on all owned faces into an Epetra vector
  // defined only on exterior faces
  virtual Import_ptr_type exterior_face_importer(void) const
  {
    return exterior_face_importer_;
  }


  //
  // Mesh Sets for ICs, BCs, Material Properties and whatever else
  //--------------------------------------------------------------
  //

  // Get number of entities of type 'category' in set
  virtual unsigned int get_set_size(const Set_ID setid, const Entity_kind kind,
                                    const Parallel_type ptype) const;

  virtual unsigned int
  get_set_size(const std::string setname, const Entity_kind kind,
               const Parallel_type ptype) const;

  // Get list of entities of type 'category' in set
  virtual void get_set_entities(const Set_ID setid, const Entity_kind kind,
                                const Parallel_type ptype,
                                Entity_ID_List& entids) const;

  virtual void
  get_set_entities_and_vofs(const std::string setname, const Entity_kind kind,
                            const Parallel_type ptype,
                            Entity_ID_List& entids,
                            std::vector<double>* vofs) const;


  // Miscellaneous functions
  virtual void write_to_exodus_file(const std::string filename) const;

 protected:
  // get faces and face dirs of a cell. This can be called by
  // cell_get_faces_and_dirs method of the base class and the data
  // cached or it can be called directly by the
  // cell_get_faces_and_dirs method of this class
  virtual void
  cell_get_faces_and_dirs_internal_(const Entity_ID cellid,
                                    Entity_ID_List& faceids,
                                    std::vector<int>& face_dirs) const;

  // Cells connected to a face - this function is implemented in each
  // mesh framework. The results are cached in the base class
  virtual void face_get_cells_internal_(
    const Entity_ID faceid, const Parallel_type ptype,
    Entity_ID_List& cellids) const;

  // edges of a face - this function is implemented in each mesh
  // framework. The results are cached in the base class
  virtual void
  face_get_edges_and_dirs_internal_(const Entity_ID faceid,
                                    Entity_ID_List& edgeids,
                                    std::vector<int>* edge_dirs,
                                    const bool ordered = true) const;

  // edges of a cell - this function is implemented in each mesh
  // framework. The results are cached in the base class.
  virtual void
  cell_get_edges_internal_(const Entity_ID cellid,
                           Entity_ID_List& edgeids) const;


  // edges and directions of a 2D cell - this function is implemented
  // in each mesh framework. The results are cached in the base class.
  virtual void
  cell_2D_get_edges_and_dirs_internal_(const Entity_ID cellid,
                                       Entity_ID_List& edgeids,
                                       std::vector<int>* edge_dirs) const;

 protected:
  Teuchos::RCP<const Mesh> parent_mesh_;

  std::vector<AmanziGeometry::Point> nodes_;
  std::map<Set_ID, bool> sets_;
  Entity_ID parent_face_;
  Cell_type cell_type_;

  Map_ptr_type cell_map_;
  Map_ptr_type face_map_;
  Import_ptr_type exterior_face_importer_;
};


} // namespace AmanziMesh
} // namespace Amanzi


#endif /* _MESH_MAPS_H_ */
