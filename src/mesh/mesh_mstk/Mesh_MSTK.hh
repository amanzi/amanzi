/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Rao Garimella, others
*/

//! Implementation of the Mesh interface leveraging MSTK.


#ifndef AMANZI_MESH_MSTK_HH_
#define AMANZI_MESH_MSTK_HH_

#include <memory>
#include <vector>
#include <sstream>

#include "AmanziComm.hh"
#include "MSTK.h"

#include "MeshFramework.hh"
#include "Point.hh"
#include "GenerationSpec.hh"
#include "dbc.hh"
#include "errors.hh"

namespace Amanzi {

namespace AmanziMesh {

// Mesh class based on the MSTK framework.
//
// Instantiating a const version of this class only guarantees that
// the underlying mesh topology and geometry does not change. For
// purposes of memory savings we use lazy initialization of face and
// edge lists (they are already present in the underlying MSTK mesh
// data structures), which means that the data structures holding the
// mesh information may still change.


class Mesh_MSTK : public MeshFramework {
 public:
  Mesh_MSTK(const Mesh_MSTK& other) = delete;

  // Constructors that read the mesh from a file
  // The request_faces and request_edges arguments have to be at the
  // end and not in the middle because if we omit them and specify a
  // pointer argument like gm or verbosity_obj, then there is implicit
  // conversion of the pointer to bool, thereby defeating the intent
  // of the call and making the pointer argument seem NULL. In C++11,
  // we could "delete" the illegal version of the call effectively
  // blocking the implicit conversion.
  Mesh_MSTK(const std::string& filename,
            const Comm_ptr_type& comm,
            const Teuchos::RCP<const AmanziGeometry::GeometricModel>& gm=Teuchos::null,
            const Teuchos::RCP<Teuchos::ParameterList>& plist=Teuchos::null);

  // Constructors that generate a mesh internally (regular hexahedral mesh only)

  // 3D
  Mesh_MSTK(const double x0, const double y0, const double z0,
            const double x1, const double y1, const double z1,
            const unsigned int nx, const unsigned int ny,
            const unsigned int nz,
            const Comm_ptr_type& comm,
            const Teuchos::RCP<const AmanziGeometry::GeometricModel>& gm=Teuchos::null,
            const Teuchos::RCP<Teuchos::ParameterList>& plist=Teuchos::null);

  // 2D
  Mesh_MSTK(const double x0, const double y0,
            const double x1, const double y1,
            const int nx, const int ny,
            const Comm_ptr_type& comm_,
            const Teuchos::RCP<const AmanziGeometry::GeometricModel>& gm=Teuchos::null,
            const Teuchos::RCP<Teuchos::ParameterList>& plist=Teuchos::null);

  // Construct a mesh by extracting a subset of entities from another
  // mesh. The subset may be specified by a setname or a list of
  // entities. In some cases like extracting a surface mesh from a
  // volume mesh, constructor can be asked to "flatten" the mesh to a
  // lower dimensional space.
  Mesh_MSTK(const Teuchos::RCP<const MeshFramework>& parent_mesh,
            const Entity_ID_List& entity_ids,
            const Entity_kind entity_kind,
            const bool flatten=false,
            const Comm_ptr_type& comm=Teuchos::null,
            const Teuchos::RCP<const AmanziGeometry::GeometricModel>& gm=Teuchos::null,
            const Teuchos::RCP<Teuchos::ParameterList>& plist=Teuchos::null);

  ~Mesh_MSTK();

  virtual void write_to_exodus_file(const std::string& filename) const override;

  //
  // Accessors/mutators
  //
  // Parent entity in the source mesh if mesh was derived from another mesh
  virtual Teuchos::RCP<const MeshFramework> get_parent() const override { return parent_mesh_; }
  void set_parent(const Teuchos::RCP<const Mesh_MSTK>& parent) { parent_mesh_ = parent; }

  virtual bool has_edges() const override { return edges_requested_; }
  virtual bool is_deformable() const override { return true; }

  //
  // Entity meta-data
  //
  // Number of entities of any kind (cell, face, node) and in a
  // particular category (OWNED, GHOST, ALL)
  virtual std::size_t getNumEntities(const Entity_kind kind, const Parallel_type ptype) const override;
  // Get parallel type of entity
  virtual Parallel_type getEntityPtype(const Entity_kind kind, const Entity_ID entid) const override;
  // Global ID of any entity
  virtual Entity_GID getEntityGID(const Entity_kind kind, const Entity_ID lid) const override;

  // corresponding entity in the parent mesh
  virtual Entity_ID getEntityParent(const Entity_kind kind, const Entity_ID entid) const override;

  // Get cell type
  Cell_type getCellType(const Entity_ID cellid) const override;

  //---------------------
  // Geometry
  //---------------------
  // Node coordinates - 3 in 3D and 2 in 2D
  virtual AmanziGeometry::Point getNodeCoordinate(const Entity_ID nodeid) const override;

  // Modify the coordinates of a node
  virtual void setNodeCoordinate(const Entity_ID nodeid,
          const AmanziGeometry::Point& coord) override;

  //
  // Downward Adjacencies
  //-------------------------
  //
  // MSTK overrides all of these.  It is a bit unclear whether that is good or
  // bad.  It may be good, because it allows a specific order to be enforced
  // (that of the standard Exodus notation).  It may be bad because there are
  // opportunities to have inconsistent ordering.

  // Get faces of a cell and directions in which the cell uses the face
  //
  // On a distributed mesh, this will return all the faces of the
  // cell, OWNED or GHOST. If ordered = true, the faces will be
  // returned in a standard order according to Exodus II convention
  // for standard cells; in all other situations (ordered = false or
  // non-standard cells), the list of faces will be in arbitrary order
  //
  // In 3D, direction is 1 if face normal points out of cell
  // and -1 if face normal points into cell
  // In 2D, direction is 1 if face/edge is defined in the same
  // direction as the cell polygon, and -1 otherwise
  virtual void getCellFacesAndDirs(const Entity_ID cellid,
          Entity_ID_List& faceids,
          Entity_Direction_List * const face_dirs) const override;

  // Get edges of a cell
  virtual void getCellEdges(const Entity_ID cellid,
                  Entity_ID_List& edgeids) const override;

  // Get edges and directions of a 2D cell
  virtual void getCell2DEdgesAndDirs(const Entity_ID cellid,
          Entity_ID_List& edgeids,
          Entity_Direction_List * const edgedirs) const override;

  // Get nodes of cell
  // On a distributed mesh, all nodes (OWNED or GHOST) of the cell
  // are returned
  // Nodes are returned in a standard order (Exodus II convention)
  // STANDARD CONVENTION WORKS ONLY FOR STANDARD CELL TYPES in 3D
  // For a general polyhedron this will return the nodes in
  // arbitrary order
  // In 2D, the nodes of the polygon will be returned in ccw order
  // consistent with the face normal
  virtual void getCellNodes(const Entity_ID cellid,
                          Entity_ID_List& nodeids) const override;


  // Edges and edge directions of a face
  virtual void getFaceEdgesAndDirs(const Entity_ID cellid,
                           Entity_ID_List& edgeids,
                           Entity_Direction_List * const edgedirs) const override;

  // Get nodes of face
  // On a distributed mesh, all nodes (OWNED or GHOST) of the face
  // are returned
  // In 3D, the nodes of the face are returned in ccw order consistent
  // with the face normal
  // In 2D, nfnodes is 2
  virtual void getFaceNodes(const Entity_ID faceid,
                          Entity_ID_List& nodeids) const override;

  // Get nodes of edge On a distributed mesh all nodes (OWNED or
  // GHOST) of the face are returned
  virtual void getEdgeNodes(const Entity_ID edgeid,
                          Entity_ID_List& nodes) const override;

  //
  // Upward adjacencies
  //-------------------
  // Cells connected to a face
  virtual void getFaceCells(const Entity_ID faceid,
                          const Parallel_type ptype,
                          Entity_ID_List& cellids) const override;

  // Cells of type 'ptype' connected to an edge
  virtual void getEdgeCells(const Entity_ID edgeid,
                          const Parallel_type ptype,
                          Entity_ID_List& cellids) const override;

  // Faces of type 'ptype' connected to an edge
  virtual void getEdgeFaces(const Entity_ID edgeid,
                          const Parallel_type ptype,
                          Entity_ID_List& faceids) const override;

  // Cells of type 'ptype' connected to a node
  virtual void getNodeCells(const Entity_ID nodeid,
                          const Parallel_type ptype,
                          Entity_ID_List& cellids) const override;

  // Faces of type 'ptype' connected to a node
  virtual void getNodeFaces(const Entity_ID nodeid,
                          const Parallel_type ptype,
                          Entity_ID_List& faceids) const override;

  // Edges of type 'ptype' connected to a node
  virtual void getNodeEdges(const Entity_ID nodeid,
                          const Parallel_type ptype,
                          Entity_ID_List& edgeids) const override;

  //
  // Get list of entities of type kind from a framework set.
  //
  // MSTK supports labeled sets
  virtual void getSetEntities(const AmanziGeometry::RegionLabeledSet& region,
          const Entity_kind kind,
          const Parallel_type ptype,
          Entity_ID_List& entids) const override;

  // Run MSTK's internal checks - meant for debugging only
  // Returns true if everything is ok, false otherwise
  bool run_internal_mstk_checks() const;

 // Private methods
 // ---------------
 private:

  void clear_internals_();
  void read_plist_();
  void pre_create_steps_(const int space_dimension);
  void post_create_steps_();

  // functions used in construction
  void init_mesh_from_file_(const std::string& filename);
  int generate_regular_mesh_(Mesh_ptr mesh, double x0, double y0, double z0,
                            double x1, double y1, double z1, int nx,
                            int ny, int nz);
  int generate_regular_mesh_(Mesh_ptr mesh, double x0, double y0,
                            double x1, double y1, int nx, int ny);
  void extract_mstk_mesh_(const List_ptr entity_ids,
                          const MType entity_dim,
                          const bool flatten = false,
                          const bool request_faces = true,
                          const bool request_edges = false);

  void collapse_degen_edges_();
  void label_celltype_();

  void init_pvert_lists_();
  void init_pedge_lists_();
  void init_pedge_dirs_();
  void init_pface_lists_();
  void init_pface_dirs_();
  void init_pface_dirs_3_();
  void init_pface_dirs_2_();
  void init_pcell_lists_();

  void init_vertex_id2handle_maps_();
  void init_edge_id2handle_maps_();
  void init_face_id2handle_maps_();
  void init_cell_id2handle_maps_();
  void init_global_ids_();

  void init_nodes_();
  void init_edges_();
  void init_faces_();
  void init_cells_();

  void inherit_labeled_sets_(MAttrib_ptr copyatt, List_ptr src_entities);
  std::string internal_name_of_set_(const AmanziGeometry::RegionLabeledSet& region,
          const Entity_kind entity_kind) const;
  std::string other_internal_name_of_set_(const AmanziGeometry::RegionLabeledSet& r,
          const Entity_kind entity_kind) const;


  // MSet_ptr build_set(const Teuchos::RCP<const AmanziGeometry::Region>& region,
  //                    const Entity_kind kind) const;

  // bool is_boundary_node_(const MEntity_ptr ment) const;


  // Map from Amanzi's mesh entity kind to MSTK's mesh type.
  MType entity_kind_to_mtype_(const Entity_kind kind) const {
    // The first index is cell dimension (0,1,2,3) and the second index
    // is the entity kind
    //
    // map order in each row is NODE, EDGE, FACE, CELL
    //
    // So, for a 1D mesh, nodes are MVERTEX type in MSTK, edges and faces
    // are also MVERTEX type, and cells are MEDGE type
    //
    // For a 2D mesh, nodes are MVERTEX type, edges and faces are MEDGE
    // type, and cells are MFACE type
    static MType const
      kind2mtype[4][4] = {{MVERTEX, MVERTEX, MVERTEX, MVERTEX},  // 0d meshes
                          {MVERTEX, MVERTEX, MVERTEX, MEDGE},    // 1d meshes
                          {MVERTEX, MEDGE,   MEDGE,   MFACE},    // 2d meshes
                          {MVERTEX, MEDGE,   MFACE,   MREGION}}; // 3d meshes
    return kind2mtype[get_manifold_dimension()][(int)kind - 1];
  }

  void getCellFacesAndDirs_ordered_(const Entity_ID cellid,
          Entity_ID_List& faceids,
          Entity_Direction_List * const face_dirs) const;

  void getCellFacesAndDirs_unordered_(const Entity_ID cellid,
          Entity_ID_List &faceids,
          Entity_Direction_List * const face_dirs) const;

 private:
  MPI_Comm mpicomm_;
  int myprocid, numprocs;
  Mesh_ptr mesh_;

  int serial_run;
  bool contiguous_gids_;
  Partitioner_type partitioner_;
  bool edges_requested_, faces_requested_;
  mutable bool edges_initialized_, faces_initialized_, cells_initialized_;

  // Local handles to entity lists (Vertices, "Faces", "Cells")
  //
  // For a surface mesh, "Faces" refers to mesh edges and "Cells"
  // refers to mesh faces
  //
  // For a solid mesh, "Faces" refers to mesh faces and "Cells"
  // refers to mesh regions
  //
  //
  // These are MSTK's definitions of types of parallel mesh entities
  // These definitions are slightly different from what Amanzi has defined
  //
  // There are 2 types of entities relevant to this code - Owned and Ghost
  //
  // 1. OWNED - owned by this processor
  //
  // 2. GHOST - not owned by this processor
  //
  // ALL = OWNED + GHOST in parallel, ALL = OWNED in serial
  MSet_ptr owned_verts_, ghost_verts_;
  mutable MSet_ptr owned_edges_, ghost_edges_;
  mutable MSet_ptr owned_faces_, ghost_faces_;
  MSet_ptr owned_cells_, ghost_cells_;

  // Deleted entity lists if some pre-processing had to be done
  // to the mesh to eliminate degenerate entities
  //
  // Likely these can be removed?  They don't seem to be used...
  bool entities_deleted_;
  List_ptr deleted_vertices_, deleted_edges_, deleted_faces_, deleted_regions_;

  // Local ID to MSTK handle map
  std::vector<MEntity_ptr> vtx_id_to_handle_;
  mutable std::vector<MEntity_ptr> edge_id_to_handle_;
  mutable std::vector<MEntity_ptr> face_id_to_handle_;
  std::vector<MEntity_ptr> cell_id_to_handle_;

  // flag whether to flip a face dir or not when returning nodes of a
  // face (relevant only on partition boundaries)
  mutable bool *faceflip_;

  // flag whether to flip an edge dir or not when returning nodes of an edge
  // (relevant only on partition boundaries)
  mutable bool *edgeflip_;

  // Attribute to precompute and store celltype
  MAttrib_ptr celltype_att_;

  // Parent entity attribute - populated if the mesh is derived from
  // another mesh
  MAttrib_ptr rparentatt_, fparentatt_, eparentatt_, vparentatt_;
  Teuchos::RCP<const Mesh_MSTK> parent_mesh_;

};


inline
Parallel_type Mesh_MSTK::getEntityPtype(const Entity_kind kind, const Entity_ID entid) const {
  MEntity_ptr ment;

  switch(kind) {
    case Entity_kind::CELL:
      ment = (MEntity_ptr) cell_id_to_handle_[entid];
      break;
    case Entity_kind::FACE:
      ment = (MEntity_ptr) face_id_to_handle_[entid];
      break;
    case Entity_kind::NODE:
      ment = (MEntity_ptr) vtx_id_to_handle_[entid];
      break;
    default:
      ment = NULL;
    }

    if (MEnt_PType(ment) == PGHOST)
      return Parallel_type::GHOST;
    else
      return Parallel_type::OWNED;
}

} // namespace AmanziMesh
} // namespace Amanzi

#endif

