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

#include "Mesh.hh"
#include "Point.hh"
#include "GeometricModel.hh"
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


class Mesh_MSTK : public Mesh {
 public:
  // copy construction not supported
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
            const Teuchos::RCP<const AmanziGeometry::GeometricModel>& gm = Teuchos::null,
            const Teuchos::RCP<const Teuchos::ParameterList>& plist = Teuchos::null,
            const bool request_faces = true,
            const bool request_edges = false);

  // Constructors that generate a mesh internally (regular hexahedral mesh only)

  // 3D
  Mesh_MSTK(const double x0,
            const double y0,
            const double z0,
            const double x1,
            const double y1,
            const double z1,
            const unsigned int nx,
            const unsigned int ny,
            const unsigned int nz,
            const Comm_ptr_type& comm,
            const Teuchos::RCP<const AmanziGeometry::GeometricModel>& gm = Teuchos::null,
            const Teuchos::RCP<const Teuchos::ParameterList>& plist = Teuchos::null,
            const bool request_faces = true,
            const bool request_edges = false);

  // 2D
  Mesh_MSTK(const double x0,
            const double y0,
            const double x1,
            const double y1,
            const int nx,
            const int ny,
            const Comm_ptr_type& comm_,
            const Teuchos::RCP<const AmanziGeometry::GeometricModel>& gm = Teuchos::null,
            const Teuchos::RCP<const Teuchos::ParameterList>& plist = Teuchos::null,
            const bool request_faces = true,
            const bool request_edges = false);

  // Construct a mesh by extracting a subset of entities from another
  // mesh. The subset may be specified by a setname or a list of
  // entities. In some cases like extracting a surface mesh from a
  // volume mesh, constructor can be asked to "flatten" the mesh to a
  // lower dimensional space.
  Mesh_MSTK(const Teuchos::RCP<const Mesh>& parent_mesh,
            const Entity_ID_List& entity_ids,
            const Entity_kind entity_kind,
            const bool flatten = false,
            const Comm_ptr_type& comm = Teuchos::null,
            const Teuchos::RCP<const AmanziGeometry::GeometricModel>& gm = Teuchos::null,
            const Teuchos::RCP<const Teuchos::ParameterList>& plist = Teuchos::null,
            const bool request_faces = true,
            const bool request_edges = false);

  ~Mesh_MSTK();


  // Get parallel type of entity
  Parallel_type entity_get_ptype(const Entity_kind kind, const Entity_ID entid) const;

  // Get cell type
  Cell_type cell_get_type(const Entity_ID cellid) const;

  // Parent entity in the source mesh if mesh was derived from another mesh
  Entity_ID entity_get_parent(const Entity_kind kind, const Entity_ID entid) const;

  virtual Teuchos::RCP<const Mesh> parent() const { return parent_mesh_; }


  // -------------------------
  // General mesh information
  // -------------------------

  // Number of entities of any kind (cell, face, node) and in a
  // particular category (OWNED, GHOST, ALL)
  unsigned int num_entities(const Entity_kind kind, const Parallel_type ptype) const;

  // Were optional edges initialized?
  virtual bool valid_edges() const { return edges_initialized; }


  // Global ID of any entity
  Entity_ID GID(const Entity_ID lid, const Entity_kind kind) const;


  //
  // Mesh Entity Adjacencies
  //-------------------------

  // Get nodes of cell
  // On a distributed mesh, all nodes (OWNED or GHOST) of the cell
  // are returned
  // Nodes are returned in a standard order (Exodus II convention)
  // STANDARD CONVENTION WORKS ONLY FOR STANDARD CELL TYPES in 3D
  // For a general polyhedron this will return the nodes in
  // arbitrary order
  // In 2D, the nodes of the polygon will be returned in ccw order
  // consistent with the face normal
  void cell_get_nodes(const Entity_ID cellid, Entity_ID_List* nodeids) const;


  // Get nodes of face
  // On a distributed mesh, all nodes (OWNED or GHOST) of the face
  // are returned
  // In 3D, the nodes of the face are returned in ccw order consistent
  // with the face normal
  // In 2D, nfnodes is 2
  void face_get_nodes(const Entity_ID faceid, Entity_ID_List* nodeids) const;


  // Get nodes of edge On a distributed mesh all nodes (OWNED or
  // GHOST) of the face are returned
  void edge_get_nodes(const Entity_ID edgeid, Entity_ID* point0, Entity_ID* point1) const;


  //
  // Upward adjacencies
  //-------------------

  // Cells of type 'ptype' connected to a node
  void
  node_get_cells(const Entity_ID nodeid, const Parallel_type ptype, Entity_ID_List* cellids) const;

  // Faces of type 'ptype' connected to a node
  void
  node_get_faces(const Entity_ID nodeid, const Parallel_type ptype, Entity_ID_List* faceids) const;

  // Edges of type 'ptype' connected to a node
  void
  node_get_edges(const Entity_ID nodeid, const Parallel_type ptype, Entity_ID_List* edgeids) const;

  // Faces of type 'ptype' connected to an edge
  void
  edge_get_faces(const Entity_ID edgeid, const Parallel_type ptype, Entity_ID_List* faceids) const;

  // Cells of type 'ptype' connected to an edge
  void
  edge_get_cells(const Entity_ID edgeid, const Parallel_type ptype, Entity_ID_List* cellids) const;


  //
  // Same level adjacencies
  //-----------------------

  // Face connected neighboring cells of given cell of a particular ptype
  // (e.g. a hex has 6 face neighbors)

  // The order in which the cellids are returned cannot be
  // guaranteed in general except when ptype = USED, in which case
  // the cellids will correcpond to cells across the respective
  // faces given by cell_get_faces
  void cell_get_face_adj_cells(const Entity_ID cellid,
                               const Parallel_type ptype,
                               Entity_ID_List* fadj_cellids) const;


  // Node connected neighboring cells of given cell
  // (a hex in a structured mesh has 26 node connected neighbors)
  // The cells are returned in no particular order
  void cell_get_node_adj_cells(const Entity_ID cellid,
                               const Parallel_type ptype,
                               Entity_ID_List* nadj_cellids) const;


  //
  // Mesh entity geometry
  //---------------------

  // Node coordinates - 3 in 3D and 2 in 2D
  void node_get_coordinates(const Entity_ID nodeid, AmanziGeometry::Point* ncoord) const;


  // Modify the coordinates of a node
  void node_set_coordinates(const Entity_ID nodeid, const AmanziGeometry::Point coords);

  void node_set_coordinates(const Entity_ID nodeid, const double* coords);


  //
  // Epetra maps
  //------------

  const Epetra_Map& cell_map(bool include_ghost) const;

  const Epetra_Map& face_map(bool include_ghost) const;

  const Epetra_Map& edge_map(bool include_ghost) const;

  const Epetra_Map& node_map(bool include_ghost) const;

  const Epetra_Map& exterior_face_map(bool include_ghost) const;

  const Epetra_Import& exterior_face_importer(void) const;

  const Epetra_Map& exterior_node_map(bool include_ghost) const;


  //
  // Boundary Conditions or Sets
  //----------------------------
  virtual bool valid_set_type(const AmanziGeometry::RegionType rtype, const Entity_kind kind) const
  {
    return true; // MSTK supports all region types
  }

  // Get list of entities of type 'category' in set
  using Mesh::get_set_entities;

  virtual void get_set_entities_and_vofs(const std::string& setname,
                                         const Entity_kind kind,
                                         const Parallel_type ptype,
                                         std::vector<Entity_ID>* entids,
                                         std::vector<double>* vofs) const;


  // this pulls back the node-based deform as well, so that it can be called when referencing a Mesh_MSTK object
  using Mesh::deform;

  // Deform a mesh so that cell volumes conform as closely as possible
  // to target volumes without dropping below the minimum volumes.  If
  // move_vertical = true, nodes will be allowed to move only in the
  // vertical direction (right now arbitrary node movement is not allowed)
  //
  // The two versions are with and without specification of the constants
  // The one without specification of constants is present only to satisfy
  // the top level Mesh function definition for now
  //
  // the deformation function is written as:
  // log( min_vol_const1 * exp (min_vol_const2 * min_vol_func) +
  //      target_vol_const1 * exp (target_vol_const2 * target_vol_func) +
  //      quality_func_const1 * exp (quality_func_const2 * quality_func) )
  //
  // *Loosely* speaking, change the '1' constants to influence the weighting
  // of different criteria and the '2' constants to influence how tightly
  // the criteria are adhered
  int deform(const std::vector<double>& target_cell_volumes_in,
             const std::vector<double>& min_cell_volumes_in,
             const Entity_ID_List& fixed_nodes,
             const bool move_vertical)
  {
    const double min_vol_const1 = 1.0e+0;
    const double min_vol_const2 = 1.0e+2;
    const double target_vol_const1 = 1.0e+0;
    const double target_vol_const2 = 1.0e+2;
    const double quality_func_const1 = 0.0e+0; // ignore quality
    const double quality_func_const2 = 1.0e+0;

    int ierr = deform(target_cell_volumes_in,
                      min_cell_volumes_in,
                      fixed_nodes,
                      move_vertical,
                      min_vol_const1,
                      min_vol_const2,
                      target_vol_const1,
                      target_vol_const2,
                      quality_func_const1,
                      quality_func_const2);
    return ierr;
  }


  int deform(const std::vector<double>& target_cell_volumes_in,
             const std::vector<double>& min_cell_volumes_in,
             const Entity_ID_List& fixed_nodes,
             const bool move_vertical,
             const double min_vol_const1,
             const double min_vol_const2,
             const double target_vol_const1,
             const double target_vol_const2,
             const double quality_func_const1,
             const double quality_func_const2);

  // Miscellaneous
  void write_to_exodus_file(const std::string filename) const;

  // Run MSTK's internal checks - meant for debugging only
  // Returns true if everything is ok, false otherwise

  bool run_internal_mstk_checks() const;

 private:
  MPI_Comm mpicomm_;
  int myprocid, numprocs;

  Mesh_ptr mesh_;

  int serial_run;
  bool contiguous_gids_;

  // Local handles to entity lists (Vertices, "Faces", "Cells")

  // For a surface mesh, "Faces" refers to mesh edges and "Cells"
  // refers to mesh faces
  //
  // For a solid mesh, "Faces" refers to mesh faces and "Cells"
  // refers to mesh regions


  // These are MSTK's definitions of types of parallel mesh entities
  // These definitions are slightly different from what Amanzi has defined
  //
  // There are 2 types of entities relevant to this code - Owned and Ghost
  //
  // 1. OWNED - owned by this processor
  //
  // 2. GHOST - not owned by this processor

  // ALL = OWNED + GHOST in parallel, ALL = OWNED in serial
  MSet_ptr OwnedVerts, NotOwnedVerts;

  mutable MSet_ptr OwnedEdges, NotOwnedEdges;

  mutable MSet_ptr OwnedFaces, NotOwnedFaces;

  MSet_ptr OwnedCells, GhostCells;

  // Flags to indicate if face and edge info is initialized
  mutable bool faces_initialized, edges_initialized;

  // Marker to indicate if an entity is not owned
  int notwownedmark;

  // Deleted entity lists if some pre-processing had to be done
  // to the mesh to eliminate degenerate entities
  bool entities_deleted;
  List_ptr deleted_vertices, deleted_edges, deleted_faces, deleted_regions;

  // Local ID to MSTK handle map
  std::vector<MEntity_ptr> vtx_id_to_handle;
  mutable std::vector<MEntity_ptr> edge_id_to_handle;
  mutable std::vector<MEntity_ptr> face_id_to_handle;
  std::vector<MEntity_ptr> cell_id_to_handle;


  // Maps
  Epetra_Map *cell_map_w_ghosts_, *cell_map_wo_ghosts_;
  mutable Epetra_Map *face_map_w_ghosts_, *face_map_wo_ghosts_;
  mutable Epetra_Map *edge_map_w_ghosts_, *edge_map_wo_ghosts_;
  Epetra_Map *node_map_w_ghosts_, *node_map_wo_ghosts_;
  Epetra_Map *extface_map_w_ghosts_,
    *extface_map_wo_ghosts_; // exterior faces (connected to only 1 cell)
  Epetra_Map *extnode_map_w_ghosts_, *extnode_map_wo_ghosts_; // exterior node

  // Epetra importer that will allow apps to import values from a Epetra
  // vector defined on all owned faces into an Epetra vector defined
  // only on exterior faces
  Epetra_Import* owned_to_extface_importer_;

  // flag whether to flip a face dir or not when returning nodes of a
  // face (relevant only on partition boundaries)
  mutable bool* faceflip;

  // flag whether to flip an edge dir or not when returning nodes of an edge
  // (relevant only on partition boundaries)
  mutable bool* edgeflip;

  // Attribute to precompute and store celltype
  MAttrib_ptr celltype_att;

  // Parent entity attribute - populated if the mesh is derived from
  // another mesh
  MAttrib_ptr rparentatt, fparentatt, eparentatt, vparentatt;

  Teuchos::RCP<const Mesh_MSTK> parent_mesh_;

  // variables needed for mesh deformation
  double* meshxyz;
  double *target_cell_volumes_, *min_cell_volumes_, *target_weights;

  // Private methods
  // ---------------

  void clear_internals_();

  void pre_create_steps_(const int space_dimension);
  void post_create_steps_(const bool request_faces, const bool request_edges);

  void init_mesh_from_file_(const std::string& filename,
                            const Partitioner_type partitioner = PARTITIONER_DEFAULT);

  void collapse_degen_edges();
  Cell_type MFace_Celltype(MFace_ptr f);
  Cell_type MRegion_Celltype(MRegion_ptr r);
  void label_celltype();

  void init_pvert_lists();
  void init_pedge_lists();
  void init_pedge_dirs();
  void init_pface_lists();
  void init_pface_dirs();
  void init_pface_dirs_3();
  void init_pface_dirs_2();
  void init_pcell_lists();

  void init_vertex_id2handle_maps();
  void init_edge_id2handle_maps();
  void init_face_id2handle_maps();
  void init_cell_id2handle_maps();
  void init_global_ids();

  void init_cell_map();
  void init_face_map();
  void init_edge_map();
  void init_node_map();

  void init_nodes();
  void init_edges();
  void init_faces();
  void init_cells();

  void init_set_info();
  void inherit_labeled_sets(MAttrib_ptr copyatt, List_ptr src_entities);
  std::string internal_name_of_set(const Teuchos::RCP<const AmanziGeometry::Region>& region,
                                   const Entity_kind entity_kind) const;
  std::string other_internal_name_of_set(const Teuchos::RCP<const AmanziGeometry::Region>& r,
                                         const Entity_kind entity_kind) const;

  int generate_regular_mesh(Mesh_ptr mesh,
                            double x0,
                            double y0,
                            double z0,
                            double x1,
                            double y1,
                            double z1,
                            int nx,
                            int ny,
                            int nz);
  int
  generate_regular_mesh(Mesh_ptr mesh, double x0, double y0, double x1, double y1, int nx, int ny);

  void extract_mstk_mesh(const List_ptr entity_ids,
                         const MType entity_dim,
                         const bool flatten = false,
                         const bool request_faces = true,
                         const bool request_edges = false);

  MSet_ptr
  build_set(const Teuchos::RCP<const AmanziGeometry::Region>& region, const Entity_kind kind) const;

  bool is_boundary_node_(const MEntity_ptr ment) const;

  //
  // Downward Adjacencies
  //---------------------

  // Get faces of a cell and directions in which the cell uses the face

  // The Amanzi coding guidelines regarding function arguments is purposely
  // violated here to allow for a default input argument

  // On a distributed mesh, this will return all the faces of the
  // cell, OWNED or GHOST. If ordered = true, the faces will be
  // returned in a standard order according to Exodus II convention
  // for standard cells; in all other situations (ordered = false or
  // non-standard cells), the list of faces will be in arbitrary order

  // In 3D, direction is 1 if face normal points out of cell
  // and -1 if face normal points into cell
  // In 2D, direction is 1 if face/edge is defined in the same
  // direction as the cell polygon, and -1 otherwise

  void cell_get_faces_and_dirs_internal_(const Entity_ID cellid,
                                         Entity_ID_List* faceids,
                                         std::vector<int>* face_dirs,
                                         const bool ordered = false) const;

  void cell_get_faces_and_dirs_ordered(const Entity_ID cellid,
                                       Entity_ID_List* faceids,
                                       std::vector<int>* face_dirs) const;

  void cell_get_faces_and_dirs_unordered(const Entity_ID cellid,
                                         Entity_ID_List* faceids,
                                         std::vector<int>* face_dirs) const;


  // Cells connected to a face
  void face_get_cells_internal_(const Entity_ID faceid,
                                const Parallel_type ptype,
                                Entity_ID_List* cellids) const;


  // Get edges of a cell
  void cell_get_edges_internal_(const Entity_ID cellid, Entity_ID_List* edgeids) const;

  // Edges and edge directions of a face
  void face_get_edges_and_dirs_internal_(const Entity_ID cellid,
                                         Entity_ID_List* edgeids,
                                         std::vector<int>* edgedirs,
                                         bool ordered = true) const;

  // Map from Amanzi's mesh entity kind to MSTK's mesh type.
  MType entity_kind_to_mtype(const Entity_kind kind) const
  {
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

    static MType const kind2mtype[4][4] = { { MVERTEX, MVERTEX, MVERTEX, MVERTEX }, // 0d meshes
                                            { MVERTEX, MVERTEX, MVERTEX, MEDGE },   // 1d meshes
                                            { MVERTEX, MEDGE, MEDGE, MFACE },       // 2d meshes
                                            { MVERTEX, MEDGE, MFACE, MREGION } };   // 3d meshes

    return kind2mtype[manifold_dimension()][(int)kind];
  }


  // Compute the value of the LOCAL component of the GLOBAL
  // deformation objective function given a new position 'nodexyz' for
  // node 'nodeid' i.e. only those terms in the global function that
  // are affected by the movement of this node.
  double deform_function(const int nodeid, double const* const nodexyz) const;

  // Finite difference gradient of deformation objective function
  void deform_gradient(const int nodeid, double const* const vxyz, double* gradient) const;

  // Finite difference hessian of deformation objective function
  void deform_hessian(const int nodeid, double const* const nodexyz, double hessian[3][3]) const;

  // Minimum eigen value of a matrix (rank 2 and rank 3)
  double mineigenvalue(const double A[3][3]) const;

  // Inverse of Hessian of rank 2 or 3
  int hessian_inverse(const double H[3][3], double iH[3][3]) const;
};


//------------
// Epetra maps
//------------

inline const Epetra_Map&
Mesh_MSTK::cell_map(bool include_ghost) const
{
  if (serial_run)
    return *cell_map_wo_ghosts_;
  else
    return (include_ghost ? *cell_map_w_ghosts_ : *cell_map_wo_ghosts_);
}


inline const Epetra_Map&
Mesh_MSTK::face_map(bool include_ghost) const
{
  if (serial_run)
    return *face_map_wo_ghosts_;
  else
    return (include_ghost ? *face_map_w_ghosts_ : *face_map_wo_ghosts_);
}


inline const Epetra_Map&
Mesh_MSTK::edge_map(bool include_ghost) const
{
  if (serial_run)
    return *edge_map_wo_ghosts_;
  else
    return (include_ghost ? *edge_map_w_ghosts_ : *edge_map_wo_ghosts_);
}


inline const Epetra_Map&
Mesh_MSTK::node_map(bool include_ghost) const
{
  if (serial_run)
    return *node_map_wo_ghosts_;
  else
    return (include_ghost ? *node_map_w_ghosts_ : *node_map_wo_ghosts_);
}


inline const Epetra_Map&
Mesh_MSTK::exterior_face_map(bool include_ghost) const
{
  if (serial_run)
    return *extface_map_wo_ghosts_;
  else
    return (include_ghost ? *extface_map_w_ghosts_ : *extface_map_wo_ghosts_);
}


inline const Epetra_Import&
Mesh_MSTK::exterior_face_importer(void) const
{
  return *owned_to_extface_importer_;
}


inline const Epetra_Map&
Mesh_MSTK::exterior_node_map(bool include_ghost) const
{
  if (serial_run)
    return *extnode_map_wo_ghosts_;
  else
    return (include_ghost ? *extnode_map_w_ghosts_ : *extnode_map_wo_ghosts_);
}


inline Parallel_type
Mesh_MSTK::entity_get_ptype(const Entity_kind kind, const Entity_ID entid) const
{
  MEntity_ptr ment;

  switch (kind) {
  case CELL:
    ment = (MEntity_ptr)cell_id_to_handle[entid];
    break;
  case FACE:
    ment = (MEntity_ptr)face_id_to_handle[entid];
    break;
  case NODE:
    ment = (MEntity_ptr)vtx_id_to_handle[entid];
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
