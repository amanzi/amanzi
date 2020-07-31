/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Rao Garimella, others
*/

//! Interface for a Mesh object.

/*!

Use the associated mesh factory to create an instance of a
derived class based on a particular mesh framework (like MSTK,
STKmesh etc.)

Assumptions:

Cells in a 2D mesh must be oriented counter clockwise. This is not
applied to meshes on manifolds.   

Design documentation:

This class is designed to be both flexible and performant (and somewhat
succeeds at both).  Many of the low-level methods here are accessed in the
innermost loop of performance-critical physics methods, so must themselves
be as simple as possible.  To meet this need, the low-level geometry ( such
as cell_volume() ) and relational topology ( such as cell_get_faces() ) are
ideally inlined to a single array access.

To accomplish this goal, this class stores a cache of data which is
accessed using the public interface.  This cache is (lazily) built up using
"cache-building" methods (i.e. compute_cell_geometry_() ) which are virtual,
but have a default implementation here.  These default implementations
simply call more virtual methods (i.e. cell_get_faces_internal() ) which
are NOT implemented here and instead use the various mesh frameworks.


There are a few exceptions to this -- Mesh_Simple and classes that inherit
from Mesh_Simple have no internal storage -- they use the cache directly as
their mesh representation.  They are not full-featured, but are useful for
some simple structured tests and mesh manipulations where a mesh is implied
(i.e. MeshLogical, MeshColumn, etc).

Note that not everything is cached, which would be memory overkill.  Only
things which have proven (or more accurately, a few were proven and the
rest were assumed based on those results) to be an issue have been cached.

NOTE: Lazy definition of the cache itself is necessarily "mutable".

*/


#ifndef AMANZI_MESH_HH_
#define AMANZI_MESH_HH_

#include <map>
#include <vector>
#include <string>

#include "Epetra_Map.h"
#include "Epetra_Import.h"
#include "nanoflann.hpp"

#include "errors.hh"
#include "AmanziComm.hh"
#include "VerboseObject.hh"
#include "Key.hh"

#include "GeometricModel.hh"
#include "Point.hh"
#include "Region.hh"

#include "CellTopology.hh"
#include "KDTree.hh"
#include "MeshDefs.hh"

// set to 0 to avoid using cache for profiling or debugging
#define AMANZI_MESH_CACHE_VARS 1

namespace Amanzi {
namespace AmanziMesh {

// friended class that uses cache
class MeshEmbeddedLogical;

class Mesh {
 protected:
  // cannot create a Mesh without type
  Mesh(const Comm_ptr_type& comm,
       const Teuchos::RCP<const AmanziGeometry::GeometricModel>& gm,
       const Teuchos::RCP<const Teuchos::ParameterList>& plist,
       const bool request_faces,
       const bool request_edges);

 public:
  
  // virtual destructor
  virtual ~Mesh() = default;

  //
  // Accessors and Mutators
  // ----------------------
  Comm_ptr_type get_comm() const {
    return comm_;
  }

  // accessor for verbose object
  Teuchos::RCP<const VerboseObject> verbosity_obj() const {
    return vo_;
  }

  // reference for vis.  Usually this is *this, but not necessarily in derived meshes
  virtual
  const Mesh& vis_mesh() const {
    return *this;
  }

  // reference for a parent mesh.  Usually this is null, but derived meshes may provide it.
  virtual
  Teuchos::RCP<const Mesh> parent() const {
    return parent_;
  }

  // Set/get the space dimension
  void set_space_dimension(const unsigned int dim) {
    space_dim_ = dim;
  }
  unsigned int space_dimension() const {
    return space_dim_;
  }

  // Set/get the manifold dimension
  void set_manifold_dimension(const unsigned int dim) {
    manifold_dim_ = dim;   // 3 is solid mesh, 2 is surface mesh
  }
  unsigned int manifold_dimension() const {
    return manifold_dim_;
  }

  // Set/Get geometric model
  void set_geometric_model(const Teuchos::RCP<const AmanziGeometry::GeometricModel>& gm) {
    geometric_model_ = gm;
  }
  Teuchos::RCP<const AmanziGeometry::GeometricModel> geometric_model() const {
    return geometric_model_;
  }

  // Get parameter list
  void set_parameter_list(const Teuchos::RCP<const Teuchos::ParameterList>& plist) {
    plist_ = plist;
    if (vo_ == Teuchos::null)
      vo_ = Teuchos::rcp(new VerboseObject(comm_,Keys::cleanPListName(plist->name()), *plist));
  }
  Teuchos::RCP<const Teuchos::ParameterList> parameter_list() const {
    return plist_;
  }
  
  // Set/Get mesh type - RECTANGULAR, GENERAL (See MeshDefs.hh)
  void set_mesh_type(const Mesh_type mesh_type) {
    mesh_type_ = mesh_type;
  }
  Mesh_type mesh_type() const {
    return mesh_type_;
  }

  bool is_logical() const { return logical_; }

  //
  // General mesh information
  // ------------------------

  // Number of entities of any kind (cell, face, node) and in a
  // particular category (OWNED, GHOST, ALL)
  virtual
  unsigned int num_entities(const Entity_kind kind,
                            const Parallel_type ptype) const = 0;

  // Were optional edges initialized?
  virtual
  bool valid_edges() const { return false; }


  //
  // Entity meta-data
  // ----------------
  // NOTE: Anything here could be cached if need be. --etc

  // Get parallel type of entity - OWNED, GHOST, ALL (See MeshDefs.hh)
  virtual
  Parallel_type entity_get_ptype(const Entity_kind kind,
                                 const Entity_ID entid) const = 0;

  // Parent entity in the source mesh if mesh was derived from another mesh
  virtual
  Entity_ID entity_get_parent(const Entity_kind kind, const Entity_ID entid) const;

  // Get cell type - UNKNOWN, TRI, QUAD, etc. See MeshDefs.hh
  virtual
  Cell_type cell_get_type(const Entity_ID cellid) const = 0;

  // Cell type name
  std::string cell_type_to_name(const Cell_type type);

  // Global ID of any entity
  virtual
  Entity_ID GID(const Entity_ID lid, const Entity_kind kind) const = 0;


  //
  // Mesh Entity Adjacencies
  //-------------------------
  // NOTE: Anything here could be cached if need be. --etc
  //
  // TODO: Things that aren't cached should be aliased, i.e. the virtual one
  // should still be the *_internal() call and this should simply call that
  // method.  That would make the deriving mesh interfaces a bit nicer
  // (i.e. currently they have cell_get_nodes() and
  // cell_get_faces_internal() ). --etc


  // Downward Adjacencies
  //---------------------

  // Get number of faces of a cell.
  //
  // On a distributed mesh, this will count all the faces of the
  // cell, OWNED or GHOST.
  unsigned int cell_get_num_faces(const Entity_ID cellid) const;
  unsigned int cell_get_max_faces() const;
  unsigned int cell_get_max_nodes() const;
  unsigned int cell_get_max_edges() const;

  // Get faces of a cell.
  //
  // On a distributed mesh, this will return all the faces of the
  // cell, OWNED or GHOST. If ordered = true, the faces will be
  // returned in a standard order according to Exodus II convention
  // for standard cells; in all other situations (ordered = false or
  // non-standard cells), the list of faces will be in arbitrary order
  //
  // EXTENSIONS: MSTK FRAMEWORK: by the way the parallel partitioning,
  // send-receive protocols and mesh query operators are designed, a side 
  // effect of this is that master and ghost entities will have the same
  // hierarchical topology. 
  void cell_get_faces(const Entity_ID cellid,
                      Entity_ID_List *faceids,
                      const bool ordered=false) const;

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
  void cell_get_faces_and_dirs(const Entity_ID cellid,
                               Entity_ID_List *faceids,
                               std::vector<int> *face_dirs,
                               const bool ordered=false) const;

  // Get the bisectors, i.e. vectors from cell centroid to face centroids.
  virtual
  void cell_get_faces_and_bisectors(const Entity_ID cellid,
                                    Entity_ID_List *faceids,
                                    std::vector<AmanziGeometry::Point> *bisectors,
                                    const bool ordered=false) const;

  // Get edges of a cell
  void cell_get_edges(const Entity_ID cellid, Entity_ID_List *edgeids) const;

  // Get edges and dirs of a 2D cell.
  //
  // This is to make the code cleaner for integrating over the cell in 2D
  // where faces and edges are identical but integrating over the cells using
  // face information is more cumbersome (one would have to take the face
  // normals, rotate them and then get a consistent edge vector)
  void cell_2D_get_edges_and_dirs(const Entity_ID cellid,
                                  Entity_ID_List *edgeids,
                                  std::vector<int> *edge_dirs) const;

  // Get nodes of a cell
  virtual
  void cell_get_nodes(const Entity_ID cellid,
                      Entity_ID_List *nodeids) const = 0;

  // Get edges of a face and directions in which the face uses the edges.
  //
  // In 3D, edge direction is 1 when it is oriented counter clockwise
  // with respect to the face natural normal.
  //
  // On a distributed mesh, this will return all the edges of the
  // face, OWNED or GHOST. If ordered = true, the edges will be
  // returned in a ccw order around the face as it is naturally defined.
  //
  // IMPORTANT NOTE IN 2D CELLS: In meshes where the cells are two
  // dimensional, faces and edges are identical. For such cells, this
  // operator will return a single edge and a direction of 1. However,
  // this direction cannot be relied upon to compute, say, a contour
  // integral around the 2D cell.
  void face_get_edges_and_dirs(const Entity_ID faceid,
                               Entity_ID_List *edgeids,
                               std::vector<int> *edge_dirs,
                               const bool ordered=false) const;

  // Get the local index of a face edge in a cell edge list
  // Example:
  //
  // face_get_edges(face=5) --> {20, 21, 35, 9, 10}
  // cell_get_edges(cell=18) --> {1, 2, 3, 5, 8, 9, 10, 13, 21, 35, 20, 37, 40}
  // face_to_cell_edge_map(face=5,cell=18) --> {10, 8, 9, 5, 6}
  void face_to_cell_edge_map(const Entity_ID faceid, const Entity_ID cellid,
                             std::vector<int> *map) const;

  // Get nodes of face
  // On a distributed mesh, all nodes (OWNED or GHOST) of the face
  // are returned
  // In 3D, the nodes of the face are returned in ccw order consistent
  // with the face normal. 
  // In 2D, nfnodes is 2.
  // IMPORTANT NOTE FOR MSTK FRAMEWORK: The first node in 3D is the 
  // starting node of the first edge oriented ccw. Thus, the ccw order
  // of vertices and edges is as follows: v0, e0, v1, e1, v2, ..., eN.
  virtual
  void face_get_nodes(const Entity_ID faceid,
                      Entity_ID_List *nodeids) const = 0;

  // Get nodes of edge
  virtual
  void edge_get_nodes(const Entity_ID edgeid,
                      Entity_ID *nodeid0, Entity_ID *nodeid1) const = 0;


  // Upward adjacencies
  //-------------------

  // Cells of type 'ptype' connected to a node
  //
  // The order of cells is not guaranteed to be the same for corresponding
  // nodes on different processors
  virtual
  void node_get_cells(const Entity_ID nodeid,
                      const Parallel_type ptype,
                      Entity_ID_List *cellids) const = 0;

  // Faces of type 'ptype' connected to a node
  //
  // The order of faces is not guarnateed to be the same for corresponding
  // nodes on different processors
  virtual
  void node_get_faces(const Entity_ID nodeid,
                      const Parallel_type ptype,
                      Entity_ID_List *faceids) const = 0;

  // Get faces of ptype of a particular cell that are connected to the
  // given node
  //
  // The order of faces is not guarnateed to be the same for corresponding
  // nodes on different processors
  virtual
  void node_get_cell_faces(const Entity_ID nodeid,
                           const Entity_ID cellid,
                           const Parallel_type ptype,
                           Entity_ID_List *faceids) const = 0;

  // Faces of type 'ptype' connected to an edge
  //
  // The order of faces is not guaranteed to be the same for corresponding
  // edges on different processors
  virtual
  void edge_get_faces(const Entity_ID edgeid,
                      const Parallel_type ptype,
                      Entity_ID_List *faceids) const { AMANZI_ASSERT(false); }

  // Cells of type 'ptype' connected to an edge
  //
  // The order of cells is not guaranteed to be the same for corresponding
  // edges on different processors
  virtual
  void edge_get_cells(const Entity_ID edgeid,
                      const Parallel_type ptype,
                      Entity_ID_List *cellids) const = 0;

  // Cells connected to a face
  //
  // The cells are returned in no particular order. Also, the order of cells
  // is not guaranteed to be the same for corresponding faces on different
  // processors
  void face_get_cells(const Entity_ID faceid,
                      const Parallel_type ptype,
                      Entity_ID_List *cellids) const;


  // Same level adjacencies
  //-----------------------

  // Face connected neighboring cells of given cell of a particular ptype
  // (e.g. a hex has 6 face neighbors)
  //
  // The order in which the cellids are returned cannot be
  // guaranteed in general except when ptype = ALL, in which case
  // the cellids will correcpond to cells across the respective
  // faces given by cell_get_faces
  virtual
  void cell_get_face_adj_cells(const Entity_ID cellid,
                               const Parallel_type ptype,
                               Entity_ID_List *fadj_cellids) const = 0;

  // Node connected neighboring cells of given cell
  // (a hex in a structured mesh has 26 node connected neighbors)
  // The cells are returned in no particular order
  virtual
  void cell_get_node_adj_cells(const Entity_ID cellid,
                               const Parallel_type ptype,
                               Entity_ID_List *nadj_cellids) const = 0;


  //
  // Column information
  //-------------------

  // Special adjacency information for geological domains with a
  // semi-structured mesh. Will return -1 if there is no suitable cell
  // The code currently makes the assumption that the "bottom" of the
  // is a flat surface in the XY plane. It then builds up information
  // about the cell above and cell below for each cell based on the
  // orientation of the face normals w.r.t the z direction. If the
  // mesh is highly warped, this could lead to ambiguities. Also,
  // intersecting columns in an unstructured mesh will lead to an
  // exception being thrown. These data structures are never populated
  // if these operators are never called. The above and below cells
  // are computed for all cells the first time one of these routines
  // is called and then cached

  //
  virtual
  int build_columns(const std::string& Msetname) const;

  // Build columns over the entire mesh. The columns are defined by
  // starting from boundary faces which have a negative-z-direction
  // normal, then collecting cells and faces while traveling downward
  // through the columns.

  virtual
  int build_columns() const;
  
  // Number of columns in mesh - must call build_columns before calling
  int num_columns(bool ghosted=false) const;

  // Given a column ID, get the cells of the column - must call build_columns before calling
  const Entity_ID_List& cells_of_column(const int columnID_) const;

  // Given a column ID, get the cells of the column - must call build_columns before calling
  const Entity_ID_List& faces_of_column(const int columnID_) const;

  // Given a cell, get its column ID - must call build_columns before calling
  int column_ID(const Entity_ID cellid) const;

  // Given a cell, get the id of the cell above it in the column - must call build_columns before calling
  Entity_ID cell_get_cell_above(const Entity_ID cellid) const;

  // Given a cell, get the id of the cell below it in the column - must call build_columns before calling
  Entity_ID cell_get_cell_below(const Entity_ID cellid) const;

  // Given a node, get the id of the node above it in the column - must call build_columns before calling
  Entity_ID node_get_node_above(const Entity_ID nodeid) const;


  //
  // Mesh entity geometry
  // --------------------

  // Node coordinates - 3 in 3D and 2 in 2D
  virtual
  void node_get_coordinates(const Entity_ID nodeid,
                            AmanziGeometry::Point *ncoord) const = 0;

  // Face coordinates - conventions same as face_to_nodes call
  // Number of nodes is the vector size.
  virtual
  void face_get_coordinates(const Entity_ID faceid,
                            std::vector<AmanziGeometry::Point> *fcoords) const = 0;

  // Coordinates of cells in standard order (Exodus II convention)
  //
  // Standard convention works only for standard cell types in 3D!  For a
  // general polyhedron this will return the node coordinates in arbitrary
  // order.
  virtual
  void cell_get_coordinates(const Entity_ID cellid,
                            std::vector<AmanziGeometry::Point> *ccoords) const = 0;

  // Volume/Area of cell
  double cell_volume(const Entity_ID cellid, const bool recompute=false) const;

  // Area/length of face
  double face_area(const Entity_ID faceid, const bool recompute=false) const;

  // Length of edge
  double edge_length(const Entity_ID edgeid, const bool recompute=false) const;

  // Centroid of cell (center of gravity not just average of node coordinates)
  //
  // The cell centroid is computed as the volume weighted average of the
  // centroids of tetrahedra from a symmetric tetrahedral
  // decomposition of the cell. The tetrahedral decomposition is
  // formed by connecting the cell center (average of cell nodes), a
  // face center (average of face nodes) and the two nodes of an edge
  // of the face
  AmanziGeometry::Point cell_centroid(const Entity_ID cellid,
                                      const bool recompute=false) const;

  // Centroid of face (center of gravity not just the average of node coordinates)
  //
  // The face centroid is computed as the area weighted average of the
  // centroids of the triangles from a symmetric triangular
  // decomposition of the face. Each triangular facet is formed by the
  // connecting the face center (average of face nodes) to the two
  // nodes of an edge of the face
  AmanziGeometry::Point face_centroid(const Entity_ID faceid,
                                      const bool recompute=false) const;

  // Centroid of edge
  AmanziGeometry::Point edge_centroid(const Entity_ID edgeid) const;


  // Normal to face
  //
  // The vector is normalized and then weighted by the area of the face
  //
  // If recompute is TRUE, then the normal is recalculated using current
  // face coordinates but not stored. (If the recomputed normal must be
  // stored, then call recompute_geometric_quantities).
  //
  // If cellid is not specified, the normal is the natural normal of
  // the face. This means that at boundaries, the normal may point in
  // or out of the domain depending on how the face is defined. On the
  // other hand, if cellid is specified, the normal is the outward
  // normal with respect to the cell. In planar and solid meshes, the
  // normal with respect to the cell on one side of the face is just
  // the negative of the normal with respect to the cell on the other
  // side. In general surfaces meshes, this will not be true at C1
  // discontinuities
  //
  // if cellid is specified, then orientation returns the direction of
  // the natural normal of the face with respect to the cell (1 is
  // pointing out of the cell and -1 pointing in)
  AmanziGeometry::Point face_normal(const Entity_ID faceid,
                                    const bool recompute=false,
                                    const Entity_ID cellid=-1,
                                    int *orientation=NULL) const;

  // Edge vector
  //
  // Vector length equals to edge length.
  //
  // If recompute is TRUE, then the vector is recalculated using current
  // edge coordinates but not stored. (If the recomputed vector must be
  // stored, then call recompute_geometric_quantities).
  //
  // If pointid is not specified, the vector is the natural direction of
  // the edge (from point0 to point1).  On the other hand, if pointid
  // is specified (has to be a point of the face), the vector is from
  // specified point to opposite point of edge.
  //
  // if pointid is specified, then orientation returns the direction of
  // the natural direction of the edge with respect to the point (1 is
  // away from the point and -1 is towards)
  AmanziGeometry::Point edge_vector(const Entity_ID edgeid,
                                    const bool recompute=false,
                                    const Entity_ID pointid=-1,
                                    int *orientation=NULL) const;

  // Is a point in a given cell?
  bool point_in_cell(const AmanziGeometry::Point &p,
                     const Entity_ID cellid) const;


  //
  // Mesh modification
  //------------------

  // Set coordinates of node
  virtual
  void node_set_coordinates(const Entity_ID nodeid,
                            const AmanziGeometry::Point ncoord) = 0;

  // Set coordinates of node
  virtual
  void node_set_coordinates(const Entity_ID nodeid,
                            const double *ncoord) = 0;

  // Just move the mesh.  Returns 0 if negative cell volumes generated, 1
  // otherwise.
  //
  // This is a rudimentary capability that requires ghosts nodes
  // also to be deformed. Amanzi does not have any built-in parallel 
  // communication capabilities, other than Trilinos Epetra object
  // communication or raw MPI. 
  virtual
  int deform(const Entity_ID_List& nodeids,
             const AmanziGeometry::Point_List& new_positions);
  
  // Deform the mesh by moving given nodes to given coordinates, one at a
  // time, ensuring validity at all points through the deformation, which is a
  // really bad idea for deformations that move large numbers of nodes more
  // than epsilon.  Really.  Just don't use this routine, it almost definitely
  // isn't what you want.
  //
  // Deform the mesh by moving given nodes to given coordinates.  If the flag
  // keep_valid is true, then the nodes are moved only as much as possible
  // without making the mesh invalid.  The final positions of the nodes is
  // returned in final_positions
  //
  // This is a rudimentary capability that requires ghosts nodes
  // also to be deformed. Amanzi does not have any built-in parallel 
  // communication capabilities, other than Trilinos Epetra object
  // communication or raw MPI. 
  virtual
  int deform(const Entity_ID_List& nodeids,
             const AmanziGeometry::Point_List& new_positions,
             const bool keep_valid,
             AmanziGeometry::Point_List *final_positions);

  // Deform a mesh so that cell volumes conform as closely as possible
  // to target volumes without dropping below the minimum volumes.  If
  // move_vertical = true, nodes will be allowed to move only in the
  // vertical direction (right now arbitrary node movement is not allowed)
  // Nodes in any set in the fixed_sets will not be permitted to move.
  virtual
  int deform(const std::vector<double>& target_cell_volumes_in,
             const std::vector<double>& min_cell_volumes_in,
             const Entity_ID_List& fixed_nodes,
             const bool move_vertical) = 0;


  // Synchronize node positions across processors
  void update_ghost_node_coordinates();


  //
  // Epetra maps
  //------------
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
  
  // Get cell map
  virtual
  const Epetra_Map& cell_map(bool include_ghost) const = 0;

  // Get face map
  virtual
  const Epetra_Map& face_map(bool include_ghost) const = 0;

  // Get edge map
  // dummy implementation so that frameworks can skip or overwrite
  virtual
  const Epetra_Map& edge_map(bool include_ghost) const
  {
    Errors::Message mesg("Edges are not implemented in this framework.");
    Exceptions::amanzi_throw(mesg);
    throw(mesg);
  }

  // Get node map
  virtual
  const Epetra_Map& node_map(bool include_ghost) const = 0;

  // Get map of only exterior faces
  virtual
  const Epetra_Map& exterior_face_map(bool include_ghost) const = 0;

  // Get Epetra importer that will allow apps to import values from a
  // Epetra vector defined on all owned faces into an Epetra vector
  // defined only on exterior faces
  virtual
  const Epetra_Import& exterior_face_importer(void) const = 0;

  // Get map of only exterior nodes
  virtual
  const Epetra_Map& exterior_node_map(bool include_ghost) const = 0;

  //
  // Mesh Sets for ICs, BCs, Material Properties and whatever else
  //--------------------------------------------------------------

  // Is this is a valid ID of a set containing entities of 'kind'
  bool valid_set_id(const Set_ID setid,
                    const Entity_kind kind) const;

  // Is this is a valid ID of a set containing entities of 'kind'
  bool valid_set_name(const std::string setname,
                      const Entity_kind kind) const;

  // Get set ID from set name - returns 0 if no match is found
  Set_ID set_id_from_name(const std::string setname) const;

  // Get set name from set ID - returns 0 if no match is found
  std::string set_name_from_id(const Set_ID setid) const;

  // Get number of entities of type 'category' in set
  virtual
  unsigned int get_set_size(const Set_ID setid,
                            const Entity_kind kind,
                            const Parallel_type ptype) const;

  virtual
  unsigned int get_set_size(const std::string setname,
                            const Entity_kind kind,
                            const Parallel_type ptype) const;

  // Get list of entities of type 'category' in set by set name.
  // -- original interface. The returned volume fractions are ignored.
  virtual
  void get_set_entities(const std::string setname,
                        const Entity_kind kind,
                        const Parallel_type ptype,
                        Entity_ID_List *entids) const;

  // -- deprecated interface
  virtual
  void get_set_entities(const Set_ID setid, 
                        const Entity_kind kind, 
                        const Parallel_type ptype, 
                        Entity_ID_List *entids) const;

  // -- new interface. Since not all regions support volume fractions 
  // (vofs), this vector is optional and could be empty.
  virtual
  void get_set_entities_and_vofs(const std::string setname,
                                 const Entity_kind kind,
                                 const Parallel_type ptype,
                                 Entity_ID_List *entids,
                                 std::vector<double> *vofs) const = 0;

  //
  // High-order mesh
  //----------------

  // Geometry of a curved edge/face is defined by a derived mesh class 
  // from the list of returned interior nodes. For a linear (default)
  // mesh, these functions return the empty list.
  virtual
  void edge_get_ho_nodes(Entity_ID edgeid,
                         AmanziGeometry::Point_List *nodes) const { nodes->clear(); }

  virtual
  void face_get_ho_nodes(Entity_ID faceid,
                         AmanziGeometry::Point_List *nodes) const { nodes->clear(); }

  //
  // Miscellaneous functions
  //------------------------

  // Write to file
  virtual
  void write_to_exodus_file(const std::string filename) const = 0;

protected:
  // Helper function to build columns
  int build_single_column_(int colnum, Entity_ID top_face) const;

  // Beginning of new interface to regions using the base mesh.
  void get_set_entities_box_vofs_(
      Teuchos::RCP<const AmanziGeometry::Region> region,
      const Entity_kind kind, 
      const Parallel_type ptype, 
      std::vector<Entity_ID>* setents,
      std::vector<double> *vofs) const;


  //
  // Cache management
  //-----------------

  // Virtual methods to fill the cache with geometric quantities.
  //
  // Default implementations use _internal() methods below.
  virtual
  int compute_cell_geometric_quantities_() const;
  virtual
  int compute_face_geometric_quantities_() const;
  virtual
  int compute_edge_geometric_quantities_() const;

  // Virtual methods to fill the cache with topological quantities.
  //
  // Default implementations use _internal() methods below.
  virtual
  void cache_cell_face_info_() const;
  virtual
  void cache_cell2edge_info_() const;
  virtual
  void cache_face2edge_info_() const;

  // Virtual methods for mesh geometry.
  //
  // These are virtual and therefore slightly expensive, so they
  // should be used once to populate the cache and not again.  They
  // have the same concepts behind them as the non- _internal()
  // versions.  Non- _internal() versions are not virtual and access
  // the cache; these do the real work and are implemented by the mesh
  // implementation.

  // Get faces of a cell and directions in which it is used.
  virtual
  void cell_get_faces_and_dirs_internal_(const Entity_ID cellid,
                                         Entity_ID_List *faceids,
                                         std::vector<int> *face_dirs,
                                         const bool ordered=false) const = 0;

  // Cells connected to a face
  virtual
  void face_get_cells_internal_(const Entity_ID faceid,
                                const Parallel_type ptype,
                                Entity_ID_List *cellids) const = 0;

  // edges of a face
  virtual
  void face_get_edges_and_dirs_internal_(
      const Entity_ID faceid,
      Entity_ID_List *edgeids,
      std::vector<int> *edge_dirs,
      const bool ordered=true) const = 0;

  // edges of a cell
  virtual
  void cell_get_edges_internal_(const Entity_ID cellid,
                                Entity_ID_List *edgeids) const = 0;

  // edges and directions of a 2D cell
  virtual
  void cell_2D_get_edges_and_dirs_internal_(
      const Entity_ID cellid,
      Entity_ID_List *edgeids,
      std::vector<int> *edge_dirs) const = 0;

  // Virtual methods to fill the cache with geometric quantities.
  //
  // Convenience methods that wrap multiple calls.
  //
  // These are declared const since they do not modify the
  // mesh but just modify cached variables declared as mutable
  virtual
  int compute_cell_geometry_(const Entity_ID cellid,
                             double *volume,
                             AmanziGeometry::Point *centroid) const;
  virtual
  int compute_face_geometry_(const Entity_ID faceid,
                             double *area,
                             AmanziGeometry::Point *centroid,
                             std::vector<AmanziGeometry::Point> *normals) const;

  virtual
  int compute_edge_geometry_(const Entity_ID edgeid,
                             double *length,
                             AmanziGeometry::Point *edge_vector) const;

 public:
  void PrintMeshStatistics() const;

 protected:
  Comm_ptr_type comm_;
  Teuchos::RCP<const AmanziGeometry::GeometricModel> geometric_model_;
  Teuchos::RCP<const Teuchos::ParameterList> plist_;
  Teuchos::RCP<const VerboseObject> vo_;

  Mesh_type mesh_type_;
  unsigned int manifold_dim_, space_dim_;

  bool logical_;
  Teuchos::RCP<const Mesh> parent_;
  
  // the cache
  // -- geometry
  mutable std::vector<double> cell_volumes_, face_areas_, edge_lengths_;
  mutable std::vector<AmanziGeometry::Point> cell_centroids_, face_centroids_;

  // -- Have to account for the fact that a "face" for a non-manifold
  // surface mesh can have more than one cell connected to
  // it. Therefore, we have to store as many normals for a face as
  // there are cells connected to it. For a given face, its normal to
  // face_get_cells()[i] is face_normals_[i]
  mutable std::vector<std::vector<AmanziGeometry::Point>> face_normals_;

  mutable std::vector<AmanziGeometry::Point> edge_vectors_;

  // -- column information, only created if columns are requested
  mutable Entity_ID_List cell_cellabove_, cell_cellbelow_, node_nodeabove_;
  mutable std::vector<Entity_ID_List> columns_cells_;
  mutable std::vector<Entity_ID_List> columns_faces_;
  mutable std::vector<Entity_ID> columnsID_;
  mutable int num_owned_cols_;
  mutable bool columns_built_;
  
  // -- topology
  mutable std::vector<Entity_ID_List> cell_face_ids_;
  mutable std::vector< std::vector<int> > cell_face_dirs_;  // 1 or -1

  // 1s complement if face is pointing out of cell; cannot use 0 as
  // cellid can be 0
  mutable std::vector<Entity_ID_List> face_cell_ids_;

  mutable std::vector< std::vector<Parallel_type> > face_cell_ptype_;
  mutable std::vector<Entity_ID_List> cell_edge_ids_;
  mutable std::vector< std::vector<int> > cell_2D_edge_dirs_;
  mutable std::vector<Entity_ID_List> face_edge_ids_;
  mutable std::vector< std::vector<int> > face_edge_dirs_;

  // -- flags to indicate what part of cache is up-to-date
  mutable bool cell2face_info_cached_, face2cell_info_cached_;
  mutable bool cell2edge_info_cached_, face2edge_info_cached_;
  mutable bool cell_geometry_precomputed_, face_geometry_precomputed_, edge_geometry_precomputed_;

  // -- region data
  mutable std::map<std::string, std::vector<int> > region_ids;
  mutable std::map<std::string, std::vector<double> > region_vofs;

  // probably should not be mutable?  these should be set by constructor and not changed! --etc
  mutable bool faces_requested_, edges_requested_;

  // friend classes change the cache?  why is this necessary? --etc
  friend class MeshEmbeddedLogical;

  // fast search tools
  mutable bool kdtree_faces_initialized_;
  mutable KDTree kdtree_faces_;
};


#include "Mesh_impl.hh"

}  // namespace AmanziMesh
}  // namespace Amanzi

#endif

