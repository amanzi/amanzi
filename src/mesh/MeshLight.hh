/*
  Mesh

  Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Konstantin Lipnikov
  Rao Garimella

  The smallest mesh class with basic connectivity and geometry
*/

#ifndef AMANZI_MESH_LIGHT_HH_
#define AMANZI_MESH_LIGHT_HH_

// Amanzi
#include "Point.hh"

// Amanzi::AmanziMesh
#include "MeshDefs.hh"

namespace Amanzi {
namespace AmanziMesh {

class MeshLight {
 public:
  MeshLight()
    : space_dim_(-1),
      faces_requested_(false),
      edges_requested_(false),
      cell2face_info_cached_(false),
      cell2edge_info_cached_(false),
      face2edge_info_cached_(false),
      face2cell_info_cached_(false),
      cell_geometry_precomputed_(false),
      face_geometry_precomputed_(false),
      edge_geometry_precomputed_(false) {};

  // initializing mesh
  void BuildCache();

  // base class functionality
  void set_space_dimension(unsigned int dim) { space_dim_ = dim; }
  unsigned int space_dimension() const { return space_dim_; }

  void set_manifold_dimension(const unsigned int dim) { manifold_dim_ = dim; }
  unsigned int manifold_dimension() const { return manifold_dim_; }

  virtual bool valid_edges() const { return false; }

  // ---------------------
  // Downward connectivity
  // ---------------------

  // Get faces of a cell
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
  //
  // new API: cache should be pre-build, e.g. in mesh constructor
  // so no additional checks is needed
  const Entity_ID_List& cell_get_faces(const Entity_ID c) const {
    return cell_face_ids_[c];
  }

  // Get directions in which the cell uses the face
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
  const std::vector<int>& cell_get_face_dirs(const Entity_ID c) const {
    return cell_face_dirs_[c];
  }

  // Get edges of a cell
  const Entity_ID_List& cell_get_edges(const Entity_ID c) const {
    return cell_edge_ids_[c];
  }

  virtual void cell_get_nodes(const Entity_ID c, Entity_ID_List *nodes) const = 0;

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
  void face_get_edges_and_dirs(
       const Entity_ID f,
       Entity_ID_List *edges,
       std::vector<int> *dirs,
       const bool ordered = false) const;

  // Get nodes of face
  // On a distributed mesh, all nodes (OWNED or GHOST) of the face
  // are returned
  // In 3D, the nodes of the face are returned in ccw order consistent
  // with the face normal. 
  // In 2D, nfnodes is 2.
  // IMPORTANT NOTE FOR MSTK FRAMEWORK: The first node in 3D is the 
  // starting node of the first edge oriented ccw. Thus, the ccw order
  // of vertices and edges is as follows: v0, e0, v1, e1, v2, ..., eN.
  virtual void face_get_nodes(const Entity_ID f, Entity_ID_List *nodes) const = 0;

  virtual void edge_get_nodes(const Entity_ID e, Entity_ID* n0, Entity_ID* n1) const = 0;

  // Get the local index of a face edge in a cell edge list
  //
  // Example:
  // face_get_edges(face=5) --> {20, 21, 35, 9, 10}
  // cell_get_edges(cell=18) --> {1, 2, 3, 5, 8, 9, 10, 13, 21, 35, 20, 37, 40}
  // face_to_cell_edge_map(face=5,cell=18) --> {10, 8, 9, 5, 6}
  void face_to_cell_edge_map(
       const Entity_ID f, const Entity_ID c, std::vector<int> *map) const;


  //--------------------
  // Upward connectivity 
  //--------------------
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

  // The cells are returned in no particular order. Also, the order of cells
  // is not guaranteed to be the same for corresponding faces on different
  // processors
  void face_get_cells(
       const Entity_ID f,
       const Parallel_type ptype,
       Entity_ID_List *cells) const;


  // --------
  // Geometry
  // --------
  // Centroid of cell (center of gravity not just average of node coordinates)
  //
  // The cell centroid is computed as the volume weighted average of the
  // centroids of tetrahedra from a symmetric tetrahedral
  // decomposition of the cell. The tetrahedral decomposition is
  // formed by connecting the cell center (average of cell nodes), a
  // face center (average of face nodes) and the two nodes of an edge
  // of the face
  AmanziGeometry::Point cell_centroid(
      const Entity_ID c, bool recompute = false) const;

  double cell_volume(
      const Entity_ID c, const bool recompute = false) const;

  AmanziGeometry::Point face_centroid(
      const Entity_ID f, const bool recompute = false) const;

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
  AmanziGeometry::Point face_normal(
      const Entity_ID f,
      const bool recompute = false,
      const Entity_ID cellid = -1,
      int *orientation = NULL) const;

  double face_area(
      const Entity_ID f, const bool recompute = false) const;

  AmanziGeometry::Point edge_centroid(const Entity_ID e) const;

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
  AmanziGeometry::Point edge_vector(
      const Entity_ID e,
      const bool recompute = false,
      const Entity_ID pointid = -1,
      int *orientation = NULL) const;

  double edge_length(
      const Entity_ID e, const bool recompute = false) const;

  // coordinates
  virtual void node_get_coordinates(
          const Entity_ID v, AmanziGeometry::Point *ncoord) const = 0;


  // ------
  // Counts
  // ------
  // Entities of kind (cell, face, node) in a particular category (OWNED, GHOST, ALL)
  virtual unsigned int num_entities(
          const Entity_kind kind, const Parallel_type ptype) const = 0;


  // ----------------
  // Entity meta-data
  // ----------------
  virtual Parallel_type entity_get_ptype(
          const Entity_kind kind, const Entity_ID entid) const = 0;

  unsigned int cell_get_num_faces(const Entity_ID c) const;
  unsigned int cell_get_max_faces() const;
  unsigned int cell_get_max_edges() const;
  unsigned int cell_get_max_nodes() const;

 protected:
  // Cache filling methods use _internal() virtual functions.
  void cache_cell2face_info_() const;
  void cache_cell2edge_info_() const;
  void cache_face2edge_info_() const;

  int compute_cell_geometric_quantities_() const;
  int compute_face_geometric_quantities_() const;
  int compute_edge_geometric_quantities_() const;

  // These are virtual and therefore slightly expensive, so they
  // should be used once to populate the cache and not again.  They
  // have the same concepts behind them as the non- _internal()
  // versions.  Non- _internal() versions are not virtual and access
  // the cache; these do the real work and are implemented by the mesh
  // implementation.

  // faces of a cell and directions in which it is used
  virtual void cell_get_faces_and_dirs_internal_(
          const Entity_ID c,
          Entity_ID_List *faces,
          std::vector<int> *fdirs,
          const bool ordered = false) const = 0;

  // edges of a cell
  virtual void cell_get_edges_internal_(
          const Entity_ID cellid,
          Entity_ID_List *edgeids) const = 0;

  // edges of a face
  virtual void face_get_edges_and_dirs_internal_(
          const Entity_ID faceid,
          Entity_ID_List *edgeids,
          std::vector<int> *edge_dirs,
          const bool ordered = true) const = 0;

  // cells of a face
  virtual void face_get_cells_internal_(
          const Entity_ID faceid,
          const Parallel_type ptype,
          Entity_ID_List *cellids) const = 0;

  // geometries
  virtual int compute_cell_geometry_(
          const Entity_ID cellid,
          double *volume,
          AmanziGeometry::Point *centroid) const = 0;

  virtual int compute_face_geometry_(
          const Entity_ID faceid,
          double *area,
          AmanziGeometry::Point *centroid,
          std::vector<AmanziGeometry::Point> *normals) const = 0;

  virtual int compute_edge_geometry_(
          const Entity_ID e,
          double *length,
          AmanziGeometry::Point *edge_vector) const = 0;

 protected:
  unsigned int space_dim_;
  unsigned int manifold_dim_;

  mutable bool faces_requested_;
  mutable bool edges_requested_;

  // cache: c -> f
  mutable bool cell2face_info_cached_;
  mutable std::vector<Entity_ID_List> cell_face_ids_;
  mutable std::vector<std::vector<int> > cell_face_dirs_;  // 1 or -1

  // cache: c -> e
  mutable bool cell2edge_info_cached_;
  mutable std::vector<Entity_ID_List> cell_edge_ids_;

  // cache: f -> e
  mutable bool face2edge_info_cached_;
  mutable std::vector<Entity_ID_List> face_edge_ids_;
  mutable std::vector<std::vector<int> > face_edge_dirs_;

  // cache: f -> c
  mutable bool face2cell_info_cached_;
  // 1s complement if face is pointing out of cell; cannot use 0 as cellid can be 0
  mutable std::vector<Entity_ID_List> face_cell_ids_;
  mutable std::vector<std::vector<Parallel_type> > face_cell_ptype_;

  // cache: cells
  mutable bool cell_geometry_precomputed_;
  mutable std::vector<double> cell_volumes_;
  mutable std::vector<AmanziGeometry::Point> cell_centroids_;

  // cache: faces
  mutable bool face_geometry_precomputed_;
  mutable std::vector<double> face_areas_;
  mutable std::vector<AmanziGeometry::Point> face_centroids_;
  // Have to account for the fact that a "face" for a non-manifold
  // surface mesh can have more than one cell connected to
  // it. Therefore, we have to store as many normals for a face as
  // there are cells connected to it. For a given face, its normal to
  // face_get_cells()[i] is face_normals_[i]
  mutable std::vector<std::vector<AmanziGeometry::Point>> face_normals_;

  // cache: edges
  mutable bool edge_geometry_precomputed_;
  mutable std::vector<double> edge_lengths_;
  mutable std::vector<AmanziGeometry::Point> edge_vectors_;
};

}  // namespace AmanziMesh
}  // namespace Amanzi

#endif
