/*
  Mesh

  Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Konstantin Lipnikov
  Rao Garimella

  The smallest abstract mesh class. 
*/

#ifndef AMANZI_MESH_STATIC_MESH_MINI_HH_
#define AMANZI_MESH_STATIC_MESH_MINI_HH_

// Amanzi
#include "Point.hh"

// Amanzi::AmanziMesh
#include "MeshDefs.hh"

namespace Amanzi {
namespace AmanziMesh {

class StaticMeshMini {
 public:
   StaticMeshMini() : space_dim_(-1) {};

  // base class functionality
  void set_space_dimension(unsigned int dim) { space_dim_ = dim; }
  unsigned int space_dimension() const { return space_dim_; }

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
  virtual void cell_get_faces(
          const Entity_ID c,
          Entity_ID_List *faces,
          const bool ordered = false) const {
    cell_get_faces_and_dirs(c, faces, NULL, ordered);
  }

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
  virtual void cell_get_faces_and_dirs(
          const Entity_ID c,
          Entity_ID_List *faces,
          std::vector<int> *dirs,
          bool ordered = false) const = 0;

  virtual void cell_get_edges(const Entity_ID c, Entity_ID_List *edges) const = 0;

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
  virtual void face_get_edges_and_dirs(
          const Entity_ID f,
          Entity_ID_List *edges,
          std::vector<int> *dirs,
          const bool ordered = false) const = 0;

  // Get nodes of face
  // On a distributed mesh, all nodes (OWNED or GHOST) of the face
  // are returned
  // In 3D, the nodes of the face are returned in ccw order consistent
  // with the face normal. 
  // In 2D, nfnodes is 2.
  // IMPORTANT NOTE FOR MSTK FRAMEWORK: The first node in 3D is the 
  // starting node of the first edge oriented ccw. Thus, the ccw order
  // of vertices and edges is as follows: v0, e0, v1, e1, v2, ..., eN.
  virtual void face_get_nodes(
          const Entity_ID f, Entity_ID_List *nodes) const = 0;

  virtual void edge_get_nodes(
          const Entity_ID e, Entity_ID* n0, Entity_ID* n1) const = 0;

  // --------------------
  // Geometric objects
  // --------------------
  // Centroid of cell (center of gravity not just average of node coordinates)
  //
  // The cell centroid is computed as the volume weighted average of the
  // centroids of tetrahedra from a symmetric tetrahedral
  // decomposition of the cell. The tetrahedral decomposition is
  // formed by connecting the cell center (average of cell nodes), a
  // face center (average of face nodes) and the two nodes of an edge
  // of the face
  virtual AmanziGeometry::Point cell_centroid(
          const Entity_ID c, const bool recompute = false) const = 0;

  virtual double cell_volume(
          const Entity_ID c, const bool recompute = false) const = 0;

  virtual AmanziGeometry::Point face_centroid(
          const Entity_ID f, const bool recompute = false) const = 0;

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
  virtual AmanziGeometry::Point face_normal(
          const Entity_ID f,
          const bool recompute = false,
          const Entity_ID cellid = -1,
          int *orientation = NULL) const = 0;

  virtual double face_area(
          const Entity_ID f, const bool recompute = false) const = 0;

  virtual AmanziGeometry::Point edge_centroid(const Entity_ID e) const = 0;

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
  virtual AmanziGeometry::Point edge_vector(
          const Entity_ID e,
          const bool recompute = false,
          const Entity_ID pointid = -1,
          int *orientation = NULL) const = 0;

  virtual double edge_length(
          const Entity_ID e, const bool recompute = false) const = 0;

  virtual void node_get_coordinates(Entity_ID v, AmanziGeometry::Point* xyz) const = 0;

 protected:
  unsigned int space_dim_;
};

}  // namespace AmanziMesh
}  // namespace Amanzi

#endif
