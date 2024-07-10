/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
           Julien Loiseau (jloiseau@lanl.gov)
           Rao Garimella (rao@lanl.gov)
*/

//! Caches mesh information for fast repeated access.

#pragma once

#include "Kokkos_Core.hpp"
#include "MeshDefs.hh"
#include "MeshCacheData.hh"
#include "MeshColumns.hh"

namespace Amanzi {
namespace AmanziMesh {

struct MeshCache {
  // view types
  static const MemSpace_kind MEM = MemSpace_kind::DEVICE;
  using Entity_ID_View = View_type<Entity_ID, MEM>;
  using cEntity_ID_View = View_type<const Entity_ID, MEM>;
  using Entity_GID_View = View_type<Entity_GID, MEM>;
  using cEntity_GID_View = View_type<const Entity_GID, MEM>;
  using Direction_View = View_type<Direction_type, MEM>;
  using cDirection_View = View_type<const Direction_type, MEM>;
  using Point_View = View_type<AmanziGeometry::Point, MEM>;
  using cPoint_View = View_type<const AmanziGeometry::Point, MEM>;
  using Double_View = View_type<double, MEM>;
  using cDouble_View = View_type<const double, MEM>;

  MeshCache() {}
  MeshCache(const MeshCacheData& data_) : data(data_) {}

  // ----------------
  // mesh properties
  // ----------------
  KOKKOS_INLINE_FUNCTION bool isOrdered() const { return data.is_ordered_; }
  KOKKOS_INLINE_FUNCTION bool isLogical() const { return data.is_logical_; }
  KOKKOS_INLINE_FUNCTION bool isSFM() const { return data.is_sfm_; } // single face mesh -- special case

  KOKKOS_INLINE_FUNCTION bool hasNodes() const { return data.has_nodes_; }
  KOKKOS_INLINE_FUNCTION bool hasEdges() const { return data.has_edges_; }
  KOKKOS_INLINE_FUNCTION bool hasNodeFaces() const { return data.has_node_faces_; }

  KOKKOS_INLINE_FUNCTION int getSpaceDimension() const { return data.space_dim_; }
  KOKKOS_INLINE_FUNCTION int getManifoldDimension() const { return data.manifold_dim_; }

  // ----------------
  // Entity meta-data
  // ----------------
  KOKKOS_INLINE_FUNCTION
  Entity_ID getNumEntities(const Entity_kind kind, const Parallel_kind ptype) const;

  KOKKOS_INLINE_FUNCTION
  Entity_ID getEntityParent(const Entity_kind kind, const Entity_ID entid) const;

  KOKKOS_INLINE_FUNCTION
  Parallel_kind getParallelKind(const Entity_kind kind, const Entity_ID id) const;

  KOKKOS_INLINE_FUNCTION
  Entity_ID getBoundaryFaceFace(const Entity_ID bf) const;

  KOKKOS_INLINE_FUNCTION
  Entity_ID getBoundaryNodeNode(const Entity_ID bf) const;

  //---------------------
  // Geometry
  //---------------------
  // node locations
  KOKKOS_INLINE_FUNCTION
  AmanziGeometry::Point getNodeCoordinate(const Entity_ID n) const;

  // coordinate views
  template <AccessPattern_kind = AccessPattern_kind::DEFAULT>
  KOKKOS_INLINE_FUNCTION cPoint_View getEdgeCoordinates(const Entity_ID e) const;

  template <AccessPattern_kind = AccessPattern_kind::DEFAULT>
  KOKKOS_INLINE_FUNCTION cPoint_View getFaceCoordinates(const Entity_ID f) const;

  template <AccessPattern_kind = AccessPattern_kind::DEFAULT>
  KOKKOS_INLINE_FUNCTION cPoint_View getCellCoordinates(const Entity_ID c) const;

  // centroids
  template <AccessPattern_kind AP = AccessPattern_kind::DEFAULT>
  KOKKOS_INLINE_FUNCTION AmanziGeometry::Point getCellCentroid(const Entity_ID c) const;

  template <AccessPattern_kind = AccessPattern_kind::DEFAULT>
  KOKKOS_INLINE_FUNCTION AmanziGeometry::Point getFaceCentroid(const Entity_ID f) const;

  template <AccessPattern_kind = AccessPattern_kind::DEFAULT>
  KOKKOS_INLINE_FUNCTION AmanziGeometry::Point getEdgeCentroid(const Entity_ID e) const;

  // at run-time, calls one of get*Centroid
  template <AccessPattern_kind = AccessPattern_kind::DEFAULT>
  KOKKOS_INLINE_FUNCTION AmanziGeometry::Point
  getCentroid(const Entity_kind kind, const Entity_ID ent) const;

  // at compile-time, calls one of get*Centroid
  template <Entity_kind, AccessPattern_kind = AccessPattern_kind::DEFAULT>
  KOKKOS_INLINE_FUNCTION AmanziGeometry::Point getCentroid(const Entity_ID ent) const;

  // extent
  template <AccessPattern_kind AP = AccessPattern_kind::DEFAULT>
  KOKKOS_INLINE_FUNCTION double getCellVolume(const Entity_ID c) const;

  template <AccessPattern_kind AP = AccessPattern_kind::DEFAULT>
  KOKKOS_INLINE_FUNCTION double getFaceArea(const Entity_ID f) const;

  template <AccessPattern_kind = AccessPattern_kind::DEFAULT>
  KOKKOS_INLINE_FUNCTION double getEdgeLength(const Entity_ID e) const;

  // at run-time, calls one of getCellVolume, getFaceArea, or getEdgeLength
  template <Entity_kind EK, AccessPattern_kind AP = AccessPattern_kind::DEFAULT>
  KOKKOS_INLINE_FUNCTION double getExtent(const Entity_ID e) const;

  // at compile-time, calls one of getCellVolume, getFaceArea, or getEdgeLength
  template <AccessPattern_kind AP = AccessPattern_kind::DEFAULT>
  KOKKOS_INLINE_FUNCTION double getExtent(const Entity_kind kind, const Entity_ID e) const;

  // Normal vector of a face
  template <AccessPattern_kind = AccessPattern_kind::DEFAULT>
  KOKKOS_INLINE_FUNCTION AmanziGeometry::Point getFaceNormal(const Entity_ID f) const;

  // Normal vector and natural direction of a face, outward with respect to a
  // cell.
  //
  // The vector is normalized and then weighted by the area of the face.
  //
  // The orientation is 1 if the outward normal is the same direction as the
  // natural normal, -1 if in opposite directions, and 0 if there is no natural
  // normal.
  template <AccessPattern_kind = AccessPattern_kind::DEFAULT>
  KOKKOS_INLINE_FUNCTION AmanziGeometry::Point
  getFaceNormal(const Entity_ID f, const Entity_ID c, int* orientation = nullptr) const;

  // Vector describing the edge, where the length is the edge length.
  //
  // Orientation is the natural orientation, e.g. that it points from node 0 to
  // node 1 with respect to edge_node adjacency information.
  template <AccessPattern_kind = AccessPattern_kind::DEFAULT>
  KOKKOS_INLINE_FUNCTION AmanziGeometry::Point getEdgeVector(const Entity_ID e) const;

  //---------------------
  // Downward adjacencies
  //---------------------
  // Get faces of a cell
  //
  // On a distributed mesh, this will return all the faces of the
  // cell, OWNED or GHOST. If the framework supports it, the faces will be
  // returned in a standard order according to Exodus II convention
  // for standard cells; in all other situations (not supported or
  // non-standard cells), the list of faces will be in arbitrary order
  //
  // In 3D manifolds and 2D manifolds of 2d space, the natural direction of a
  // face is defined by the right-hand-rule of the node ordering.  In 2D
  // manifolds of 3D space, there is no natural direction.
  //
  // Use of bisectors instead of coordinate geometry enables use of logical
  // meshes which may not have geometric coordinates.
  //
  // New interfaces that return by const reference should be preferred rather
  // than those that return void and expect the return value as an argument --
  // this new-style interface works better with Kokkos and should be more
  // efficient in all cases.
  KOKKOS_INLINE_FUNCTION size_type getCellNumFaces(const Entity_ID c) const;


  // note, no AccessPattern_kind -- only works on cached
  KOKKOS_INLINE_FUNCTION
  const Entity_ID& getCellFace(const Entity_ID c, const size_type i) const;

  KOKKOS_INLINE_FUNCTION
  Kokkos::pair<cEntity_ID_View, cDirection_View> getCellFacesAndDirections(const Entity_ID c) const;

  KOKKOS_INLINE_FUNCTION
  Kokkos::pair<cEntity_ID_View, cPoint_View> getCellFacesAndBisectors(const Entity_ID c) const;

  KOKKOS_INLINE_FUNCTION
  cEntity_ID_View getCellFaces(const Entity_ID c) const;

  //
  // Downward adjacency
  //
  // Get edges of a cell
  KOKKOS_INLINE_FUNCTION size_type getCellNumEdges(const Entity_ID c) const;

  KOKKOS_INLINE_FUNCTION
  const Entity_ID& getCellEdge(const Entity_ID c, const size_type i) const;

  template <AccessPattern_kind AP = AccessPattern_kind::DEFAULT>
  KOKKOS_INLINE_FUNCTION cEntity_ID_View getCellEdges(const Entity_ID c) const;

  // Get nodes of a cell.
  KOKKOS_INLINE_FUNCTION size_type getCellNumNodes(const Entity_ID c) const;

  KOKKOS_INLINE_FUNCTION
  Entity_ID getCellNode(const Entity_ID c, const size_type i) const;

  template <AccessPattern_kind = AccessPattern_kind::DEFAULT>
  KOKKOS_INLINE_FUNCTION
  cEntity_ID_View getCellNodes(const Entity_ID c) const;


  // Get edges of a face and directions in which the face uses the edges.
  //
  // In 3D, edge direction is 1 when it is oriented counter clockwise
  // with respect to the face natural normal.
  //
  // On a distributed mesh, this will return all the edges of the
  // face, OWNED or GHOST. If the framework supports it, the edges will be
  // returned in a ccw order around the face as it is naturally defined.
  //
  // IMPORTANT NOTE IN 2D: In meshes where the cells are two
  // dimensional, faces and edges are identical. For such cells, this
  // operator will return a single edge and a direction of 1. However,
  // this direction cannot be relied upon to compute, say, a contour
  // integral around the 2D cell.
  // template<AccessPattern_kind = AccessPattern_kind::DEFAULT>
  KOKKOS_INLINE_FUNCTION
  size_type getFaceNumEdges(const Entity_ID f) const;

  KOKKOS_INLINE_FUNCTION
  const Entity_ID& getFaceEdge(const Entity_ID f, const size_type i) const;

  KOKKOS_INLINE_FUNCTION
  cEntity_ID_View getFaceEdges(const Entity_ID f) const;

  KOKKOS_INLINE_FUNCTION
  Kokkos::pair<cEntity_ID_View, cDirection_View> getFaceEdgesAndDirections(const Entity_ID f) const;


  // Get nodes of face
  //
  // In 3D, the nodes of the face are returned in ccw order consistent
  // with the face normal.
  KOKKOS_INLINE_FUNCTION
  size_type getFaceNumNodes(const Entity_ID f) const;

  KOKKOS_INLINE_FUNCTION
  const Entity_ID& getFaceNode(const Entity_ID f, const size_type i) const;

  KOKKOS_INLINE_FUNCTION
  cEntity_ID_View getFaceNodes(const Entity_ID f) const;

  // Get nodes of edge
  KOKKOS_INLINE_FUNCTION
  size_type getEdgeNumNodes(const Entity_ID e) const;

  KOKKOS_INLINE_FUNCTION
  const Entity_ID& getEdgeNode(const Entity_ID e, const size_type i) const;

  KOKKOS_INLINE_FUNCTION
  cEntity_ID_View getEdgeNodes(const Entity_ID e) const;


  //-------------------
  // Upward adjacencies
  //-------------------
  // Face --> Cell
  // The cells are returned in no particular order. Also, the order of cells
  // is not guaranteed to be the same for corresponding faces on different
  // processors
  KOKKOS_INLINE_FUNCTION
  size_type getFaceNumCells(const Entity_ID f) const;

  KOKKOS_INLINE_FUNCTION
  const Entity_ID& getFaceCell(const Entity_ID f, const size_type i) const;

  KOKKOS_INLINE_FUNCTION
  cEntity_ID_View getFaceCells(const Entity_ID f) const;

  // Edge --> Cell
  //
  // The order of cells is not guaranteed to be the same for corresponding
  // edges on different processors
  KOKKOS_INLINE_FUNCTION
  size_type getEdgeNumCells(const Entity_ID f) const;

  KOKKOS_INLINE_FUNCTION
  const Entity_ID& getEdgeCell(const Entity_ID f, const size_type i) const;

  template <AccessPattern_kind AP = AccessPattern_kind::DEFAULT>
  KOKKOS_INLINE_FUNCTION
  cEntity_ID_View getEdgeCells(const Entity_ID e) const;

  // Edge --> Face
  // NOTE: The order of faces is not guaranteed to be the same for
  // corresponding edges on different processors
  KOKKOS_INLINE_FUNCTION
  size_type getEdgeNumFaces(const Entity_ID f) const;

  KOKKOS_INLINE_FUNCTION
  const Entity_ID& getEdgeFace(const Entity_ID f, const size_type i) const;

  KOKKOS_INLINE_FUNCTION
  cEntity_ID_View getEdgeFaces(const Entity_ID e) const;

  // Node --> Cell
  // NOTE: The order of cells is not guaranteed to be the same for
  // corresponding nodes on different processors
  template <AccessPattern_kind AP = AccessPattern_kind::DEFAULT>
  KOKKOS_INLINE_FUNCTION cEntity_ID_View
  getNodeCells(const Entity_ID n) const;

  // Node --> Face
  // Faces of type parallel 'ptype' connected to a node
  // NOTE: The order of faces is not guarnateed to be the same for
  // corresponding nodes on different processors
  KOKKOS_INLINE_FUNCTION cEntity_ID_View getNodeFaces(const Entity_ID n) const;

  // Node --> Edge
  //
  // The order of edges is not guaranteed to be the same for corresponding
  // node on different processors
  KOKKOS_INLINE_FUNCTION cEntity_ID_View getNodeEdges(const Entity_ID n) const;

  // all the dual views
  MeshCacheData data;

  // columnar structure dual views
  MeshColumns columns;
};


} // namespace AmanziMesh
} // namespace Amanzi
