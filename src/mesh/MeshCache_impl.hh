/*
  Copyright 2010-201x held jointly by LANL, ORNL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
           Julien Loiseau (jloiseau@lanl.gov)
           Rao Garimella (rao@lanl.gov)
*/
//! Caches mesh information for fast repeated access.

#pragma once

#include "MeshDefs.hh"
#include "MeshAlgorithms.hh"
#include "MeshFramework.hh"
#include "MeshCache_decl.hh"

namespace Amanzi {
namespace AmanziMesh {

//
// This is an example where all specializations are implemented
//
template<AccessPattern AP>
AmanziGeometry::Point
MeshCache::getNodeCoordinate(const Entity_ID n) const
{
  if (node_coordinates_cached) return node_coordinates[n];
  if (framework_mesh_.get()) return framework_mesh_->getNodeCoordinate(n);
  throwAccessError_("getNodeCoordinate");
  return AmanziGeometry::Point();
}

template<> inline
AmanziGeometry::Point
MeshCache::getNodeCoordinate<AccessPattern::CACHE>(const Entity_ID n) const
{
  return node_coordinates[n];
}

template<> inline
AmanziGeometry::Point
MeshCache::getNodeCoordinate<AccessPattern::FRAMEWORK>(const Entity_ID n) const
{
  return framework_mesh_->getNodeCoordinate(n);
}

template<> inline
AmanziGeometry::Point
MeshCache::getNodeCoordinate<AccessPattern::RECOMPUTE>(const Entity_ID n) const
{
  throwAccessError_("getNodeCoordinate");
  return AmanziGeometry::Point();
}

//
// This is an example where only the default is currently implemented because
// there is no real need to use the framework in this case.  We COULD implement
// all of these in the future, and the interface would allow this.
//
template<AccessPattern AP>
Point_View
MeshCache::getEdgeCoordinates(const Entity_ID n) const
{
  if (edge_coordinates_cached) return edge_coordinates[n];
  return MeshAlgorithms::getEdgeCoordinates(*this, n);
}

template<AccessPattern AP>
Point_View
MeshCache::getFaceCoordinates(const Entity_ID n) const
{
  if (face_coordinates_cached) return face_coordinates[n];
  return MeshAlgorithms::getFaceCoordinates(*this, n);
}

template<AccessPattern AP>
Point_View
MeshCache::getCellCoordinates(const Entity_ID n) const
{
  if (cell_coordinates_cached) return cell_coordinates[n];
  return MeshAlgorithms::getCellCoordinates(*this, n);
}

template<AccessPattern AP>
AmanziGeometry::Point
MeshCache::getCellCentroid(const Entity_ID c) const
{
  if (cell_geometry_cached) return cell_centroids[c];
  auto vol_cent = MeshAlgorithms::computeCellGeometry(*this, c);
  return vol_cent.second;
}

template<AccessPattern AP>
AmanziGeometry::Point
MeshCache::getFaceCentroid(const Entity_ID f) const
{
  if (face_geometry_cached) return face_centroids[f];
  auto area_cent_norms = MeshAlgorithms::computeFaceGeometry(*this, f);
  return std::get<1>(area_cent_norms);
}

template<AccessPattern AP>
AmanziGeometry::Point
MeshCache::getEdgeCentroid(const Entity_ID e) const
{
  if (edge_geometry_cached) return edge_centroids[e];
  auto vec_cent = MeshAlgorithms::computeEdgeGeometry(*this, e);
  return vec_cent.second;
}

template<AccessPattern AP>
AmanziGeometry::Point
MeshCache::getCentroid(const Entity_kind kind, const Entity_ID ent)
{
  switch(kind) {
    case (Entity_kind::CELL) :
      return getCellCentroid<AP>(ent);
      break;
    case (Entity_kind::FACE) :
      return getFaceCentroid<AP>(ent);
      break;
    case (Entity_kind::EDGE) :
      return getEdgeCentroid<AP>(ent);
      break;
    case (Entity_kind::NODE) :
      return getNodeCoordinate<AP>(ent);
      break;
    default :
      Errors::Message msg("Invalid argument kind to getCentroid");
      Exceptions::amanzi_throw(msg);
  }
  return AmanziGeometry::Point();
}

// extent
template<AccessPattern AP>
double
MeshCache::getCellVolume(const Entity_ID c) const
{
  if (cell_geometry_cached) return cell_volumes[c];
  auto vol_cent = MeshAlgorithms::computeCellGeometry(*this, c);
  return vol_cent.first;
}

template<AccessPattern AP>
double
MeshCache::getFaceArea(const Entity_ID f) const
{
  if (face_geometry_cached) return AmanziGeometry::norm(face_normals[f][0]);
  auto area_cent_norms = MeshAlgorithms::computeFaceGeometry(*this, f);
  return std::get<0>(area_cent_norms);
}

template<AccessPattern AP>
double
MeshCache::getEdgeLength(const Entity_ID e) const
{
  if (edge_geometry_cached) return AmanziGeometry::norm(edge_vectors[e]);
  auto vec_cent = MeshAlgorithms::computeEdgeGeometry(*this, e);
  return AmanziGeometry::norm(vec_cent.first);
}

// Normal vector of a face
template<AccessPattern AP>
AmanziGeometry::Point
MeshCache::getFaceNormal(const Entity_ID f) const
{
  // if (face_geometry_cached) return face_normals[f][0];
  // if (framework_mesh_.get()) return framework_mesh_->getFaceNormal(f);
  // auto area_cent_norms = MeshAlgorithms::computeFaceGeometry(*this, f);
  // return std::get<2>(area_cent_norms)[0];
  // optimize this later! --etc
  return getFaceNormal(f, -1, nullptr);
}

// Normal vector and natural direction of a face, outward with respect to a
// cell.
//
// The vector is normalized and then weighted by the area of the face.
//
// The orientation is 1 if the outward normal is the same direction as the
// natural normal, -1 if in opposite directions, and 0 if there is no natural
// normal.
template<AccessPattern AP>
AmanziGeometry::Point
MeshCache::getFaceNormal(const Entity_ID f,
        const Entity_ID c, int* orientation) const
{
  AmanziGeometry::Point normal;

  if (face_geometry_cached) {
    auto fcells = getFaceCells(f, Parallel_type::ALL);
    if (orientation) *orientation = 0;
    Entity_ID cc = (c < 0) ? fcells[0] : c;

    int i = std::find(fcells.begin(), fcells.end(), cc) - fcells.begin();
    normal = face_normals[f][i];

    if (getSpaceDimension() == getManifoldDimension()) {
      if (c < 0) {
        normal *= MeshAlgorithms::getFaceDirectionInCell(*this, f, cc);
      } else if (orientation) {
        *orientation = MeshAlgorithms::getFaceDirectionInCell(*this, f, cc);
      }
    } else if (c < 0) {
      Errors::Message msg("MeshFramework: asking for the natural normal of a submanifold mesh is not valid.");
      Exceptions::amanzi_throw(msg);
    }

  } else if (framework_mesh_.get()) {
    normal = framework_mesh_->getFaceNormal(f, c, orientation);

  } else {
    auto geom = MeshAlgorithms::computeFaceGeometry(*this, f);

    auto fcells = getFaceCells(f, Parallel_type::ALL);
    if (orientation) *orientation = 0;
    Entity_ID cc = (c < 0) ? fcells[0] : c;
    int i = std::find(fcells.begin(), fcells.end(), cc) - fcells.begin();
    normal = std::get<2>(geom)[i];

    if (getSpaceDimension() == getManifoldDimension()) {
      if (c < 0) {
        normal *= MeshAlgorithms::getFaceDirectionInCell(*this, f, cc);
      } else if (orientation) {
        *orientation = MeshAlgorithms::getFaceDirectionInCell(*this, f, cc);
      }
    } else if (c < 0) {
      Errors::Message msg("MeshFramework: asking for the natural normal of a submanifold mesh is not valid.");
      Exceptions::amanzi_throw(msg);
    }
  }
  return normal;
}

// Vector describing the edge, where the length is the edge length.
//
// Orientation is the natural orientation, e.g. that it points from node 0 to
// node 1 with respect to edge_node adjacency information.
template<AccessPattern AP>
AmanziGeometry::Point
MeshCache::getEdgeVector(const Entity_ID e) const
{
  if (edge_geometry_cached) return edge_vectors[e];
  auto vec_cent = MeshAlgorithms::computeEdgeGeometry(*this, e);
  return vec_cent.first;
}

//---------------------
// Downward adjacencies
//---------------------
// Get faces of a cell
//
// On a distributed mesh, this will return all the faces of the
// cell, OWNED or GHOST. If the framework supports it, the faces will be
// returned in a standard order according to Exodus II convention
// for standard cells in all other situations (not supported or
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
template<AccessPattern AP>
std::size_t
MeshCache::getCellNumFaces(const Entity_ID c) const
{
  return getCellFaces(c).size();
}

template<AccessPattern AP>
const Entity_ID_View
MeshCache::getCellFaces(const Entity_ID c) const
{
  if (cell_geometry_cached) return cell_faces[c];
  if (framework_mesh_.get()) {
    Entity_ID_View cfaces;
    framework_mesh_->getCellFaces(c, cfaces);
    return cfaces;
  }
  throwAccessError_("getCellNumFaces");
  return cell_faces[c]; // NOTE: silences compiler warnings
}

template<AccessPattern AP>
const Entity_Direction_View
MeshCache::getCellFaceDirections(const Entity_ID c) const
{
  if (cell_geometry_cached) return cell_face_directions[c];
  if (framework_mesh_.get()) {
    Entity_ID_View faces;
    Entity_Direction_View dirs;
    framework_mesh_->getCellFacesAndDirs(c, faces, &dirs);
    return dirs;
  }
  throwAccessError_("getCellNumFaces");
  return cell_face_directions[c]; // NOTE: silences compiler warnings
}


template<AccessPattern AP>
std::pair<const Entity_ID_View, const Entity_Direction_View>
MeshCache::getCellFacesAndDirections(const Entity_ID c) const
{
  if (cell_geometry_cached) return std::make_pair(cell_faces[c], cell_face_directions[c]);
  if (framework_mesh_.get()) {
    Entity_ID_View faces;
    Entity_Direction_View dirs;
    framework_mesh_->getCellFacesAndDirs(c, faces, &dirs);
    return std::make_pair(faces, dirs);
  }
  throwAccessError_("getCellNumFaces");
  return std::make_pair(cell_faces[c], cell_face_directions[c]); // NOTE: silences compiler warnings
}

template<AccessPattern AP>
std::pair<const Entity_ID_View, const Point_View>
MeshCache::getCellFacesAndBisectors(const Entity_ID c) const
{
  if (face_geometry_cached) return std::make_pair(cell_faces[c], cell_face_bisectors[c]);
  Entity_ID_View faces = getCellFaces(c);
  Point_View bisectors;
  MeshAlgorithms::computeBisectors(*this, c, faces, bisectors);
  return std::make_pair(faces, bisectors);
}

// NOTE: all deprecated pragmas should go back in after finished refactoring
// [[deprecated("Prefer to use non-void variant that returns faces directly")]]
template<AccessPattern AP>
void
MeshCache::getCellFaces(const Entity_ID c,
                        Entity_ID_View& faces) const
{
  faces = getCellFaces(c);
}


//[[deprecated("Prefer to use non-void variant that returns faces directly")]]
template<AccessPattern AP>
void
MeshCache::getCellFacesAndDirs(const Entity_ID c,
        Entity_ID_View& faces,
        Entity_Direction_View * const dirs) const
{
  if (dirs) std::tie(faces, *dirs) = getCellFacesAndDirections(c);
  else faces = getCellFaces(c);
}

//[[deprecated("Prefer to use non-void variant that returns faces directly")]]
template<AccessPattern AP>
void
MeshCache::getCellFacesAndBisectors(
  const Entity_ID c,
  Entity_ID_View& faces,
  Point_View * const bisectors) const
{
  if (bisectors) std::tie(faces, *bisectors) = getCellFacesAndBisectors(c);
  else faces = getCellFaces(c);
}

// Get edges of a cell.
template<AccessPattern AP>
std::size_t
MeshCache::getCellNumEdges(const Entity_ID c) const
{
  return getCellEdges(c).size();
}

template<AccessPattern AP>
const Entity_ID_View
MeshCache::getCellEdges(const Entity_ID c) const
{
  if (cell_edges_cached) return cell_edges[c];
  return MeshAlgorithms::computeCellEdges(*this, c);
}

//[[deprecated("Prefer to use non-void variant that returns edges directly")]]
template<AccessPattern AP>
void
MeshCache::getCellEdges(const Entity_ID c, Entity_ID_View& edges) const
{
  edges = getCellEdges(c);
}

// Get nodes of a cell.
template<AccessPattern AP>
std::size_t
MeshCache::getCellNumNodes(const Entity_ID c) const
{
  return getCellNodes(c).size();
}

template<AccessPattern AP>
const Entity_ID_View
MeshCache::getCellNodes(const Entity_ID c) const
{
  if (cell_nodes_cached) return cell_nodes[c];
  if (framework_mesh_) {
    // note we prefer the framework to the algorithm, because the framework may
    // provide these in a special order, while the algorithm does not!
    //
    // Note this opens up to possible bugs in 2D -- todo -- order nodes in the
    // algorithm in 2D at least.
    Entity_ID_View cnodes;
    framework_mesh_->getCellNodes(c, cnodes);
    return cnodes;
  }
  return MeshAlgorithms::computeCellNodes(*this, c);
}

//[[deprecated("Prefer to use non-void variant that returns nodes directly")]]
template<AccessPattern AP>
void
MeshCache::getCellNodes(const Entity_ID c, Entity_ID_View& nodes) const
{
  nodes = getCellNodes(c);
}

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
template<AccessPattern AP>
std::size_t
MeshCache::getFaceNumEdges(const Entity_ID f) const
{
  return getFaceEdges(f).size();
}

template<AccessPattern AP>
const Entity_ID_View
MeshCache::getFaceEdges(const Entity_ID f) const
{
  if (face_edges_cached) return face_edges[f];
  if (framework_mesh_.get()) {
    Entity_ID_List edges;
    framework_mesh_->getFaceEdges(f, edges);
    return edges;
  }
  throwAccessError_("getFaceEdges");
  return face_edges[f]; // silences warnings
}

template<AccessPattern AP>
std::pair<const Entity_ID_View, const Entity_Direction_View>
MeshCache::getFaceEdgesAndDirections(const Entity_ID f) const
{
  if (face_edges_cached) return std::make_pair(face_edges[f], face_edge_directions[f]);
  if (framework_mesh_.get()) {
    Entity_ID_View edges;
    Entity_Direction_View dirs;
    framework_mesh_->getFaceEdgesAndDirs(f, edges, &dirs);
    return std::make_pair(edges, dirs);
  }
  throwAccessError_("getFaceEdgesAndDirections");
  return std::make_pair(Entity_ID_View(), Entity_Direction_View());
}

//[[deprecated("Prefer to use non-void variant that returns edges directly")]]
template<AccessPattern AP>
void
MeshCache::getFaceEdges(const Entity_ID f,
                        Entity_ID_View& edges) const
{
  edges = getFaceEdges(f);
}

//[[deprecated("Prefer to use non-void variant that returns edges directly")]]
template<AccessPattern AP>
void
MeshCache::getFaceEdgesAndDirs(const Entity_ID f,
        Entity_ID_View& edges,
        Entity_Direction_View * const dirs) const
{
  if (dirs) std::tie(edges, *dirs) = getFaceEdgesAndDirections(f);
  else edges = getFaceEdges(f);
}

// Get nodes of face
//
// In 3D, the nodes of the face are returned in ccw order consistent
// with the face normal.
template<AccessPattern AP>
std::size_t
MeshCache::getFaceNumNodes(const Entity_ID f) const
{
  return getFaceNodes(f).size();
}

template<AccessPattern AP>
const Entity_ID_View
MeshCache::getFaceNodes(const Entity_ID f) const
{
  if (face_nodes_cached) return face_nodes[f];
  if (framework_mesh_.get()) {
    Entity_ID_View fnodes;
    framework_mesh_->getFaceNodes(f, fnodes);
    return fnodes;
  }
  throwAccessError_("getFaceNodes");
  return face_nodes[f]; // silence warnings
}

//[[deprecated("Prefer to use non-void variant that returns nodes directly")]]
template<AccessPattern AP>
void
MeshCache::getFaceNodes(const Entity_ID f, Entity_ID_View& nodes) const
{
  nodes = getFaceNodes(f);
}

// Get nodes of edge
template<AccessPattern AP>
std::size_t
MeshCache::getEdgeNumNodes(const Entity_ID e) const
{
  return getEdgeNodes(e).size();
}

template<AccessPattern AP>
const Entity_ID_View
MeshCache::getEdgeNodes(const Entity_ID e) const
{
  if (edge_nodes_cached) return edge_nodes[e];
  if (framework_mesh_.get()) {
    Entity_ID_View enodes;
    framework_mesh_->getEdgeNodes(e, enodes);
    return enodes;
  }
  throwAccessError_("getEdgeNodes");
  return edge_nodes[e];
}

//[[deprecated("Prefer to use non-void variant that returns nodes directly")]]
template<AccessPattern AP>
void
MeshCache::getEdgeNodes(const Entity_ID e, Entity_ID_View& nodes) const
{
  nodes = getEdgeNodes(e);
}

//[[deprecated("Prefer to use non-void variant that returns nodes directly")]]
template<AccessPattern AP>
void
MeshCache::getEdgeNodes(const Entity_ID e, Entity_ID* n0, Entity_ID* n1) const
{
  auto nodes = getEdgeNodes(e);
  if (n0) *n0 = nodes[0];
  if (n1) *n1 = nodes[1];
}


//-------------------
// Upward adjacencies
//-------------------
// The cells are returned in no particular order. Also, the order of cells
// is not guaranteed to be the same for corresponding faces on different
// processors
template<AccessPattern AP>
std::size_t
MeshCache::getFaceNumCells(const Entity_ID f, const Parallel_type ptype) const
{
  return getFaceCells(f, ptype).size();
}

template<AccessPattern AP>
const Entity_ID_View
MeshCache::getFaceCells(const Entity_ID f, const Parallel_type ptype) const
{
  if (face_cells_cached) return face_cells[f];
  if (framework_mesh_.get()) {
    Entity_ID_List fcells;
    framework_mesh_->getFaceCells(f, ptype, fcells);
    return fcells;
  }
  throwAccessError_("getFaceCells");
  return face_cells[f]; // silences warnings
}

template<AccessPattern AP>
void
MeshCache::getFaceCells(const Entity_ID f,
                        const Parallel_type ptype,
                        Entity_ID_View& cells) const
{
  cells = getFaceCells(f, ptype);
}



} // namespace AmanziMesh
} // namespace Amanzi
