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
// See MeshCache_decl.hh for documentation!

#pragma once

#include "MeshCache_decl.hh"
#include "MeshHelpers.hh"

namespace Amanzi {
namespace AmanziMesh {

// ----------------
// Entity meta-data
// ----------------
KOKKOS_INLINE_FUNCTION
Entity_ID MeshCache::getNumEntities(const Entity_kind kind, const Parallel_kind ptype) const
{
  Entity_ID nowned, nall;
  switch (kind) {
    case (Entity_kind::CELL):
      nowned = data.ncells_owned;
      nall = data.ncells_all;
      break;
    case (Entity_kind::FACE):
      nowned = data.nfaces_owned;
      nall = data.nfaces_all;
      break;
    case (Entity_kind::EDGE):
      nowned = data.nedges_owned;
      nall = data.nedges_all;
      break;
    case (Entity_kind::NODE):
      nowned = data.nnodes_owned;
      nall = data.nnodes_all;
      break;
    case (Entity_kind::BOUNDARY_FACE):
      nowned = data.nboundary_faces_owned;
      nall = data.nboundary_faces_all;
      break;
    case (Entity_kind::BOUNDARY_NODE):
      nowned = data.nboundary_nodes_owned;
      nall = data.nboundary_nodes_all;
      break;
    default:
      nowned = -1;
      nall = -1;
  }

  switch (ptype) {
    case (Parallel_kind::OWNED):
      return nowned;
      break;
    case (Parallel_kind::ALL):
      return nall;
      break;
    case Parallel_kind::GHOST:
      return nall - nowned;
      break;
    default:
      return 0;
  }
}


// corresponding entity in the parent mesh
//
// Note the kind refers to the kind in _this_ mesh -- for some lifted meshes,
// this may not be the same as the entity kind in the parent mesh.  That
// logic is left to the user of this class -- we simply store the IDs.
KOKKOS_INLINE_FUNCTION
Entity_ID MeshCache::getEntityParent(const Entity_kind kind, const Entity_ID entid) const
{
  if (data.parent_entities_cached) {
    switch (kind) {
      case Entity_kind::CELL:
        return view<MEM>(data.parent_cells)[entid];
        break;
      case Entity_kind::FACE:
        return view<MEM>(data.parent_faces)[entid];
        break;
      case Entity_kind::EDGE:
        return view<MEM>(data.parent_edges)[entid];
        break;
      case Entity_kind::NODE:
        return view<MEM>(data.parent_nodes)[entid];
      default: {
      }
    }
  }
  assert(false);
  return -1;
}


KOKKOS_INLINE_FUNCTION
Parallel_kind MeshCache::getParallelKind(const Entity_kind kind, const Entity_ID id) const
{
  if (id < getNumEntities(kind, Parallel_kind::OWNED)) {
    return Parallel_kind::OWNED;
  } else if (id < getNumEntities(kind, Parallel_kind::ALL)) {
    return Parallel_kind::GHOST;
  }
  return Parallel_kind::UNKNOWN;
}


KOKKOS_INLINE_FUNCTION
Entity_ID
MeshCache::getBoundaryFaceFace(const Entity_ID bf) const
{
  return view<MEM>(data.boundary_faces)(bf);
}

KOKKOS_INLINE_FUNCTION
Entity_ID
MeshCache::getBoundaryNodeNode(const Entity_ID bf) const
{
  return view<MEM>(data.boundary_nodes)(bf);
}


KOKKOS_INLINE_FUNCTION AmanziGeometry::Point
MeshCache::getNodeCoordinate(const Entity_ID n) const
{
  return Impl::Getter<MEM, AccessPattern_kind::CACHE>::get(
    data.node_coordinates_cached, data.node_coordinates, nullptr, n);
}


template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION typename MeshCache::cPoint_View
MeshCache::getEdgeCoordinates(const Entity_ID c) const
{
  return Impl::RaggedGetter<MEM, AP>::get(
    data.edge_coordinates_cached, data.edge_coordinates, nullptr, c);
}


template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION typename MeshCache::cPoint_View
MeshCache::getFaceCoordinates(const Entity_ID f) const
{
  return Impl::RaggedGetter<MEM, AP>::get(
    data.face_coordinates_cached, data.face_coordinates, nullptr, f);
}


template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION typename MeshCache::cPoint_View
MeshCache::getCellCoordinates(const Entity_ID c) const
{
  return Impl::RaggedGetter<MEM, AP>::get(
    data.cell_coordinates_cached, data.cell_coordinates, nullptr, c);
}


template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION AmanziGeometry::Point
MeshCache::getCellCentroid(const Entity_ID c) const
{
  return Impl::Getter<MEM, AP>::get(data.cell_geometry_cached, data.cell_centroids, nullptr, c);
}


template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION AmanziGeometry::Point
MeshCache::getFaceCentroid(const Entity_ID f) const
{
  return Impl::Getter<MEM, AP>::get(data.face_geometry_cached, data.face_centroids, nullptr, f);
}


template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION AmanziGeometry::Point
MeshCache::getEdgeCentroid(const Entity_ID c) const
{
  return Impl::Getter<MEM, AP>::get(data.edge_geometry_cached, data.edge_centroids, nullptr, c);
}


template <AccessPattern_kind AP>
AmanziGeometry::Point
MeshCache::getCentroid(const Entity_kind kind, const Entity_ID ent) const
{
  switch (kind) {
  case (Entity_kind::CELL):
    return getCellCentroid<AP>(ent);
  case (Entity_kind::FACE):
    return getFaceCentroid<AP>(ent);
  case (Entity_kind::BOUNDARY_FACE):
    return getFaceCentroid<AP>(getBoundaryFaceFace(ent));
  case (Entity_kind::EDGE):
    return getEdgeCentroid<AP>(ent);
  case (Entity_kind::NODE):
    return getNodeCoordinate(ent);
  default:
    assert(false);
  }
  return AmanziGeometry::Point();
}


template <Entity_kind EK, AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION AmanziGeometry::Point
MeshCache::getCentroid(const Entity_ID ent) const
{
  static_assert(EK != Entity_kind::UNKNOWN, "No centroid information for Entity_kind::UNKNOWN");
  if constexpr (EK == Entity_kind::CELL) {
    return getCellCentroid<AP>(ent);
  } else if constexpr (EK == Entity_kind::FACE) {
    return getFaceCentroid<AP>(ent);
  } else if constexpr (EK == Entity_kind::BOUNDARY_FACE) {
    return getFaceCentroid<AP>(getBoundaryFaceFace(ent));
  } else if constexpr (EK == Entity_kind::EDGE) {
    return getEdgeCentroid<AP>(ent);
  } else if constexpr (EK == Entity_kind::NODE) {
    return getNodeCoordinate(ent);
  }
  assert(false);
  return AmanziGeometry::Point();
}


template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION double
MeshCache::getCellVolume(const Entity_ID c) const
{
  return Impl::Getter<MEM, AP>::get(data.cell_geometry_cached, data.cell_volumes, nullptr, c);
}


template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION double
MeshCache::getFaceArea(const Entity_ID f) const
{
  return Impl::Getter<MEM, AP>::get(data.face_geometry_cached, data.face_areas, nullptr, f);
}


template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION double
MeshCache::getEdgeLength(const Entity_ID e) const
{
  return Impl::Getter<MEM, AP>::get(data.edge_lengths_cached, data.edge_lengths, nullptr, e);
}


template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION double
MeshCache::getExtent(const Entity_kind kind, const Entity_ID ent) const
{
  switch (kind) {
  case (Entity_kind::CELL):
    return getCellVolume<AP>(ent);
  case (Entity_kind::FACE):
    return getFaceArea<AP>(ent);
  case (Entity_kind::BOUNDARY_FACE):
    return getFaceArea<AP>(getBoundaryFaceFace(ent));
  case (Entity_kind::EDGE):
    return getEdgeLength<AP>(ent);
  default:
    assert(false);
  }
  return -1.0;
}


template <Entity_kind EK, AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION double
MeshCache::getExtent(const Entity_ID ent) const
{
  static_assert(EK == Entity_kind::CELL ||
                EK == Entity_kind::FACE ||
                EK == Entity_kind::BOUNDARY_FACE ||
                EK == Entity_kind::EDGE,
                "getExtent only supports CELL, FACE and EDGE");
  if constexpr (EK == Entity_kind::CELL) {
    return getCellVolume<AP>(ent);
  } else if constexpr (EK == Entity_kind::FACE) {
    return getFaceArea<AP>(ent);
  } else if constexpr (EK == Entity_kind::BOUNDARY_FACE) {
    return getFaceArea<AP>(getBoundaryFaceFace(ent));
  } else if constexpr (EK == Entity_kind::EDGE) {
    return getEdgeLength<AP>(ent);
  }
  return -1;
}


// Normal vector of a face
template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION AmanziGeometry::Point
MeshCache::getFaceNormal(const Entity_ID f) const
{
  return getFaceNormal<AP>(f, -1, nullptr);
}


template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION
AmanziGeometry::Point
MeshCache::getFaceNormal(const Entity_ID f, const Entity_ID c, int* orientation) const
{
  if (this->isSFM()) {
    AmanziGeometry::Point normal = data.face_normals.get<MemSpace_kind::HOST>(f, 0);
    return normal;
  }

  if constexpr (MEM == MemSpace_kind::DEVICE) { assert(data.face_geometry_cached); }

  if (data.face_geometry_cached) {
    if (c < 0) {
      assert(orientation == nullptr);
      auto normal = data.face_normals.get<MEM>(f, 0);
      auto nat_dir = data.face_normal_orientations.get<MEM>(f, 0);
      if (abs(nat_dir) == 1) {
        return nat_dir * normal;
      } else {
        auto fcells = getFaceCells(f);
        assert(fcells.size() == 2);
        assert(data.cell_global_indices_cached);

        // average normals oriented from lower to higher GIDs
        int pos_i = view<MemSpace_kind::DEVICE>(data.cell_global_indices)(fcells[0]) >
                        view<MemSpace_kind::HOST>(data.cell_global_indices)(fcells[1]) ?
                      0 :
                      1;
        return (data.face_normals.get<MEM>(f, 1 - pos_i) - data.face_normals.get<MEM>(f, pos_i)) /
               2;
      }
    }

    auto fcells = getFaceCells(f);
    std::size_t i;
    for (i = 0; i != fcells.size(); ++i)
      if (fcells[i] == c) break;

    if (orientation) *orientation = data.face_normal_orientations.get<MEM>(f, i) > 0 ? 1 : -1;
    return data.face_normals.get<MEM>(f, i);
  }

  assert(false && "No access to cache/framework/compute available");
  return AmanziGeometry::Point();
}


//---------------------
// Downward adjacencies
//---------------------
//
// Cell --> Face
//
KOKKOS_INLINE_FUNCTION
size_type
MeshCache::getCellNumFaces(const Entity_ID c) const
{
  assert(data.cell_faces_cached);
  return data.cell_faces.size<MEM>(c);
}


KOKKOS_INLINE_FUNCTION
const Entity_ID&
MeshCache::getCellFace(const Entity_ID c, const size_type i) const
{
  assert(data.cell_faces_cached);
  return data.cell_faces.get<MEM>(c, i);
}


KOKKOS_INLINE_FUNCTION
MeshCache::cEntity_ID_View
MeshCache::getCellFaces(const Entity_ID c) const
{
  assert(data.cell_faces_cached);
  return data.cell_faces.getRow<MEM>(c);
}


KOKKOS_INLINE_FUNCTION
Kokkos::pair<MeshCache::cEntity_ID_View,
             MeshCache::cDirection_View>
MeshCache::getCellFacesAndDirections(const Entity_ID c) const
{
  assert(data.cell_faces_cached);
  return Kokkos::make_pair(data.cell_faces.getRow<MEM>(c),
                           data.cell_face_directions.getRow<MEM>(c));
}


KOKKOS_INLINE_FUNCTION
Kokkos::pair<MeshCache::cEntity_ID_View,
             MeshCache::cPoint_View>
MeshCache::getCellFacesAndBisectors(const Entity_ID c) const
{
  assert(data.cell_faces_cached);
  return Kokkos::make_pair(data.cell_faces.getRow<MEM>(c),
                           data.cell_face_bisectors.getRow<MEM>(c));
}


//
// Cell --> Edge
//
KOKKOS_INLINE_FUNCTION
size_type
MeshCache::getCellNumEdges(const Entity_ID c) const
{
  assert(data.cell_edges_cached);
  return data.cell_edges.size<MEM>(c);
}


KOKKOS_INLINE_FUNCTION
const Entity_ID&
MeshCache::getCellEdge(const Entity_ID c, const size_type i) const
{
  assert(data.cell_edges_cached);
  return data.cell_edges.get<MEM>(c, i);
}


template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION
MeshCache::cEntity_ID_View
MeshCache::getCellEdges(const Entity_ID c) const
{
  return Impl::RaggedGetter<MEM, AP>::get(data.cell_edges_cached, data.cell_edges, nullptr, c);
}


//
// Cell --> Node
//
KOKKOS_INLINE_FUNCTION
size_type
MeshCache::getCellNumNodes(const Entity_ID c) const
{
  assert(data.cell_nodes_cached);
  return data.cell_nodes.size<MEM>(c);
}


KOKKOS_INLINE_FUNCTION
Entity_ID
MeshCache::getCellNode(const Entity_ID c, const size_type i) const
{
  assert(data.cell_nodes_cached);
  return data.cell_nodes.get<MEM>(c,i);
}


template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION
MeshCache::cEntity_ID_View
MeshCache::getCellNodes(const Entity_ID c) const
{
  return Impl::RaggedGetter<MEM, AP>::get(data.cell_nodes_cached, data.cell_nodes, nullptr, c);
}


//
// Face --> Edge
//
KOKKOS_INLINE_FUNCTION
size_type
MeshCache::getFaceNumEdges(const Entity_ID f) const
{
  assert(data.face_edges_cached);
  return data.face_edges.size<MEM>(f);
}


KOKKOS_INLINE_FUNCTION
const Entity_ID&
MeshCache::getFaceEdge(const Entity_ID f, const size_type i) const
{
  assert(data.face_edges_cached);
  return data.face_edges.get<MEM>(f,i);
}


KOKKOS_INLINE_FUNCTION
MeshCache::cEntity_ID_View
MeshCache::getFaceEdges(const Entity_ID f) const
{
  assert(data.face_edges_cached);
  return data.face_edges.getRow<MEM>(f);
}


KOKKOS_INLINE_FUNCTION
Kokkos::pair<MeshCache::cEntity_ID_View, MeshCache::cDirection_View>
MeshCache::getFaceEdgesAndDirections(const Entity_ID f) const
{
  assert(data.face_edges_cached);
  return Kokkos::make_pair(data.face_edges.getRow<MEM>(f),
                           data.face_edge_directions.getRow<MEM>(f));
}


//
// Face --> Node
//
KOKKOS_INLINE_FUNCTION
size_type
MeshCache::getFaceNumNodes(const Entity_ID f) const
{
  assert(data.face_nodes_cached);
  return data.face_nodes.size<MEM>(f);
}


KOKKOS_INLINE_FUNCTION
const Entity_ID&
MeshCache::getFaceNode(const Entity_ID f, const size_type i) const
{
  assert(data.face_nodes_cached);
  return data.face_nodes.get<MEM>(f,i);
}


KOKKOS_INLINE_FUNCTION
MeshCache::cEntity_ID_View
MeshCache::getFaceNodes(const Entity_ID f) const
{
  assert(data.face_nodes_cached);
  return data.face_nodes.getRow<MEM>(f);
}


//
// Edge --> Node
//
KOKKOS_INLINE_FUNCTION
size_type
MeshCache::getEdgeNumNodes(const Entity_ID f) const
{
  assert(data.edge_nodes_cached);
  return data.edge_nodes.size<MEM>(f);
}


KOKKOS_INLINE_FUNCTION
const Entity_ID&
MeshCache::getEdgeNode(const Entity_ID f, const size_type i) const
{
  assert(data.edge_nodes_cached);
  return data.edge_nodes.get<MEM>(f,i);
}


KOKKOS_INLINE_FUNCTION
MeshCache::cEntity_ID_View
MeshCache::getEdgeNodes(const Entity_ID f) const
{
  assert(data.edge_nodes_cached);
  return data.edge_nodes.getRow<MEM>(f);
}


//
// Face --> Cell
//
KOKKOS_INLINE_FUNCTION
size_type
MeshCache::getFaceNumCells(const Entity_ID f) const
{
  assert(data.face_cells_cached);
  return data.face_cells.size<MEM>(f);
}


KOKKOS_INLINE_FUNCTION
const Entity_ID&
MeshCache::getFaceCell(const Entity_ID f, const size_type i) const
{
  assert(data.face_cells_cached);
  return data.face_cells.get<MEM>(f,i);
}


KOKKOS_INLINE_FUNCTION
MeshCache::cEntity_ID_View
MeshCache::getFaceCells(const Entity_ID f) const
{
  assert(data.face_cells_cached);
  return data.face_cells.getRow<MEM>(f);
}


//
// Edge --> Cell
//
KOKKOS_INLINE_FUNCTION
size_type
MeshCache::getEdgeNumCells(const Entity_ID e) const
{
  assert(data.edge_cells_cached);
  return data.edge_cells.size<MEM>(e);
}


KOKKOS_INLINE_FUNCTION
const Entity_ID&
MeshCache::getEdgeCell(const Entity_ID e, const size_type i) const
{
  assert(data.edge_cells_cached);
  return data.edge_cells.get<MEM>(e,i);
}


template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION
MeshCache::cEntity_ID_View
MeshCache::getEdgeCells(const Entity_ID e) const
{
  return Impl::RaggedGetter<MEM, AP>::get(data.edge_cells_cached, data.edge_cells, nullptr, e);
}


//
// Edge --> Face
//
KOKKOS_INLINE_FUNCTION
size_type
MeshCache::getEdgeNumFaces(const Entity_ID e) const
{
  assert(data.edge_faces_cached);
  return data.edge_faces.size<MEM>(e);
}


KOKKOS_INLINE_FUNCTION
const Entity_ID&
MeshCache::getEdgeFace(const Entity_ID e, const size_type i) const
{
  assert(data.edge_faces_cached);
  return data.edge_faces.get<MEM>(e,i);
}


KOKKOS_INLINE_FUNCTION
MeshCache::cEntity_ID_View
MeshCache::getEdgeFaces(const Entity_ID e) const
{
  assert(data.edge_faces_cached);
  return data.edge_faces.getRow<MEM>(e);
}


//
// Node --> Cell
//
template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION
MeshCache::cEntity_ID_View
MeshCache::getNodeCells(const Entity_ID n) const
{
  return Impl::RaggedGetter<MEM, AP>::get(data.node_cells_cached, data.node_cells, nullptr, n);
}


//
// Node --> Face
//
KOKKOS_INLINE_FUNCTION
MeshCache::cEntity_ID_View
MeshCache::getNodeFaces(const Entity_ID n) const
{
  assert(data.node_faces_cached);
  return data.node_faces.getRow<MEM>(n);
}


//
// Node --> Edge
//
KOKKOS_INLINE_FUNCTION
MeshCache::cEntity_ID_View
MeshCache::getNodeEdges(const Entity_ID n) const
{
  assert(data.node_edges_cached);
  return data.node_edges.getRow<MEM>(n);
}


} // namespace AmanziMesh
} // namespace Amanzi
