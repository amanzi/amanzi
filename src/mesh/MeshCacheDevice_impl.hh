#pragma once

#include "Teuchos_CommHelpers.hpp"

#include "Geometry.hh"
#include "Iterators.hh"
#include "MeshUtils.hh"
#include "MeshFramework.hh"
#include "MeshCacheDevice_decl.hh"
#include "MeshAlgorithms.hh"
#include "MeshSets.hh"

namespace Amanzi::AmanziMesh {

template <AccessPattern_kind AP>
AmanziGeometry::Point
MeshCacheDevice::getCentroid(const Entity_kind kind, const Entity_ID ent) const
{
  switch (kind) {
  case (Entity_kind::CELL):
    return getCellCentroid<AP>(ent);
  case (Entity_kind::FACE):
    return getFaceCentroid<AP>(ent);
  case (Entity_kind::BOUNDARY_FACE):
    return getFaceCentroid<AP>(getBoundaryFaces()(ent));
  case (Entity_kind::EDGE):
    return getEdgeCentroid<AP>(ent);
  case (Entity_kind::NODE):
    return getNodeCoordinate<AP>(ent);
  default:
    assert(false);
  }
  return AmanziGeometry::Point();
}

template <Entity_kind EK, AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION AmanziGeometry::Point
MeshCacheDevice::getCentroid(const Entity_ID ent) const
{
  static_assert(EK != Entity_kind::UNKNOWN && EK != Entity_kind::BOUNDARY_NODE, "No centroid information for Entity_kind::UNKNOWN or Entity_kind::BOUNDARY_NODE"); 
  if constexpr (EK == Entity_kind::CELL) {
    return getCellCentroid<AP>(ent);
  } else if constexpr (EK == Entity_kind::FACE) {
    return getFaceCentroid<AP>(ent);
  } else if constexpr (EK == Entity_kind::BOUNDARY_FACE) {
    return getFaceCentroid<AP>(getBoundaryFaces()(ent));
  } else if constexpr (EK == Entity_kind::EDGE) {
    return getEdgeCentroid<AP>(ent);
  } else if constexpr (EK == Entity_kind::NODE) {
    return getNodeCoordinate<AP>(ent);
  }  
}


template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION double
MeshCacheDevice::getExtent(const Entity_kind kind, const Entity_ID ent) const
{
  switch (kind) {
  case (Entity_kind::CELL):
    return getCellVolume<AP>(ent);
  case (Entity_kind::FACE):
    return getFaceArea<AP>(ent);
  case (Entity_kind::EDGE):
    return getEdgeLength<AP>(ent);
  default:
    AMANZI_ASSERT(false);
  }
  return -1.0;
}


template <Entity_kind EK, AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION double
MeshCacheDevice::getExtent(const Entity_ID ent) const
{
  static_assert(EK == Entity_kind::CELL || EK == Entity_kind::FACE || EK == Entity_kind::EDGE, "getExtent only supports CELL, FACE and EDGE"); 
  if constexpr (EK == Entity_kind::CELL) {
    return getCellVolume<AP>(ent);
  } else if constexpr (EK == Entity_kind::FACE) {
    return getFaceArea<AP>(ent);
  } else if constexpr (EK == Entity_kind::EDGE) {
    return getEdgeLength<AP>(ent);
  }
}

template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION size_type
MeshCacheDevice::getFaceNumCells(const Entity_ID f, const Parallel_kind ptype) const
{
  static_assert(AP != AccessPattern_kind::COMPUTE);
  static_assert(AP != AccessPattern_kind::FRAMEWORK);
  // this is where a generic function would probably help?
  if constexpr (AP == AccessPattern_kind::CACHE) {
    assert(data_.face_cells_cached);
    if (ptype == Parallel_kind::ALL) {
      return data_.face_cells.size<MEM>(f);
    } else {
      int count = 0;
      int n_all = data_.face_cells.size<MEM>(f);
      for (int j = 0; j != n_all; ++j) {
        if (getFaceCell(f, j) < ncells_owned)
          ++count;
        else
          break;
      }
      return count;
    }
  } else {
    if (data_.face_cells_cached) return getFaceNumCells<AccessPattern_kind::CACHE>(f, ptype);
    return getFaceCells(f).size();
  }
}

// Normal vector of a face
template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION AmanziGeometry::Point
MeshCacheDevice::getFaceNormal(const Entity_ID f) const
{
  return getFaceNormal<AP>(f, -1, nullptr);
}


template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION AmanziGeometry::Point
MeshCacheDevice::getFaceNormal(const Entity_ID f, const Entity_ID c, int* orientation) const
{
  if (this->isSFM()) {
    AmanziGeometry::Point normal = data_.face_normals.get<MemSpace_kind::HOST>(f, 0);
    return normal;
  }

  if constexpr (MEM == MemSpace_kind::DEVICE) { assert(data_.face_geometry_cached); }

  if (data_.face_geometry_cached) {
    if (c < 0) {
      assert(orientation == nullptr);
      auto normal = data_.face_normals.get<MEM>(f, 0);
      auto nat_dir = data_.face_normal_orientations.get<MEM>(f, 0);
      if (abs(nat_dir) == 1) {
        return nat_dir * normal;
      } else {
        auto fcells = getFaceCells(f);
        assert(fcells.size() == 2);
        assert(data_.cell_global_indices_cached);

        // average normals oriented from lower to higher GIDs
        int pos_i = view<MemSpace_kind::DEVICE>(data_.cell_global_indices)(fcells[0]) >
                        view<MemSpace_kind::HOST>(data_.cell_global_indices)(fcells[1]) ?
                      0 :
                      1;
        return (data_.face_normals.get<MEM>(f, 1 - pos_i) - data_.face_normals.get<MEM>(f, pos_i)) /
               2;
      }
    }

    auto fcells = getFaceCells(f);
    std::size_t i;
    for (i = 0; i != fcells.size(); ++i)
      if (fcells[i] == c) break;

    if (orientation) *orientation = data_.face_normal_orientations.get<MEM>(f, i) > 0 ? 1 : -1;
    return data_.face_normals.get<MEM>(f, i);
  }

  assert(false && "No access to cache/framework/compute available");
  return AmanziGeometry::Point();
}

template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION AmanziGeometry::Point
MeshCacheDevice::getFaceCentroid(const Entity_ID f) const
{
  return Impl::Getter<MEM, AP>::get(data_.face_geometry_cached, data_.face_centroids, nullptr, f);
}


template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION double
MeshCacheDevice::getFaceArea(const Entity_ID f) const
{
  return Impl::Getter<MEM, AP>::get(data_.face_geometry_cached, data_.face_areas, nullptr, f);
}

template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION typename MeshCacheDevice::cPoint_View
MeshCacheDevice::getFaceCoordinates(const Entity_ID f) const
{
  return Impl::RaggedGetter<MEM, AP>::get(
    data_.face_coordinates_cached, data_.face_coordinates, nullptr, f);
}


template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION double
MeshCacheDevice::getCellVolume(const Entity_ID c) const
{
  return Impl::Getter<MEM, AP>::get(data_.cell_geometry_cached, data_.cell_volumes, nullptr, c);
}

template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION size_type
MeshCacheDevice::getCellNumNodes(const Entity_ID c) const
{
  static_assert(AP != AccessPattern_kind::COMPUTE);
  static_assert(AP != AccessPattern_kind::FRAMEWORK);
  // this is where a generic function would probably help?
  if constexpr (AP == AccessPattern_kind::CACHE) {
    assert(data_.cell_nodes_cached);
    return data_.cell_nodes.size<MEM>(c);
  } else {
    if (data_.cell_nodes_cached) return getCellNumNodes<AccessPattern_kind::CACHE>(c);
    return getCellNodes(c).size();
  }
}

template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION AmanziGeometry::Point
MeshCacheDevice::getCellCentroid(const Entity_ID c) const
{
  return Impl::Getter<MEM, AP>::get(data_.cell_geometry_cached, data_.cell_centroids, nullptr, c);
}

template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION typename MeshCacheDevice::cPoint_View
MeshCacheDevice::getCellCoordinates(const Entity_ID c) const
{
  return Impl::RaggedGetter<MEM, AP>::get(
    data_.cell_coordinates_cached, data_.cell_coordinates, nullptr, c);
}


template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION size_type
MeshCacheDevice::getCellNumEdges(const Entity_ID c) const
{
  static_assert(AP != AccessPattern_kind::COMPUTE);
  static_assert(AP != AccessPattern_kind::FRAMEWORK);
  // this is where a generic function would probably help?
  if constexpr (AP == AccessPattern_kind::CACHE) {
    assert(data_.cell_edges_cached);
    return data_.cell_edges.size<MEM>(c);
  } else {
    if (data_.cell_edges_cached) return getCellNumEdges<AccessPattern_kind::CACHE>(c);
    return getCellEdges(c).size();
  }
}

//===================
//    getNode*
//===================

template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION AmanziGeometry::Point
MeshCacheDevice::getNodeCoordinate(const Entity_ID n) const
{
  return Impl::Getter<MEM, AP>::get(
    data_.node_coordinates_cached, data_.node_coordinates, nullptr, n);
}

template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION typename MeshCacheDevice::cEntity_ID_View
MeshCacheDevice::getNodeCells(const Entity_ID n, const Parallel_kind ptype) const
{
  cEntity_ID_View cells;
  getNodeCells<AP>(n, ptype, cells);
  return cells;
}

template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION void
MeshCacheDevice::getNodeCells(const Entity_ID n,
                              const Parallel_kind ptype,
                              cEntity_ID_View& cells) const
{
  cells = Impl::RaggedGetter<MEM, AP>::get(data_.node_cells_cached, data_.node_cells, nullptr, n);
}


template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION typename MeshCacheDevice::cEntity_ID_View
MeshCacheDevice::getNodeFaces(const Entity_ID n) const
{
  cEntity_ID_View faces;
  getNodeFaces<AP>(n, faces);
  return faces;
}

template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION void
MeshCacheDevice::getNodeFaces(const Entity_ID n, cEntity_ID_View& faces) const
{
  faces = Impl::RaggedGetter<MEM, AP>::get(data_.node_faces_cached, data_.node_faces, nullptr, n);
}


//==================
//    getEdge*
//==================

template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION typename MeshCacheDevice::cEntity_ID_View
MeshCacheDevice::getEdgeNodes(const Entity_ID e) const
{
  return Impl::RaggedGetter<MEM, AP>::get(data_.edge_nodes_cached, data_.edge_nodes, nullptr, e);
}


template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION Entity_ID
MeshCacheDevice::getEdgeNode(const Entity_ID e, const size_type i) const
{
  // Compute list and use only one?
  cEntity_ID_View nodes;
  getEdgeNodes(e, nodes);
  return nodes[i];
}


template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION void
MeshCacheDevice::getEdgeNodes(const Entity_ID e, cEntity_ID_View& nodes) const
{
  nodes = Impl::RaggedGetter<MEM, AP>::get(data_.edge_nodes_cached, data_.edge_nodes, nullptr, e);
}


template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION typename MeshCacheDevice::cPoint_View
MeshCacheDevice::getEdgeCoordinates(const Entity_ID c) const
{
  return Impl::RaggedGetter<MEM, AP>::get(
    data_.edge_coordinates_cached, data_.edge_coordinates, nullptr, c);
}


template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION AmanziGeometry::Point
MeshCacheDevice::getEdgeCentroid(const Entity_ID c) const
{
  return Impl::Getter<MEM, AP>::get(data_.edge_geometry_cached, data_.edge_centroids, nullptr, c);
}


template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION AmanziGeometry::Point
MeshCacheDevice::getEdgeVector(const Entity_ID e) const
{
  return Impl::Getter<MEM, AP>::get(data_.edge_geometry_cached, data_.edge_vectors, nullptr, e);
}


template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION double
MeshCacheDevice::getEdgeLength(const Entity_ID e) const
{
  return Impl::Getter<MEM, AP>::get(data_.edge_lengths_cached, data_.edge_lengths, nullptr, e);
}


template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION typename MeshCacheDevice::cEntity_ID_View
MeshCacheDevice::getEdgeCells(const Entity_ID e) const
{
  cEntity_ID_View cells;
  getEdgeCells<AP>(e, cells);
  return cells;
}


template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION void
MeshCacheDevice::getEdgeCells(const Entity_ID e, cEntity_ID_View& cells) const
{
  cells = Impl::RaggedGetter<MEM, AP>::get(data_.edge_cells_cached, data_.edge_cells, nullptr, e);
}


template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION typename MeshCacheDevice::cEntity_ID_View
MeshCacheDevice::getEdgeFaces(const Entity_ID e) const
{
  cEntity_ID_View faces;
  getEdgeFaces<AP>(e, faces);
  return faces;
}


template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION void
MeshCacheDevice::getEdgeFaces(const Entity_ID e, cEntity_ID_View& faces) const
{
  faces = Impl::RaggedGetter<MEM, AP>::get(data_.edge_faces_cached, data_.edge_faces, nullptr, e);
}


template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION typename MeshCacheDevice::cEntity_ID_View
MeshCacheDevice::getNodeEdges(const Entity_ID n) const
{
  cEntity_ID_View edges;
  getNodeEdges<AP>(n, edges);
  return edges;
}


template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION void
MeshCacheDevice::getNodeEdges(const Entity_ID n, cEntity_ID_View& edges) const
{
  edges = Impl::RaggedGetter<MEM, AP>::get(data_.node_edges_cached, data_.node_edges, nullptr, n);
}

} // namespace Amanzi::AmanziMesh
