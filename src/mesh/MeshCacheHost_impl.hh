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

#include "Teuchos_CommHelpers.hpp"

#include "Geometry.hh"
#include "Iterators.hh"
#include "MeshUtils.hh"
#include "MeshFramework.hh"
#include "MeshCacheHost_decl.hh"
#include "MeshAlgorithms.hh"
#include "MeshSets.hh"

namespace Amanzi {
namespace AmanziMesh {

class SingleFaceMesh;


template <AccessPattern_kind AP>
AmanziGeometry::Point
MeshCacheHost::getCentroid(const Entity_kind kind, const Entity_ID ent) const
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
AmanziGeometry::Point
MeshCacheHost::getCentroid(const Entity_ID ent) const
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
double
MeshCacheHost::getExtent(const Entity_kind kind, const Entity_ID ent) const
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
double
MeshCacheHost::getExtent(const Entity_ID ent) const
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


//===================
//    getFace*
//===================

template <AccessPattern_kind AP>
size_type
MeshCacheHost::getFaceNumCells(const Entity_ID f, const Parallel_kind ptype) const
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
template <AccessPattern_kind AP>
typename MeshCacheHost::cEntity_ID_View
MeshCacheHost::getFaceCells(const Entity_ID f) const
{
  cEntity_ID_View fcells;
  getFaceCells<AP>(f, fcells);
  return fcells;
}

template <AccessPattern_kind AP>
void
MeshCacheHost::getFaceCells(const Entity_ID f, cEntity_ID_View& fcells) const
{
  fcells = Impl::RaggedGetter<MEM, AP>::get(
    data_.face_cells_cached,
    data_.face_cells,
    framework_mesh_,
    [&](const int i) {
      MeshFramework::cEntity_ID_View cells;
      framework_mesh_->getFaceCells(f, cells);
      return cells;
    },
    nullptr,
    f);
}


// Normal vector of a face
template <AccessPattern_kind AP>
AmanziGeometry::Point
MeshCacheHost::getFaceNormal(const Entity_ID f) const
{
  return getFaceNormal<AP>(f, -1, nullptr);
}

template <AccessPattern_kind AP>
AmanziGeometry::Point
MeshCacheHost::getFaceNormal(const Entity_ID f, const Entity_ID c, int* orientation) const
{
  if (this->isSFM()) {
    AmanziGeometry::Point normal = data_.face_normals.get<MemSpace_kind::HOST>(f, 0);
    return normal;
  }

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

        // average normals oriented from lower to higher GIDs
        int pos_i =
          getEntityGID(Entity_kind::CELL, fcells[0]) > getEntityGID(Entity_kind::CELL, fcells[1]) ?
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

  return algorithms_->computeFaceNormal(*this, f, c, orientation);
}

template <AccessPattern_kind AP>
AmanziGeometry::Point
MeshCacheHost::getFaceCentroid(const Entity_ID f) const
{
  auto cf = [&](const int i) { return algorithms_->computeFaceCentroid(*this, i); };
  return Impl::Getter<MEM, AP>::get(
    data_.face_geometry_cached, data_.face_centroids, framework_mesh_, nullptr, cf, f);
}

template <AccessPattern_kind AP>
double
MeshCacheHost::getFaceArea(const Entity_ID f) const
{
  auto cf = [&](const Entity_ID i) { return algorithms_->computeFaceArea(*this, i); };

  return Impl::Getter<MEM, AP>::get(
    data_.face_geometry_cached, data_.face_areas, framework_mesh_, nullptr, cf, f);
}

template <AccessPattern_kind AP>
typename MeshCacheHost::cPoint_View
MeshCacheHost::getFaceCoordinates(const Entity_ID f) const
{
  return Impl::RaggedGetter<MEM, AP>::get(
    data_.face_coordinates_cached,
    data_.face_coordinates,
    framework_mesh_,
    [&](const Entity_ID i) { return framework_mesh_->getFaceCoordinates(i); },
    nullptr,
    f);
}

//===================
//    getCell*
//===================

// extent
template <AccessPattern_kind AP>
double
MeshCacheHost::getCellVolume(const Entity_ID c) const
{
  auto cf = [&](const int i) { return algorithms_->computeCellVolume(*this, i); };

  return Impl::Getter<MEM, AP>::get(
    data_.cell_geometry_cached, data_.cell_volumes, framework_mesh_, nullptr, cf, c);
}

template <AccessPattern_kind AP>
size_type
MeshCacheHost::getCellNumFaces(const Entity_ID c) const
{
  static_assert(AP != AccessPattern_kind::COMPUTE);
  static_assert(AP != AccessPattern_kind::FRAMEWORK);
  // this is where a generic function would probably help?
  if constexpr (AP == AccessPattern_kind::CACHE) {
    assert(data_.cell_faces_cached);
    return data_.cell_faces.size<MEM>(c);
  } else {
    if (data_.cell_faces_cached) return getCellNumFaces<AccessPattern_kind::CACHE>(c);
    return getCellFaces(c).size();
  }
}


template <AccessPattern_kind AP>
typename MeshCacheHost::cEntity_ID_View
MeshCacheHost::getCellFaces(const Entity_ID c) const
{
  cEntity_ID_View cfaces;
  getCellFaces<AP>(c, cfaces);
  return cfaces;
}

template <AccessPattern_kind AP>
void
MeshCacheHost::getCellFaces(const Entity_ID c, cEntity_ID_View& cfaces) const
{
  cfaces = Impl::RaggedGetter<MEM, AP>::get(
    data_.cell_faces_cached,
    data_.cell_faces,
    framework_mesh_,
    [&](const int i) {
      MeshFramework::cEntity_ID_View cf;
      framework_mesh_->getCellFaces(i, cf);
      return cf;
    },
    nullptr,
    c);
}

template <AccessPattern_kind AP>
AmanziGeometry::Point
MeshCacheHost::getCellCentroid(const Entity_ID c) const
{
  auto cf = [&](const int i) { return algorithms_->computeCellCentroid(*this, i); };

  return Impl::Getter<MEM, AP>::get(
    data_.cell_geometry_cached, data_.cell_centroids, framework_mesh_, nullptr, cf, c);
}

template <AccessPattern_kind AP>
size_type
MeshCacheHost::getCellNumNodes(const Entity_ID c) const
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
typename MeshCacheHost::cPoint_View
MeshCacheHost::getCellCoordinates(const Entity_ID c) const
{
  auto cf = [&](const int i) { return Impl::computeCellCoordinates(*this, c); };

  return Impl::RaggedGetter<MEM, AP>::get(
    data_.cell_coordinates_cached,
    data_.cell_coordinates,
    framework_mesh_,
    [&](const int i) { return framework_mesh_->getCellCoordinates(i); },
    cf,
    c);
}

template <AccessPattern_kind AP>
size_type
MeshCacheHost::getCellNumEdges(const Entity_ID c) const
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


template <AccessPattern_kind AP>
typename MeshCacheHost::cEntity_ID_View
MeshCacheHost::getCellEdges(const Entity_ID c) const
{
  return Impl::RaggedGetter<MEM, AP>::get(
    data_.cell_edges_cached,
    data_.cell_edges,
    framework_mesh_,
    [&](const int i) {
      MeshFramework::cEntity_ID_View ce;
      framework_mesh_->getCellEdges(i, ce);
      return ce;
    },
    nullptr,
    c);
}


template <AccessPattern_kind AP>
void
MeshCacheHost::getCellEdges(const Entity_ID c, cEntity_ID_View& cedges) const
{
  cedges = Impl::RaggedGetter<MEM, AP>::get(
    data_.cell_edges_cached,
    data_.cell_edges,
    framework_mesh_,
    [&](const int i) {
      MeshFramework::cEntity_ID_View ce;
      framework_mesh_->getCellEdges(i, ce);
      return ce;
    },
    nullptr,
    c);
}

//===================
//    getNode*
//===================


template <AccessPattern_kind AP>
AmanziGeometry::Point
MeshCacheHost::getNodeCoordinate(const Entity_ID n) const
{
  return Impl::Getter<MEM, AP>::get(
    data_.node_coordinates_cached,
    data_.node_coordinates,
    framework_mesh_,
    [&](const int i) { return framework_mesh_->getNodeCoordinate(i); },
    nullptr,
    n);
}

template <AccessPattern_kind AP>
typename MeshCacheHost::cEntity_ID_View
MeshCacheHost::getNodeCells(const Entity_ID n, const Parallel_kind ptype) const
{
  cEntity_ID_View cells;
  getNodeCells<AP>(n, ptype, cells);
  return cells;
}

template <AccessPattern_kind AP>
void
MeshCacheHost::getNodeCells(const Entity_ID n,
                            const Parallel_kind ptype,
                            cEntity_ID_View& cells) const
{
  cells = Impl::RaggedGetter<MEM, AP>::get(
    data_.node_cells_cached,
    data_.node_cells,
    framework_mesh_,
    [&](const int i) {
      MeshFramework::cEntity_ID_View lcells;
      framework_mesh_->getNodeCells(i, lcells);
      return lcells;
    },
    nullptr,
    n);
}

template <AccessPattern_kind AP>
typename MeshCacheHost::cEntity_ID_View
MeshCacheHost::getNodeFaces(const Entity_ID n) const
{
  cEntity_ID_View faces;
  getNodeFaces<AP>(n, faces);
  return faces;
}

template <AccessPattern_kind AP>
void
MeshCacheHost::getNodeFaces(const Entity_ID n, cEntity_ID_View& faces) const
{
  faces = Impl::RaggedGetter<MEM, AP>::get(
    data_.node_faces_cached,
    data_.node_faces,
    framework_mesh_,
    [&](const int i) {
      MeshFramework::cEntity_ID_View lfaces;
      framework_mesh_->getNodeFaces(i, lfaces);
      return lfaces;
    },
    nullptr,
    n);
}

//==================
//    getEdge*
//==================

template <AccessPattern_kind AP>
typename MeshCacheHost::cEntity_ID_View
MeshCacheHost::getEdgeNodes(const Entity_ID e) const
{
  return Impl::RaggedGetter<MEM, AP>::get(
    data_.edge_nodes_cached,
    data_.edge_nodes,
    framework_mesh_,
    [&](const int i) {
      MeshFramework::cEntity_ID_View nodes;
      framework_mesh_->getEdgeNodes(i, nodes);
      return nodes;
    },
    nullptr,
    e);
}

template <AccessPattern_kind AP>
Entity_ID
MeshCacheHost::getEdgeNode(const Entity_ID e, const size_type i) const
{
  // Compute list and use only one?
  cEntity_ID_View nodes;
  getEdgeNodes(e, nodes);
  return nodes[i];
}

template <AccessPattern_kind AP>
void
MeshCacheHost::getEdgeNodes(const Entity_ID e, cEntity_ID_View& nodes) const
{
  nodes = Impl::RaggedGetter<MEM, AP>::get(
    data_.edge_nodes_cached,
    data_.edge_nodes,
    framework_mesh_,
    [&](const int i) {
      MeshFramework::cEntity_ID_View lnodes;
      framework_mesh_->getEdgeNodes(i, lnodes);
      return lnodes;
    },
    nullptr,
    e);
}

template <AccessPattern_kind AP>
typename MeshCacheHost::cPoint_View
MeshCacheHost::getEdgeCoordinates(const Entity_ID c) const
{
  return Impl::RaggedGetter<MEM, AP>::get(
    data_.edge_coordinates_cached,
    data_.edge_coordinates,
    framework_mesh_,
    [&](const int i) { return framework_mesh_->getEdgeCoordinates(i); },
    nullptr,
    c);
}

template <AccessPattern_kind AP>
AmanziGeometry::Point
MeshCacheHost::getEdgeCentroid(const Entity_ID c) const
{
  auto cf = [&](const int i) { return algorithms_->computeEdgeCentroid(*this, i); };

  return Impl::Getter<MEM, AP>::get(
    data_.edge_geometry_cached, data_.edge_centroids, framework_mesh_, nullptr, cf, c);
}

template <AccessPattern_kind AP>
AmanziGeometry::Point
MeshCacheHost::getEdgeVector(const Entity_ID e) const
{
  auto cf = [&](const int i) { return algorithms_->computeEdgeVector(*this, i, -1, nullptr); };

  return Impl::Getter<MEM, AP>::get(
    data_.edge_geometry_cached, data_.edge_vectors, framework_mesh_, nullptr, cf, e);
}


template <AccessPattern_kind AP>
double
MeshCacheHost::getEdgeLength(const Entity_ID e) const
{
  auto cf = [&](const int i) { return algorithms_->computeEdgeLength(*this, i); };

  return Impl::Getter<MEM, AP>::get(
    data_.edge_lengths_cached, data_.edge_lengths, framework_mesh_, nullptr, cf, e);
}


template <AccessPattern_kind AP>
typename MeshCacheHost::cEntity_ID_View
MeshCacheHost::getEdgeCells(const Entity_ID e) const
{
  cEntity_ID_View cells;
  getEdgeCells<AP>(e, cells);
  return cells;
}


template <AccessPattern_kind AP>
void
MeshCacheHost::getEdgeCells(const Entity_ID e, cEntity_ID_View& cells) const
{
  cells = Impl::RaggedGetter<MEM, AP>::get(
    data_.edge_cells_cached,
    data_.edge_cells,
    framework_mesh_,
    [&](const int i) {
      MeshFramework::cEntity_ID_View lcells;
      framework_mesh_->getEdgeCells(i, lcells);
      return lcells;
    },
    nullptr,
    e);
}

template <AccessPattern_kind AP>
typename MeshCacheHost::cEntity_ID_View
MeshCacheHost::getEdgeFaces(const Entity_ID e) const
{
  cEntity_ID_View faces;
  getEdgeFaces<AP>(e, faces);
  return faces;
}

template <AccessPattern_kind AP>
void
MeshCacheHost::getEdgeFaces(const Entity_ID e, cEntity_ID_View& faces) const
{
  faces = Impl::RaggedGetter<MEM, AP>::get(
    data_.edge_faces_cached,
    data_.edge_faces,
    framework_mesh_,
    [&](const int i) {
      MeshFramework::cEntity_ID_View lfaces;
      framework_mesh_->getEdgeFaces(i, lfaces);
      return lfaces;
    },
    nullptr,
    e);
}


template <AccessPattern_kind AP>
typename MeshCacheHost::cEntity_ID_View
MeshCacheHost::getNodeEdges(const Entity_ID n) const
{
  cEntity_ID_View edges;
  getNodeEdges<AP>(n, edges);
  return edges;
}


template <AccessPattern_kind AP>
void
MeshCacheHost::getNodeEdges(const Entity_ID n, cEntity_ID_View& edges) const
{
  edges = Impl::RaggedGetter<MEM, AP>::get(
    data_.node_edges_cached,
    data_.node_edges,
    framework_mesh_,
    [&](const int i) {
      MeshFramework::cEntity_ID_View ledges;
      framework_mesh_->getNodeEdges(i, ledges);
      return ledges;
    },
    nullptr,
    n);
}


} // namespace AmanziMesh
} // namespace Amanzi
