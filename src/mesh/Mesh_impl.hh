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

#include "Mesh_decl.hh"

#include "GeometricModel.hh"
#include "MeshFramework.hh"
#include "MeshAlgorithms.hh"
#include "MeshSets.hh"
#include "MeshColumns.hh"


namespace Amanzi {
namespace AmanziMesh {

// -----------------
// Set Entitites
// -----------------

//
// TODO --etc
// This should be updated -- only cache the ALL set, then construct the
// OWNED or GHOST set on demand.  No need to save both OWNED and ALL.
//
// TODO --etc
// Don't cache "ALL" or "BOUNDARY" labeled sets.  Are there others that are
// just as fast to create as to cache?
//
template<MemSpace_kind MEM>
auto // cEntity_ID_View on MEM
Mesh::getSetEntities(const std::string& region_name,
               const Entity_kind kind,
               const Parallel_kind ptype) const
{
  auto key = std::make_tuple(region_name, kind, ptype);
  auto key_all = std::make_tuple(region_name, kind, Parallel_kind::ALL);

  if (!sets_->count(key_all)) {
    auto region = getGeometricModel()->FindRegion(region_name);
    if (region == Teuchos::null) {
      Errors::Message msg;
      msg << "Cannot find region of name \"" << region_name << "\" in the geometric model.";
      Exceptions::amanzi_throw(msg);
    }

    (*sets_)[key_all] = asDualView(resolveMeshSet(*region, kind, Parallel_kind::ALL, *this));
  }

  if (!sets_->count(key)) {
    if (ptype == Parallel_kind::OWNED) {
      auto v_all = view<MemSpace_kind::HOST>(sets_->at(key_all));
      size_t n_ents = getNumEntities(kind, Parallel_kind::OWNED);
      size_t i;
      for (i = 0; i != v_all.size(); ++i)
        if (v_all(i) >= n_ents) break;
      auto v_owned = Kokkos::subview(v_all, Kokkos::make_pair((size_t)0, i));
      (*sets_)[key] = asDualView(v_owned);

    } else if (ptype == Parallel_kind::GHOST) {
      auto v_all = view<MemSpace_kind::HOST>(sets_->at(key_all));
      size_t n_ents = getNumEntities(kind, Parallel_kind::OWNED);
      size_t i;
      for (i = 0; i != v_all.size(); ++i)
        if (v_all(i) >= n_ents) break;
      auto v_ghosted = Kokkos::subview(v_all, Kokkos::make_pair(i, v_all.size()));
      (*sets_)[key] = asDualView(v_ghosted);
    } else {
      AMANZI_ASSERT(false);
    }
  }

  return view<MEM>(sets_->at(key));
}


template<MemSpace_kind MEM>
auto // Kokkos::pair<cEntity_ID_View, cDouble_View> on MEM
Mesh::getSetEntitiesAndVolumeFractions(const std::string& region_name,
        const Entity_kind kind,
        const Parallel_kind ptype) const
{
  auto key = std::make_tuple(region_name, kind, ptype);

  if (!set_vol_fracs_->count(key)) {
    auto region = getGeometricModel()->FindRegion(region_name);
    if (region == Teuchos::null) {
      Errors::Message msg;
      msg << "Cannot find region of name \"" << region_name << "\" in the geometric model.";
      Exceptions::amanzi_throw(msg);
    }

    Double_View vol_fracs_list;
    (*sets_)[key] =
      asDualView(resolveMeshSetVolumeFractions(*region, kind, ptype, vol_fracs_list, *this));
    (*set_vol_fracs_)[key] = asDualView<double>(vol_fracs_list);
  }

  return Kokkos::pair(view<MEM>(sets_->at(key)), view<MEM>(set_vol_fracs_->at(key)));
}


// ----------------------
// Entity meta-data
// ----------------------
template<MemSpace_kind MEM>
auto
Mesh::getEntityParents(const Entity_kind kind) const
{
  AMANZI_ASSERT(data_.parent_entities_cached);
  switch (kind) {
  case Entity_kind::CELL:
    return view<MEM>(data_.parent_cells);
    break;
  case Entity_kind::FACE:
    return view<MEM>(data_.parent_faces);
    break;
  case Entity_kind::EDGE:
    return view<MEM>(data_.parent_edges);
    break;
  case Entity_kind::NODE:
    return view<MEM>(data_.parent_nodes);
    break;
  default: {
  }
  }
  AMANZI_ASSERT(false);
  return view<MEM>(data_.parent_cells);
}



// ----------------------
// Entity relations
// ----------------------
// a list of ALL face LIDs that are on the boundary
inline typename Mesh::cEntity_ID_View
Mesh::getBoundaryFaces() const
{
  return maps_->getBoundaryFaces<MEM>();
}

// a list of ALL node LIDs that are on the boundary
inline typename Mesh::cEntity_ID_View
Mesh::getBoundaryNodes() const
{
  return maps_->getBoundaryNodes<MEM>();
}

inline Entity_ID
Mesh::getBoundaryFaceFace(const Entity_ID bf) const
{
  return getBoundaryFaces()(bf);
}

inline Entity_ID
Mesh::getBoundaryNodeNode(const Entity_ID bf) const
{
  return getBoundaryNodes()(bf);
}



// ----------------------
// Geometry: Coordinates
// ----------------------
template <AccessPattern_kind AP>
typename Mesh::cPoint_View
Mesh::getEdgeCoordinates(const Entity_ID c) const
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
typename Mesh::cPoint_View
Mesh::getFaceCoordinates(const Entity_ID f) const
{
  return Impl::RaggedGetter<MEM, AP>::get(
    data_.face_coordinates_cached,
    data_.face_coordinates,
    framework_mesh_,
    [&](const Entity_ID i) { return framework_mesh_->getFaceCoordinates(i); },
    nullptr,
    f);
}


template <AccessPattern_kind AP>
typename Mesh::cPoint_View
Mesh::getCellCoordinates(const Entity_ID c) const
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


// ----------------------
// Geometry: Centroids
// ----------------------
template <AccessPattern_kind AP>
AmanziGeometry::Point
Mesh::getCellCentroid(const Entity_ID c) const
{
  auto cf = [&](const int i) { return algorithms_->computeCellCentroid(*this, i); };

  return Impl::Getter<MEM, AP>::get(
    data_.cell_geometry_cached, data_.cell_centroids, framework_mesh_, nullptr, cf, c);
}


template <AccessPattern_kind AP>
AmanziGeometry::Point
Mesh::getFaceCentroid(const Entity_ID f) const
{
  auto cf = [&](const int i) { return algorithms_->computeFaceCentroid(*this, i); };
  return Impl::Getter<MEM, AP>::get(
    data_.face_geometry_cached, data_.face_centroids, framework_mesh_, nullptr, cf, f);
}


template <AccessPattern_kind AP>
AmanziGeometry::Point
Mesh::getEdgeCentroid(const Entity_ID c) const
{
  auto cf = [&](const int i) { return algorithms_->computeEdgeCentroid(*this, i); };

  return Impl::Getter<MEM, AP>::get(
    data_.edge_geometry_cached, data_.edge_centroids, framework_mesh_, nullptr, cf, c);
}


template <AccessPattern_kind AP>
AmanziGeometry::Point
Mesh::getCentroid(const Entity_kind kind, const Entity_ID ent) const
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
AmanziGeometry::Point
Mesh::getCentroid(const Entity_ID ent) const
{
  static_assert(EK != Entity_kind::UNKNOWN &&
                EK != Entity_kind::BOUNDARY_NODE,
                "No centroid information for Entity_kind::UNKNOWN or Entity_kind::BOUNDARY_NODE");

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
  return AmanziGeometry::Point();
}


// ----------------------
// Geometry: Extents
// ----------------------
template <AccessPattern_kind AP>
double
Mesh::getCellVolume(const Entity_ID c) const
{
  auto cf = [&](const int i) { return algorithms_->computeCellVolume(*this, i); };
  return Impl::Getter<MEM, AP>::get(
    data_.cell_geometry_cached, data_.cell_volumes, framework_mesh_, nullptr, cf, c);
}


template <AccessPattern_kind AP>
double
Mesh::getFaceArea(const Entity_ID f) const
{
  auto cf = [&](const Entity_ID i) { return algorithms_->computeFaceArea(*this, i); };
  return Impl::Getter<MEM, AP>::get(
    data_.face_geometry_cached, data_.face_areas, framework_mesh_, nullptr, cf, f);
}


template <AccessPattern_kind AP>
double
Mesh::getEdgeLength(const Entity_ID e) const
{
  auto cf = [&](const int i) { return algorithms_->computeEdgeLength(*this, i); };
  return Impl::Getter<MEM, AP>::get(
    data_.edge_lengths_cached, data_.edge_lengths, framework_mesh_, nullptr, cf, e);
}


template <AccessPattern_kind AP>
double
Mesh::getExtent(const Entity_kind kind, const Entity_ID ent) const
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
    AMANZI_ASSERT(false);
  }
  return -1.0;
}


template <Entity_kind EK, AccessPattern_kind AP>
double
Mesh::getExtent(const Entity_ID ent) const
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
  return -1.0;
}


//-----------------------
// Geometry: other
//-----------------------
template <AccessPattern_kind AP>
AmanziGeometry::Point
Mesh::getFaceNormal(const Entity_ID f) const
{
  return getFaceNormal<AP>(f, -1, nullptr);
}


// Normal vector of a face
template <AccessPattern_kind AP>
AmanziGeometry::Point
Mesh::getFaceNormal(const Entity_ID f, const Entity_ID c, int* orientation) const
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
Mesh::getEdgeVector(const Entity_ID e) const
{
  auto cf = [&](const int i) { return algorithms_->computeEdgeVector(*this, i, -1, nullptr); };

  return Impl::Getter<MEM, AP>::get(
    data_.edge_geometry_cached, data_.edge_vectors, framework_mesh_, nullptr, cf, e);
}


// ------------------------
// Downward adjacencies
// ------------------------
template <AccessPattern_kind AP>
typename Mesh::cEntity_ID_View
Mesh::getCellEdges(const Entity_ID c) const
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
typename Mesh::cEntity_ID_View
Mesh::getCellNodes(const Entity_ID c) const
{
  return Impl::RaggedGetter<MEM, AP>::get(
    data_.cell_nodes_cached,
    data_.cell_nodes,
    framework_mesh_,
    [&](const int i) {
      MeshFramework::cEntity_ID_View ce;
      framework_mesh_->getCellNodes(i, ce);
      return ce;
    },
    nullptr,
    c);
}


template <AccessPattern_kind AP>
void
Mesh::getCellNodes(const Entity_ID c, typename Mesh::cEntity_ID_View& nodes) const
{
  nodes = getCellNodes<AP>(c);
}


// ------------------------
// Upward adjacencies
// ------------------------
template <AccessPattern_kind AP>
typename Mesh::cEntity_ID_View
Mesh::getEdgeCells(const Entity_ID c) const
{
  return Impl::RaggedGetter<MEM, AP>::get(
    data_.edge_cells_cached,
    data_.edge_cells,
    framework_mesh_,
    [&](const int i) {
      MeshFramework::cEntity_ID_View ce;
      framework_mesh_->getEdgeCells(i, ce);
      return ce;
    },
    nullptr,
    c);
}


template <AccessPattern_kind AP>
typename Mesh::cEntity_ID_View
Mesh::getNodeCells(const Entity_ID c) const
{
  return Impl::RaggedGetter<MEM, AP>::get(
    data_.node_cells_cached,
    data_.node_cells,
    framework_mesh_,
    [&](const int i) {
      MeshFramework::cEntity_ID_View ce;
      framework_mesh_->getNodeCells(i, ce);
      return ce;
    },
    nullptr,
    c);
}

} // namespace AmanziMesh
} // namespace Amanzi
