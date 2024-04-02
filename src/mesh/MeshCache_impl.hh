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

#include "Geometry.hh"
#include "Iterators.hh"
#include "MeshUtils.hh"
#include "MeshFramework.hh"
#include "MeshCache_decl.hh"
#include "MeshAlgorithms.hh"
#include "MeshSets.hh"

namespace Amanzi {
namespace AmanziMesh {

class SingleFaceMesh;

template <MemSpace_kind MEM>
template <MemSpace_kind MEM_OTHER>
MeshCache<MEM>::MeshCache(MeshCache<MEM_OTHER>& other)
  : MeshCacheBase(other),
    parent_(other.getParentMesh() == Teuchos::null ? Teuchos::null :
                                                     onMemSpace<MEM>(other.getParentMesh())),
    vis_mesh_(other.getVisMeshPtr() == Teuchos::null ? Teuchos::null :
                                                       onMemSpace<MEM>(other.getVisMeshPtr()))
{}


template <MemSpace_kind MEM>
void
MeshCache<MEM>::setParentMesh(const Teuchos::RCP<const MeshCache>& parent)
{
  if (parent_ != Teuchos::null && parent_ != parent) {
    Errors::Message msg("MeshCache::setParentMesh given conflicting parent mesh.");
    Exceptions::amanzi_throw(msg);
  } else {
    parent_ = parent;
    cacheParentEntities();
  }
}

template <MemSpace_kind MEM>
void
MeshCache<MEM>::setMeshFramework(const Teuchos::RCP<MeshFramework>& framework_mesh)
{
  framework_mesh_ = framework_mesh;
  // always save the algorithms, so we can throw away the data
  has_edges_ = framework_mesh->hasEdges();
  has_nodes_ = framework_mesh->hasNodes();
  has_node_faces_ = framework_mesh->hasNodeFaces();
  comm_ = framework_mesh_->getComm();
  gm_ = framework_mesh_->getGeometricModel();
  space_dim_ = framework_mesh_->getSpaceDimension();
  manifold_dim_ = framework_mesh_->getManifoldDimension();
  is_logical_ = framework_mesh_->isLogical();
  is_ordered_ = framework_mesh_->isOrdered();
  has_edges_ = framework_mesh_->hasEdges();
  ncells_owned = framework_mesh_->getNumEntities(Entity_kind::CELL, Parallel_kind::OWNED);
  ncells_all = framework_mesh_->getNumEntities(Entity_kind::CELL, Parallel_kind::ALL);
  nfaces_owned = framework_mesh_->getNumEntities(Entity_kind::FACE, Parallel_kind::OWNED);
  nfaces_all = framework_mesh_->getNumEntities(Entity_kind::FACE, Parallel_kind::ALL);
  if (has_edges_) {
    nedges_owned = framework_mesh_->getNumEntities(Entity_kind::EDGE, Parallel_kind::OWNED);
    nedges_all = framework_mesh_->getNumEntities(Entity_kind::EDGE, Parallel_kind::ALL);
  }
  if (has_nodes_) {
    nnodes_owned = framework_mesh_->getNumEntities(Entity_kind::NODE, Parallel_kind::OWNED);
    nnodes_all = framework_mesh_->getNumEntities(Entity_kind::NODE, Parallel_kind::ALL);
  }
  bool natural_ordered_maps = plist_->get<bool>("natural map ordering", false);
  maps_.initialize(*framework_mesh_, natural_ordered_maps);

  nboundary_faces_owned = maps_.getNBoundaryFaces(Parallel_kind::OWNED);
  nboundary_faces_all = maps_.getNBoundaryFaces(Parallel_kind::ALL);
  if (hasNodes()) {
    nboundary_nodes_owned = maps_.getNBoundaryNodes(Parallel_kind::OWNED);
    nboundary_nodes_all = maps_.getNBoundaryNodes(Parallel_kind::ALL);
  }

  cacheDefault(*this);
}


template <MemSpace_kind MEM>
bool
MeshCache<MEM>::isValidSetType(const AmanziGeometry::RegionType rtype, const Entity_kind kind) const
{
  if (rtype == AmanziGeometry::RegionType::LABELEDSET && framework_mesh_.get()) {
    return framework_mesh_->isValidSetType(rtype, kind);
  }
  return true;
}


template <MemSpace_kind MEM>
int
MeshCache<MEM>::getSetSize(const std::string& region_name,
                           const Entity_kind kind,
                           const Parallel_kind ptype) const
{
  return getSetEntities(region_name, kind, ptype).size();
}


//
// TODO --etc
// This should be updated -- only cache the ALL set, then construct the
// OWNED or GHOST set on demand.  No need to save both OWNED and ALL.
//
// TODO --etc
// Don't cache "ALL" or "BOUNDARY" labeled sets.  Are there others that are
// just as fast to create as to cache?
//
template <MemSpace_kind MEM>
typename MeshCache<MEM>::cEntity_ID_View
MeshCache<MEM>::getSetEntities(const std::string& region_name,
                               const Entity_kind kind,
                               const Parallel_kind ptype) const
{
  auto key = std::make_tuple(region_name, kind, ptype);
  auto key_all = std::make_tuple(region_name, kind, Parallel_kind::ALL);
  bool check_new = false;

  if (!sets_.count(key_all)) {
    auto region = getGeometricModel()->FindRegion(region_name);
    if (region == Teuchos::null) {
      Errors::Message msg;
      msg << "Cannot find region of name \"" << region_name << "\" in the geometric model.";
      Exceptions::amanzi_throw(msg);
    }

    auto this_on_host = onMemSpace<MemSpace_kind::HOST>(*this);
    sets_[key_all] = asDualView(resolveMeshSet(*region, kind, Parallel_kind::ALL, this_on_host));
    check_new = true;
  }

  if (!sets_.count(key)) {
    if (ptype == Parallel_kind::OWNED) {
      auto v_all = view<MemSpace_kind::HOST>(sets_.at(key_all));
      size_t n_ents = getNumEntities(kind, Parallel_kind::OWNED);
      size_t i;
      for (i = 0; i != v_all.size(); ++i)
        if (v_all(i) >= n_ents) break;
      auto v_owned = Kokkos::subview(v_all, Kokkos::make_pair((size_t)0, i));
      sets_[key] = asDualView(v_owned);

    } else if (ptype == Parallel_kind::GHOST) {
      auto v_all = view<MemSpace_kind::HOST>(sets_.at(key_all));
      size_t n_ents = getNumEntities(kind, Parallel_kind::OWNED);
      size_t i;
      for (i = 0; i != v_all.size(); ++i)
        if (v_all(i) >= n_ents) break;
      auto v_ghosted = Kokkos::subview(v_all, Kokkos::make_pair(i, v_all.size()));
      sets_[key] = asDualView(v_ghosted);
    } else {
      AMANZI_ASSERT(false);
    }
    check_new = true;
  }

  if (check_new) {
    //   auto region = getGeometricModel()->FindRegion(region_name);
    //   // Error on zero -- some zeros already error internally (at the framework
    //   // level) but others don't.  This is the highest level we can catch these at.
    //   int lsize = view<MEM>(sets_.at(key)).size();
    //   int gsize = 0;
    //   getComm()->SumAll(&lsize, &gsize, 1);
    //   if (gsize == 0 && getComm()->MyPID() == 0) {
    //     Errors::Message msg;
    //     msg << "AmanziMesh::getSetEntities: Region \"" << region->get_name() << "\" of type \""
    //         << to_string(region->get_type()) << "\" and kind \"" << to_string(kind)
    //         << "\" is empty (globally).\n";
    //     std::cout << msg.what();
    //     // Exceptions::amanzi_throw(msg);
    //   }
  }

  return view<MEM>(sets_.at(key));
}

template <MemSpace_kind MEM>
Kokkos::pair<typename MeshCache<MEM>::cEntity_ID_View, typename MeshCache<MEM>::cDouble_View>
MeshCache<MEM>::getSetEntitiesAndVolumeFractions(const std::string& region_name,
                                                 const Entity_kind kind,
                                                 const Parallel_kind ptype) const
{
  MeshSets::key_type key{ region_name, kind, ptype };
  if (!set_vol_fracs_.count(key)) {
    auto region = getGeometricModel()->FindRegion(region_name);
    if (region == Teuchos::null) {
      Errors::Message msg;
      msg << "Cannot find region of name \"" << region_name << "\" in the geometric model.";
      Exceptions::amanzi_throw(msg);
    }

    Double_View vol_fracs_list;
    const MeshCache<MemSpace_kind::HOST> this_on_host = onMemSpace<MemSpace_kind::HOST>(*this);
    sets_[key] =
      asDualView(resolveMeshSetVolumeFractions(*region, kind, ptype, vol_fracs_list, this_on_host));
    set_vol_fracs_[key] = asDualView<double>(vol_fracs_list);
  }
  return Kokkos::pair(view<MEM>(sets_.at(key)), view<MEM>(set_vol_fracs_.at(key)));
}


template <MemSpace_kind MEM>
Entity_ID
MeshCache<MEM>::getNumEntities(const Entity_kind kind, const Parallel_kind ptype) const
{
  Entity_ID nowned, nall;
  switch (kind) {
  case (Entity_kind::CELL):
    nowned = ncells_owned;
    nall = ncells_all;
    break;
  case (Entity_kind::FACE):
    nowned = nfaces_owned;
    nall = nfaces_all;
    break;
  case (Entity_kind::EDGE):
    nowned = nedges_owned;
    nall = nedges_all;
    break;
  case (Entity_kind::NODE):
    nowned = nnodes_owned;
    nall = nnodes_all;
    break;
  case (Entity_kind::BOUNDARY_FACE):
    nowned = nboundary_faces_owned;
    nall = nboundary_faces_all;
    break;
  case (Entity_kind::BOUNDARY_NODE):
    nowned = nboundary_nodes_owned;
    nall = nboundary_nodes_all;
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

template <MemSpace_kind MEM>
template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION AmanziGeometry::Point
MeshCache<MEM>::getCentroid(const Entity_kind kind, const Entity_ID ent) const
{
  switch (kind) {
  case (Entity_kind::CELL):
    return getCellCentroid<AP>(ent);
  case (Entity_kind::FACE):
    return getFaceCentroid<AP>(ent);
  case (Entity_kind::EDGE):
    return getEdgeCentroid<AP>(ent);
  case (Entity_kind::NODE):
    return getNodeCoordinate<AP>(ent);
  default:
    assert(false);
  }
  return AmanziGeometry::Point();
}


template <MemSpace_kind MEM>
template <Entity_kind EK, AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION AmanziGeometry::Point
MeshCache<MEM>::getCentroid(const Entity_ID ent) const
{
  if constexpr (EK == Entity_kind::CELL) {
    return getCellCentroid<AP>(ent);
  } else if constexpr (EK == Entity_kind::FACE) {
    return getFaceCentroid<AP>(ent);
  } else if constexpr (EK == Entity_kind::EDGE) {
    return getEdgeCentroid<AP>(ent);
  } else if constexpr (EK == Entity_kind::NODE) {
    return getNodeCoordinate<AP>(ent);
  }
  AMANZI_ASSERT(false);
  return AmanziGeometry::Point();
}


template <MemSpace_kind MEM>
template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION double
MeshCache<MEM>::getExtent(const Entity_kind kind, const Entity_ID ent) const
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


template <MemSpace_kind MEM>
template <Entity_kind EK, AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION double
MeshCache<MEM>::getExtent(const Entity_ID ent) const
{
  if constexpr (EK == Entity_kind::CELL) {
    return getCellVolume<AP>(ent);
  } else if constexpr (EK == Entity_kind::FACE) {
    return getFaceArea<AP>(ent);
  } else if constexpr (EK == Entity_kind::EDGE) {
    return getEdgeLength<AP>(ent);
  }
  AMANZI_ASSERT(false);
  return -1.0;
}


//===================
//    getFace*
//===================

template <MemSpace_kind MEM>
KOKKOS_INLINE_FUNCTION typename MeshCache<MEM>::cEntity_ID_View
MeshCache<MEM>::getFaceEdges(const Entity_ID f) const
{
  return Impl::RaggedGetter<MEM>::get(
    data_.face_edges_cached,
    data_.face_edges,
    framework_mesh_,
    [&](const int i) {
      cEntity_ID_View res;
      framework_mesh_->getFaceEdges(f, res);
      return res;
    },
    nullptr,
    f);
}

template <MemSpace_kind MEM>
KOKKOS_INLINE_FUNCTION const Entity_ID&
MeshCache<MEM>::getFaceEdge(const Entity_ID f, const size_type i) const
{
  assert(data_.face_edges_cached);
  return data_.face_edges.get<MEM>(f, i);
}

template <MemSpace_kind MEM>
KOKKOS_INLINE_FUNCTION
  Kokkos::pair<typename MeshCache<MEM>::cEntity_ID_View, typename MeshCache<MEM>::cDirection_View>
  MeshCache<MEM>::getFaceEdgesAndDirections(const Entity_ID f) const
{
  cEntity_ID_View edges;
  cDirection_View dirs;
  framework_mesh_->getFaceEdgesAndDirs(f, edges, &dirs);
  return Kokkos::make_pair(edges, dirs);
}

template <MemSpace_kind MEM>
KOKKOS_INLINE_FUNCTION void
MeshCache<MEM>::getFaceEdges(const Entity_ID f, cEntity_ID_View& edges) const
{
  auto [fedges, dirs] = getFaceEdgesAndDirections(f);
  edges = fedges;
}

template <MemSpace_kind MEM>
KOKKOS_INLINE_FUNCTION void
MeshCache<MEM>::getFaceEdgesAndDirs(const Entity_ID f,
                                    cEntity_ID_View& edges,
                                    cDirection_View* const dirs) const
{
  if constexpr (MEM == MemSpace_kind::DEVICE) {
    if (data_.face_edges_cached) {
      edges = data_.face_edges.getRowUnmanaged<MEM>(f);
      if (dirs) *dirs = data_.face_edge_directions.getRowUnmanaged<MEM>(f);
      return;
    }

  } else {
    if (data_.face_edges_cached) {
      edges = data_.face_edges.getRowUnmanaged<MEM>(f);
      if (dirs) *dirs = data_.face_edge_directions.getRowUnmanaged<MEM>(f);
      return;
    }

    if (framework_mesh_.get()) {
      framework_mesh_->getFaceEdgesAndDirs(f, edges, dirs);
      return;
    }
  }
  assert(false);
}


template <MemSpace_kind MEM>
template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION size_type
MeshCache<MEM>::getFaceNumCells(const Entity_ID f, const Parallel_kind ptype) const
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


template <MemSpace_kind MEM>
KOKKOS_INLINE_FUNCTION typename MeshCache<MEM>::cEntity_ID_View
MeshCache<MEM>::getFaceCells(const Entity_ID f) const
{
  cEntity_ID_View fcells;
  getFaceCells(f, fcells);
  return fcells;
}


template <MemSpace_kind MEM>
KOKKOS_INLINE_FUNCTION const Entity_ID&
MeshCache<MEM>::getFaceCell(const Entity_ID f, const size_type i) const
{
  assert(data_.face_cells_cached);
  return data_.face_cells.get<MEM>(f, i);
}

template <MemSpace_kind MEM>
template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION void
MeshCache<MEM>::getFaceCells(const Entity_ID f, cEntity_ID_View& fcells) const
{
  fcells = Impl::RaggedGetter<MEM, AP>::get(
    data_.face_cells_cached,
    data_.face_cells,
    framework_mesh_,
    [&](const int i) {
      cEntity_ID_View cells;
      framework_mesh_->getFaceCells(f, cells);
      return cells;
    },
    nullptr,
    f);
}


// Normal vector of a face
template <MemSpace_kind MEM>
template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION AmanziGeometry::Point
MeshCache<MEM>::getFaceNormal(const Entity_ID f) const
{
  return getFaceNormal<AP>(f, -1, nullptr);
}

template <MemSpace_kind MEM>
template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION AmanziGeometry::Point
MeshCache<MEM>::getFaceNormal(const Entity_ID f, const Entity_ID c, int* orientation) const
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

  if constexpr (MEM == MemSpace_kind::HOST) {
    return algorithms_->computeFaceNormal(*this, f, c, orientation);
  }

  assert(false && "No access to cache/framework/compute available");
  return AmanziGeometry::Point();
}


template <MemSpace_kind MEM>
template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION AmanziGeometry::Point
MeshCache<MEM>::getFaceCentroid(const Entity_ID f) const
{
  auto cf = Impl::ComputeFunction<MEM>::template hostOnly<
    std::function<AmanziGeometry::Point(const Entity_ID)>>();
  if constexpr (MEM == MemSpace_kind::HOST) {
    cf = [&](const int i) { return algorithms_->computeFaceCentroid(*this, i); };
  }

  return Impl::Getter<MEM, AP>::get(
    data_.face_geometry_cached, data_.face_centroids, framework_mesh_, nullptr, cf, f);
}


template <MemSpace_kind MEM>
KOKKOS_INLINE_FUNCTION typename MeshCache<MEM>::cEntity_ID_View
MeshCache<MEM>::getFaceNodes(const Entity_ID f) const
{
  cEntity_ID_View fcells;
  getFaceNodes(f, fcells);
  return fcells;
}


template <MemSpace_kind MEM>
KOKKOS_INLINE_FUNCTION const Entity_ID&
MeshCache<MEM>::getFaceNode(const Entity_ID f, const size_type i) const
{
  assert(data_.face_nodes_cached);
  return data_.face_nodes.get<MEM>(f, i);
}


template <MemSpace_kind MEM>
KOKKOS_INLINE_FUNCTION void
MeshCache<MEM>::getFaceNodes(const Entity_ID f, cEntity_ID_View& fcells) const
{
  if constexpr (MEM == MemSpace_kind::DEVICE) {
    if (data_.face_nodes_cached) {
      fcells = data_.face_nodes.getRowUnmanaged<MEM>(f);
      return;
    }

  } else {
    if (data_.face_nodes_cached) {
      fcells = data_.face_nodes.getRowUnmanaged<MEM>(f);
      return;
    }
    if (framework_mesh_.get()) {
      framework_mesh_->getFaceNodes(f, fcells);
      return;
    }
  }
  assert(false);
}

template <MemSpace_kind MEM>
template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION double
MeshCache<MEM>::getFaceArea(const Entity_ID f) const
{
  auto cf = Impl::ComputeFunction<MEM>::template hostOnly<std::function<double(const Entity_ID)>>();
  if constexpr (MEM == MemSpace_kind::HOST) {
    cf = [&](const Entity_ID i) { return algorithms_->computeFaceArea(*this, i); };
  }

  return Impl::Getter<MEM, AP>::get(
    data_.face_geometry_cached, data_.face_areas, framework_mesh_, nullptr, cf, f);
}

template <MemSpace_kind MEM>
template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION typename MeshCache<MEM>::cPoint_View
MeshCache<MEM>::getFaceCoordinates(const Entity_ID f) const
{
  return Impl::RaggedGetter<MEM, AP>::get(
    data_.face_coordinates_cached,
    data_.face_coordinates,
    framework_mesh_,
    [&](const Entity_ID i) { return framework_mesh_->getFaceCoordinates(i); },
    nullptr,
    f);
}

template <MemSpace_kind MEM>
KOKKOS_INLINE_FUNCTION std::vector<int>
MeshCache<MEM>::getFaceCellEdgeMap(const Entity_ID faceid, const Entity_ID cellid) const
{
  std::vector<int> map;
  cEntity_ID_View fedgeids;
  cDirection_View fedgedirs;

  getFaceEdgesAndDirs(faceid, fedgeids, &fedgedirs);
  auto cedgeids = getCellEdges(cellid);

  map.resize(fedgeids.size(), -1);

  for (int f = 0; f < fedgeids.size(); ++f) {
    Entity_ID fedge = fedgeids[f];

    for (int c = 0; c < cedgeids.size(); ++c) {
      if (fedge == cedgeids[c]) {
        map[f] = c;
        break;
      }
    }
  }
  return map;
}

//===================
//    getCell*
//===================

template <MemSpace_kind MEM>
KOKKOS_INLINE_FUNCTION std::size_t
MeshCache<MEM>::getCellMaxFaces() const
{
  std::size_t n(0);
  int ncells = getNumEntities(Entity_kind::CELL, Parallel_kind::OWNED);
  for (int c = 0; c < ncells; ++c) {
    auto v = getCellNumFaces(c);
    if (n < v) n = v;
  }
  return n;
}

template <MemSpace_kind MEM>
KOKKOS_INLINE_FUNCTION std::size_t
MeshCache<MEM>::getCellMaxEdges() const
{
  std::size_t n(0);
  if (hasEdges()) {
    int ncells = getNumEntities(Entity_kind::CELL, Parallel_kind::OWNED);
    for (int c = 0; c < ncells; ++c) {
      const auto& edges = getCellEdges(c);
      if (n < edges.size()) n = edges.size();
    }
  }
  return n;
}

// extent
template <MemSpace_kind MEM>
template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION double
MeshCache<MEM>::getCellVolume(const Entity_ID c) const
{
  auto cf = Impl::ComputeFunction<MEM>::template hostOnly<std::function<double(const Entity_ID)>>();
  if constexpr (MEM == MemSpace_kind::HOST) {
    cf = [&](const int i) { return algorithms_->computeCellVolume(*this, i); };
  }

  return Impl::Getter<MEM, AP>::get(
    data_.cell_geometry_cached, data_.cell_volumes, framework_mesh_, nullptr, cf, c);
}

template <MemSpace_kind MEM>
template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION size_type
MeshCache<MEM>::getCellNumFaces(const Entity_ID c) const
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

template <MemSpace_kind MEM>
template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION typename MeshCache<MEM>::cEntity_ID_View
MeshCache<MEM>::getCellFaces(const Entity_ID c) const
{
  cEntity_ID_View cfaces;
  getCellFaces<AP>(c, cfaces);
  return cfaces;
}


template <MemSpace_kind MEM>
template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION void
MeshCache<MEM>::getCellFaces(const Entity_ID c, cEntity_ID_View& cfaces) const
{
  cfaces = Impl::RaggedGetter<MEM, AP>::get(
    data_.cell_faces_cached,
    data_.cell_faces,
    framework_mesh_,
    [&](const int i) {
      cEntity_ID_View cf;
      framework_mesh_->getCellFaces(i, cf);
      return cf;
    },
    nullptr,
    c);
}


template <MemSpace_kind MEM>
KOKKOS_INLINE_FUNCTION const Entity_ID&
MeshCache<MEM>::getCellFace(const Entity_ID c, const size_type i) const
{
  assert(data_.cell_faces_cached);
  return data_.cell_faces.get<MEM>(c, i);
}


template <MemSpace_kind MEM>
KOKKOS_INLINE_FUNCTION
  Kokkos::pair<typename MeshCache<MEM>::cEntity_ID_View, typename MeshCache<MEM>::cDirection_View>
  MeshCache<MEM>::getCellFacesAndDirections(const Entity_ID c) const
{
  cEntity_ID_View cfaces;
  cDirection_View dirs;
  getCellFacesAndDirs(c, cfaces, &dirs);
  return Kokkos::pair(cfaces, dirs);
}


template <MemSpace_kind MEM>
KOKKOS_INLINE_FUNCTION void
MeshCache<MEM>::getCellFacesAndDirs(const Entity_ID c,
                                    cEntity_ID_View& faces,
                                    cDirection_View* const dirs) const
{
  if constexpr (MEM == MemSpace_kind::DEVICE) {
    if (data_.cell_faces_cached) {
      faces = data_.cell_faces.getRowUnmanaged<MEM>(c);
      if (dirs) *dirs = data_.cell_face_directions.getRowUnmanaged<MEM>(c);
      return;
    }

  } else {
    if (data_.cell_faces_cached) {
      faces = data_.cell_faces.getRowUnmanaged<MEM>(c);
      if (dirs) *dirs = data_.cell_face_directions.getRowUnmanaged<MEM>(c);
      return;
    }

    if (framework_mesh_.get()) {
      framework_mesh_->getCellFacesAndDirs(c, faces, dirs);
      return;
    }
  }
  assert(false);
}


template <MemSpace_kind MEM>
KOKKOS_INLINE_FUNCTION
  Kokkos::pair<typename MeshCache<MEM>::cEntity_ID_View, typename MeshCache<MEM>::cPoint_View>
  MeshCache<MEM>::getCellFacesAndBisectors(const Entity_ID c) const
{
  cEntity_ID_View cfaces;
  cPoint_View bisectors;
  getCellFacesAndBisectors(c, cfaces, &bisectors);
  return Kokkos::make_pair(cfaces, bisectors);
}


template <MemSpace_kind MEM>
KOKKOS_INLINE_FUNCTION void
MeshCache<MEM>::getCellFacesAndBisectors(const Entity_ID c,
                                         cEntity_ID_View& faces,
                                         cPoint_View* const bisectors) const
{
  if constexpr (MEM == MemSpace_kind::DEVICE) {
    if (data_.cell_faces_cached) {
      faces = data_.cell_faces.getRowUnmanaged<MEM>(c);
      if (bisectors) *bisectors = data_.cell_face_bisectors.getRowUnmanaged<MEM>(c);
      return;
    }

  } else {
    if (data_.cell_faces_cached) {
      faces = data_.cell_faces.getRowUnmanaged<MEM>(c);
      if (bisectors) *bisectors = data_.cell_face_bisectors.getRowUnmanaged<MEM>(c);
      return;
    }
    if (algorithms_.get()) {
      algorithms_->computeCellFacesAndBisectors(*this, c, faces, bisectors);
      return;
    }
  }
  assert(false);
}


template <MemSpace_kind MEM>
template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION AmanziGeometry::Point
MeshCache<MEM>::getCellCentroid(const Entity_ID c) const
{
  auto cf = Impl::ComputeFunction<MEM>::template hostOnly<
    std::function<AmanziGeometry::Point(const Entity_ID)>>();
  if constexpr (MEM == MemSpace_kind::HOST) {
    cf = [&](const int i) { return algorithms_->computeCellCentroid(*this, i); };
  }

  return Impl::Getter<MEM, AP>::get(
    data_.cell_geometry_cached, data_.cell_centroids, framework_mesh_, nullptr, cf, c);
}


template <MemSpace_kind MEM>
template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION size_type
MeshCache<MEM>::getCellNumNodes(const Entity_ID c) const
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


template <MemSpace_kind MEM>
KOKKOS_INLINE_FUNCTION typename MeshCache<MEM>::cEntity_ID_View
MeshCache<MEM>::getCellNodes(const Entity_ID c) const
{
  cEntity_ID_View nodes;
  getCellNodes(c, nodes);
  return nodes;
}

template <MemSpace_kind MEM>
KOKKOS_INLINE_FUNCTION Entity_ID
MeshCache<MEM>::getCellNode(const Entity_ID c, const size_type i) const
{
  // Compute list and use only one?
  cEntity_ID_View nodes;
  getCellNodes(c, nodes);
  return nodes[i];
}

template <MemSpace_kind MEM>
KOKKOS_INLINE_FUNCTION void
MeshCache<MEM>::getCellNodes(const Entity_ID c, cEntity_ID_View& nodes) const
{
  nodes = Impl::RaggedGetter<MEM, AccessPattern_kind::DEFAULT>::get(
    data_.cell_nodes_cached,
    data_.cell_nodes,
    framework_mesh_,
    [&](const int i) {
      cEntity_ID_View cf;
      framework_mesh_->getCellNodes(i, cf);
      return cf;
    },
    nullptr,
    c);

#if 0
  if constexpr (MEM == MemSpace_kind::DEVICE){
    assert(data_.cell_faces_cached);
    assert(data_.face_nodes_cached);
    for(int i = 0 ; i < nodes.size(); ++i)
      nodes[i] = -1;

    auto faces = getCellFaces(c);
    int i = 0;
    for(int f = 0 ; f < faces.size() ; ++f){
        auto fnodes = getFaceNodes(faces[f]);
        for(int n = 0; n < fnodes.size(); ++n){
          if(i >= nodes.size()){
            printf("Increase shared memory size\n");
            assert(false);
          }
          if(!is_present(fnodes[n],nodes)){
            nodes[i++] = fnodes[n];
          }
        }
      }
      nodes = nodes.subview(Kokkos::make_pair(0,i));
  }else{
    if constexpr(std::is_same_v<ViewType,Span<typename ViewType::value_type>>){
      auto v = Impl::computeCellNodes(*this,c);
      nodes = ViewType{v.data(),v.size()};
    } else {
      nodes = Impl::computeCellNodes(*this,c);
    }
  }
#endif
}

template <MemSpace_kind MEM>
template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION typename MeshCache<MEM>::cPoint_View
MeshCache<MEM>::getCellCoordinates(const Entity_ID c) const
{
  return Impl::RaggedGetter<MEM, AP>::get(
    data_.cell_coordinates_cached,
    data_.cell_coordinates,
    framework_mesh_,
    [&](const int i) { return framework_mesh_->getCellCoordinates(i); },
    nullptr,
    c);
}


template <MemSpace_kind MEM>
template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION size_type
MeshCache<MEM>::getCellNumEdges(const Entity_ID c) const
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


template <MemSpace_kind MEM>
template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION typename MeshCache<MEM>::cEntity_ID_View
MeshCache<MEM>::getCellEdges(const Entity_ID c) const
{
  return Impl::RaggedGetter<MEM, AP>::get(
    data_.cell_edges_cached,
    data_.cell_edges,
    framework_mesh_,
    [&](const int i) {
      cEntity_ID_View ce;
      framework_mesh_->getCellEdges(i, ce);
      return ce;
    },
    nullptr,
    c);
}


template <MemSpace_kind MEM>
template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION void
MeshCache<MEM>::getCellEdges(const Entity_ID c, cEntity_ID_View& cedges) const
{
  cedges = Impl::RaggedGetter<MEM, AP>::get(
    data_.cell_edges_cached,
    data_.cell_edges,
    framework_mesh_,
    [&](const int i) {
      cEntity_ID_View ce;
      framework_mesh_->getCellEdges(i, ce);
      return ce;
    },
    nullptr,
    c);
}


template <MemSpace_kind MEM>
KOKKOS_INLINE_FUNCTION const Entity_ID&
MeshCache<MEM>::getCellEdge(const Entity_ID c, const size_type i) const
{
  assert(data_.cell_edges_cached);
  return data_.cell_edges.get<MEM>(c, i);
}

template <MemSpace_kind MEM>
std::size_t
MeshCache<MEM>::getCellMaxNodes() const
{
  std::size_t n(0);
  int ncells = getNumEntities(Entity_kind::CELL, Parallel_kind::OWNED);
  for (int c = 0; c < ncells; ++c) {
    auto nodes = getCellNodes(c);
    if (n < nodes.size()) n = nodes.size();
  }
  return n;
}

//===================
//    getNode*
//===================

template <MemSpace_kind MEM>
template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION AmanziGeometry::Point
MeshCache<MEM>::getNodeCoordinate(const Entity_ID n) const
{
  return Impl::Getter<MEM, AP>::get(
    data_.node_coordinates_cached,
    data_.node_coordinates,
    framework_mesh_,
    [&](const int i) { return framework_mesh_->getNodeCoordinate(i); },
    nullptr,
    n);
}


template <MemSpace_kind MEM>
void
MeshCache<MEM>::setNodeCoordinate(const Entity_ID n, const AmanziGeometry::Point& p)
{
  if (data_.node_coordinates_cached) view<MemSpace_kind::HOST>(data_.node_coordinates)(n) = p;
  if (framework_mesh_.get()) framework_mesh_->setNodeCoordinate(n, p);
}

template <MemSpace_kind MEM>
void
MeshCache<MEM>::setNodeCoordinates(const cEntity_ID_View& nodes, const cPoint_View& new_coords)
{
  auto bf = view<MEM>(data_.node_coordinates);
  if (data_.node_coordinates_cached) {
    if constexpr (MEM == MemSpace_kind::HOST) {
      for (int i = 0; i != nodes.size(); ++i) { bf(nodes(i)) = new_coords(i); }
    } else {
      Kokkos::parallel_for(
        "setNodeCoordinates", nodes.size(), KOKKOS_LAMBDA(const int& i) {
          bf(nodes(i)) = new_coords(i);
        });
    }
  }

  if (framework_mesh_.get()) {
    View_type<const Entity_ID, MemSpace_kind::HOST> nodes_on_host;
    View_type<const AmanziGeometry::Point, MemSpace_kind::HOST> coords_on_host;
    if constexpr (MEM == MemSpace_kind::HOST) {
      nodes_on_host = nodes;
      coords_on_host = new_coords;
    } else {
      View_type<Entity_ID, MemSpace_kind::HOST> nc_nodes_on_host;
      View_type<AmanziGeometry::Point, MemSpace_kind::HOST> nc_coords_on_host;
      Kokkos::resize(nc_nodes_on_host, nodes.size());
      Kokkos::deep_copy(nc_nodes_on_host, nodes);
      nodes_on_host = nc_nodes_on_host;

      Kokkos::resize(nc_coords_on_host, nodes.size());
      Kokkos::deep_copy(nc_coords_on_host, new_coords);
      coords_on_host = nc_coords_on_host;
    }

    for (int i = 0; i != nodes.size(); ++i) {
      framework_mesh_->setNodeCoordinate(nodes_on_host(i), coords_on_host(i));
    }
  }
  recacheGeometry();
}


template <MemSpace_kind MEM>
template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION typename MeshCache<MEM>::cEntity_ID_View
MeshCache<MEM>::getNodeCells(const Entity_ID n, const Parallel_kind ptype) const
{
  cEntity_ID_View cells;
  getNodeCells<AP>(n, ptype, cells);
  return cells;
}

template <MemSpace_kind MEM>
template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION void
MeshCache<MEM>::getNodeCells(const Entity_ID n,
                             const Parallel_kind ptype,
                             cEntity_ID_View& cells) const
{
  cells = Impl::RaggedGetter<MEM, AP>::get(
    data_.node_cells_cached,
    data_.node_cells,
    framework_mesh_,
    [&](const int i) {
      cEntity_ID_View lcells;
      framework_mesh_->getNodeCells(i, lcells);
      return lcells;
    },
    nullptr,
    n);
}

template <MemSpace_kind MEM>
template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION typename MeshCache<MEM>::cEntity_ID_View
MeshCache<MEM>::getNodeFaces(const Entity_ID n) const
{
  cEntity_ID_View faces;
  getNodeFaces<AP>(n, faces);
  return faces;
}

template <MemSpace_kind MEM>
template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION void
MeshCache<MEM>::getNodeFaces(const Entity_ID n, cEntity_ID_View& faces) const
{
  faces = Impl::RaggedGetter<MEM, AP>::get(
    data_.node_faces_cached,
    data_.node_faces,
    framework_mesh_,
    [&](const int i) {
      cEntity_ID_View lfaces;
      framework_mesh_->getNodeFaces(i, lfaces);
      return lfaces;
    },
    nullptr,
    n);
}

template <MemSpace_kind MEM>
template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION typename MeshCache<MEM>::cEntity_ID_View
MeshCache<MEM>::getNodeEdges(const Entity_ID n) const
{
  cEntity_ID_View edges;
  getNodeEdges<AP>(n, edges);
  return edges;
}

template <MemSpace_kind MEM>
template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION void
MeshCache<MEM>::getNodeEdges(const Entity_ID n, cEntity_ID_View& edges) const
{
  edges = Impl::RaggedGetter<MEM, AP>::get(
    data_.node_edges_cached,
    data_.node_edges,
    framework_mesh_,
    [&](const int i) {
      cEntity_ID_View ledges;
      framework_mesh_->getNodeEdges(i, ledges);
      return ledges;
    },
    nullptr,
    n);
}

template <MemSpace_kind MEM>
KOKKOS_INLINE_FUNCTION Cell_kind
MeshCache<MEM>::getCellType(const Entity_ID c) const
{
  return Impl::getCellType(*this, c);
}

template <MemSpace_kind MEM>
Parallel_kind
MeshCache<MEM>::getParallelType(const Entity_kind& kind, const Entity_ID id) const
{
  if (id < getNumEntities(kind, Parallel_kind::OWNED)) {
    return Parallel_kind::OWNED;
  } else if (id < getNumEntities(kind, Parallel_kind::ALL)) {
    return Parallel_kind::GHOST;
  }
  return Parallel_kind::UNKNOWN;
}

//==================
//    getEdge*
//==================

template <MemSpace_kind MEM>
template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION typename MeshCache<MEM>::cEntity_ID_View
MeshCache<MEM>::getEdgeNodes(const Entity_ID e) const
{
  return Impl::RaggedGetter<MEM, AP>::get(
    data_.edge_nodes_cached,
    data_.edge_nodes,
    framework_mesh_,
    [&](const int i) {
      cEntity_ID_View nodes;
      framework_mesh_->getEdgeNodes(i, nodes);
      return nodes;
    },
    nullptr,
    e);
}

template <MemSpace_kind MEM>
template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION Entity_ID
MeshCache<MEM>::getEdgeNode(const Entity_ID e, const size_type i) const
{
  // Compute list and use only one?
  cEntity_ID_View nodes;
  getEdgeNodes(e, nodes);
  return nodes[i];
}


template <MemSpace_kind MEM>
template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION void
MeshCache<MEM>::getEdgeNodes(const Entity_ID e, cEntity_ID_View& nodes) const
{
  nodes = Impl::RaggedGetter<MEM, AP>::get(
    data_.edge_nodes_cached,
    data_.edge_nodes,
    framework_mesh_,
    [&](const int i) {
      cEntity_ID_View lnodes;
      framework_mesh_->getEdgeNodes(i, lnodes);
      return lnodes;
    },
    nullptr,
    e);
}

template <MemSpace_kind MEM>
template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION typename MeshCache<MEM>::cPoint_View
MeshCache<MEM>::getEdgeCoordinates(const Entity_ID c) const
{
  return Impl::RaggedGetter<MEM, AP>::get(
    data_.edge_coordinates_cached,
    data_.edge_coordinates,
    framework_mesh_,
    [&](const int i) { return framework_mesh_->getEdgeCoordinates(i); },
    nullptr,
    c);
}

template <MemSpace_kind MEM>
template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION AmanziGeometry::Point
MeshCache<MEM>::getEdgeCentroid(const Entity_ID c) const
{
  auto cf = Impl::ComputeFunction<MEM>::template hostOnly<
    std::function<AmanziGeometry::Point(const Entity_ID)>>();
  if constexpr (MEM == MemSpace_kind::HOST) {
    cf = [&](const int i) { return algorithms_->computeEdgeCentroid(*this, i); };
  }

  return Impl::Getter<MEM, AP>::get(
    data_.edge_geometry_cached, data_.edge_centroids, framework_mesh_, nullptr, cf, c);
}

template <MemSpace_kind MEM>
template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION AmanziGeometry::Point
MeshCache<MEM>::getEdgeVector(const Entity_ID e) const
{
  auto cf = Impl::ComputeFunction<MEM>::template hostOnly<
    std::function<AmanziGeometry::Point(const Entity_ID)>>();
  if constexpr (MEM == MemSpace_kind::HOST) {
    cf = [&](const int i) { return algorithms_->computeEdgeVector(*this, i, -1, nullptr); };
  }

  return Impl::Getter<MEM, AP>::get(
    data_.edge_geometry_cached, data_.edge_vectors, framework_mesh_, nullptr, cf, e);
}


template <MemSpace_kind MEM>
template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION double
MeshCache<MEM>::getEdgeLength(const Entity_ID e) const
{
  auto cf = Impl::ComputeFunction<MEM>::template hostOnly<std::function<double(const Entity_ID)>>();
  if constexpr (MEM == MemSpace_kind::HOST) {
    cf = [&](const int i) { return algorithms_->computeEdgeLength(*this, i); };
  }

  return Impl::Getter<MEM, AP>::get(
    data_.edge_lengths_cached, data_.edge_lengths, framework_mesh_, nullptr, cf, e);
}

template <MemSpace_kind MEM>
template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION typename MeshCache<MEM>::cEntity_ID_View
MeshCache<MEM>::getEdgeCells(const Entity_ID e) const
{
  cEntity_ID_View cells;
  getEdgeCells<AP>(e, cells);
  return cells;
}

template <MemSpace_kind MEM>
template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION void
MeshCache<MEM>::getEdgeCells(const Entity_ID e, cEntity_ID_View& cells) const
{
  cells = Impl::RaggedGetter<MEM, AP>::get(
    data_.edge_cells_cached,
    data_.edge_cells,
    framework_mesh_,
    [&](const int i) {
      cEntity_ID_View lcells;
      framework_mesh_->getEdgeCells(i, lcells);
      return lcells;
    },
    nullptr,
    e);
}

template <MemSpace_kind MEM>
template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION typename MeshCache<MEM>::cEntity_ID_View
MeshCache<MEM>::getEdgeFaces(const Entity_ID e) const
{
  cEntity_ID_View faces;
  getEdgeFaces<AP>(e, faces);
  return faces;
}

template <MemSpace_kind MEM>
template <AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION void
MeshCache<MEM>::getEdgeFaces(const Entity_ID e, cEntity_ID_View& faces) const
{
  faces = Impl::RaggedGetter<MEM, AP>::get(
    data_.edge_faces_cached,
    data_.edge_faces,
    framework_mesh_,
    [&](const int i) {
      cEntity_ID_View lfaces;
      framework_mesh_->getEdgeFaces(i, lfaces);
      return lfaces;
    },
    nullptr,
    e);
}

//
// Build the cache, fine grained control
// =============================================
template <MemSpace_kind MEM>
void
MeshCache<MEM>::cacheCellGeometry()
{
  assert(framework_mesh_.get());
  if (data_.cell_geometry_cached) return;
  auto lambda = [this](Entity_ID c, Double_View cvol, Point_View ccent) {
    std::tie(cvol[c], ccent[c]) = algorithms_->computeCellGeometry(*this, c);
  };
  std::tie(data_.cell_volumes, data_.cell_centroids) = asDualView(lambda, ncells_all);
  data_.cell_geometry_cached = true;
}

// cell-face adjacencies
template <MemSpace_kind MEM>
void
MeshCache<MEM>::cacheCellFaces()
{
  assert(framework_mesh_.get());
  if (data_.cell_faces_cached) return;
  int num_cells =
    framework_mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::ALL);
  auto lambda1 = [this](Entity_ID c, View_type<const Entity_ID, MemSpace_kind::HOST>& cfaces) {
    this->framework_mesh_->getCellFaces(c, cfaces);
  };
  data_.cell_faces = asRaggedArray_DualView<Entity_ID>(lambda1, num_cells);
  auto lambda2 = [this](Entity_ID c, View_type<const Direction_type, MemSpace_kind::HOST>& dirs) {
    this->framework_mesh_->getCellFaceDirs(c, dirs);
  };
  data_.cell_face_directions = asRaggedArray_DualView<int>(lambda2, num_cells);
  data_.cell_faces_cached = true;
}

template <MemSpace_kind MEM>
void
MeshCache<MEM>::cacheCellCoordinates()
{
  assert(framework_mesh_.get());
  if (data_.cell_coordinates_cached) return;
  int num_cells =
    framework_mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::ALL);
  auto lambda = [this](Entity_ID c,
                       View_type<const AmanziGeometry::Point, MemSpace_kind::HOST>& ccoords) {
    ccoords = this->framework_mesh_->getCellCoordinates(c);
  };
  data_.cell_coordinates = asRaggedArray_DualView<Amanzi::AmanziGeometry::Point>(lambda, num_cells);
  data_.cell_coordinates_cached = true;
}

template <MemSpace_kind MEM>
void
MeshCache<MEM>::cacheCellEdges()
{
  assert(framework_mesh_.get());
  if (data_.cell_edges_cached) return;
  framework_mesh_->hasEdgesOrThrow();
  int num_cells =
    framework_mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::ALL);
  auto lambda = [this](Entity_ID c, View_type<const Entity_ID, MemSpace_kind::HOST>& cedges) {
    this->framework_mesh_->getCellEdges(c, cedges);
  };
  data_.cell_edges = asRaggedArray_DualView<Entity_ID>(lambda, num_cells);
  data_.cell_edges_cached = true;
}

template <MemSpace_kind MEM>
void
MeshCache<MEM>::cacheCellNodes()
{
  assert(framework_mesh_.get());
  if (data_.cell_nodes_cached) return;
  int num_cells =
    framework_mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::ALL);
  auto lambda = [this](Entity_ID c, View_type<const Entity_ID, MemSpace_kind::HOST>& cnodes) {
    this->framework_mesh_->getCellNodes(c, cnodes);
  };
  data_.cell_nodes = asRaggedArray_DualView<Entity_ID>(lambda, num_cells);
  data_.cell_nodes_cached = true;
}

template <MemSpace_kind MEM>
void
MeshCache<MEM>::cacheFaceGeometry()
{
  if (data_.face_geometry_cached) return;
  assert(framework_mesh_.get());
  data_.face_areas.resize(nfaces_all);
  data_.face_centroids.resize(nfaces_all);

  // slurp down the RaggedArray for normals using a lambda that, as a side
  // effect, captures area and centroid too.
  auto area_view = view<MemSpace_kind::HOST>(data_.face_areas);
  auto centroid_view = view<MemSpace_kind::HOST>(data_.face_centroids);
  auto lambda = [&, this](const Entity_ID& f, cPoint_View& normals) {
    auto area_cent_normal = algorithms_->computeFaceGeometry(*this, f);
    area_view[f] = std::get<0>(area_cent_normal);
    centroid_view[f] = std::get<1>(area_cent_normal);
    normals = std::get<2>(area_cent_normal);
  };
  data_.face_normals = asRaggedArray_DualView<AmanziGeometry::Point>(lambda, nfaces_all);

  // still must sync areas/centroids
  Kokkos::deep_copy(view<MemSpace_kind::DEVICE>(data_.face_areas),
                    view<MemSpace_kind::HOST>(data_.face_areas));
  Kokkos::deep_copy(view<MemSpace_kind::DEVICE>(data_.face_centroids),
                    view<MemSpace_kind::HOST>(data_.face_centroids));
  data_.face_geometry_cached = true;

  // cache normal orientations
  //
  // This is complicated... it encodes two pieces of info:
  //   * a sign, which says whether the face's natural normal is inward or outward relative to the cell
  //   * a magnitude, which is 1 if the natural normal is defined as
  //     positive/negative of this or a 2 if the natural normal is defined as
  //     an average of this and the other normal.  Note 2 is only used when
  //     space_dim != manifold_dim and ncells == 2
  auto lambda2 = [&, this](const Entity_ID& f,
                           View_type<const Direction_type, MemSpace_kind::HOST>& dirs) {
    Direction_View ldirs;
    // This NEEDS to call the framework or be passed an host mesh to call the function on the host.
    View_type<const Entity_ID, MemSpace_kind::HOST> fcells;
    framework_mesh_->getFaceCells(f, fcells);
    Kokkos::resize(ldirs, fcells.size());
    for (int i = 0; i != fcells.size(); ++i) {
      if ((getSpaceDimension() == getManifoldDimension()) || (fcells.size() != 2)) {
        ldirs(i) = Impl::getFaceDirectionInCell(*this, f, fcells(i));
      } else {
        ldirs(i) = 2 * Impl::getFaceDirectionInCell(*this, f, fcells(i));
      }
    }
    dirs = ldirs;
  };
  data_.face_normal_orientations = asRaggedArray_DualView<int>(lambda2, nfaces_all);

  // cache cell-face-bisectors -- make this a separate call?  Think about
  // granularity here.
  auto lambda3 = [&, this](const Entity_ID& c, cPoint_View& bisectors) {
    View_type<const Entity_ID, MemSpace_kind::HOST> cfaces;
    algorithms_->computeCellFacesAndBisectors(*this, c, cfaces, &bisectors);
  };
  data_.cell_face_bisectors = asRaggedArray_DualView<AmanziGeometry::Point>(lambda3, ncells_all);
}


// face-cell adjacencies
template <MemSpace_kind MEM>
void
MeshCache<MEM>::cacheFaceCells()
{
  if (data_.face_cells_cached) return;
  auto lambda = [this](Entity_ID f, View_type<const Entity_ID, MemSpace_kind::HOST>& fcells) {
    this->framework_mesh_->getFaceCells(f, fcells);
  };
  data_.face_cells = asRaggedArray_DualView<Entity_ID>(lambda, nfaces_all);
  data_.face_cells_cached = true;
}

// face-cell adjacencies
template <MemSpace_kind MEM>
void
MeshCache<MEM>::cacheFaceEdges()
{
  if (data_.face_edges_cached) return;
  framework_mesh_->hasEdgesOrThrow();
  auto lambda = [this](Entity_ID f, View_type<const Entity_ID, MemSpace_kind::HOST>& fedges) {
    this->framework_mesh_->getFaceEdges(f, fedges);
  };
  auto lambda_dir = [this](Entity_ID f,
                           View_type<const Direction_type, MemSpace_kind::HOST>& dirs) {
    View_type<const Entity_ID, MemSpace_kind::HOST> l;
    framework_mesh_->getFaceEdgesAndDirs(f, l, &dirs);
  };
  data_.face_edges = asRaggedArray_DualView<Entity_ID>(lambda, nfaces_all);
  data_.face_edge_directions = asRaggedArray_DualView<int>(lambda_dir, nfaces_all);
  data_.face_edges_cached = true;
}

// face-node adjacencies
template <MemSpace_kind MEM>
void
MeshCache<MEM>::cacheFaceNodes()
{
  if (data_.face_nodes_cached) return;
  auto lambda = [this](Entity_ID f, View_type<const Entity_ID, MemSpace_kind::HOST>& fnodes) {
    this->framework_mesh_->getFaceNodes(f, fnodes);
  };
  data_.face_nodes = asRaggedArray_DualView<Entity_ID>(lambda, nfaces_all);
  data_.face_nodes_cached = true;
}

template <MemSpace_kind MEM>
void
MeshCache<MEM>::cacheFaceCoordinates()
{
  assert(framework_mesh_.get());
  if (data_.face_coordinates_cached) return;
  auto lambda = [this](Entity_ID f, cPoint_View& fcoords) {
    fcoords = this->framework_mesh_->getFaceCoordinates(f);
  };
  data_.face_coordinates =
    asRaggedArray_DualView<Amanzi::AmanziGeometry::Point>(lambda, nfaces_all);

  data_.face_coordinates_cached = true;
}

// edge centroid, length, vector
template <MemSpace_kind MEM>
void
MeshCache<MEM>::cacheEdgeGeometry()
{
  assert(framework_mesh_.get());
  if (data_.edge_geometry_cached) return;
  framework_mesh_->hasEdgesOrThrow();
  auto lambda = [this](Entity_ID e, Point_View ev, Point_View ec) {
    std::tie(ev[e], ec[e]) = algorithms_->computeEdgeGeometry(*this, e);
  };
  std::tie(data_.edge_vectors, data_.edge_centroids) = asDualView(lambda, nedges_all);
  data_.edge_geometry_cached = true;
}

// edge-cell adjacencies
template <MemSpace_kind MEM>
void
MeshCache<MEM>::cacheEdgeCells()
{
  assert(framework_mesh_.get());
  if (data_.edge_cells_cached) return;
  framework_mesh_->hasEdgesOrThrow();
  auto lambda = [this](Entity_ID f, View_type<const Entity_ID, MemSpace_kind::HOST>& ecells) {
    this->framework_mesh_->getEdgeCells(f, ecells);
  };
  data_.edge_cells = asRaggedArray_DualView<Entity_ID>(lambda, nedges_all);
  data_.edge_cells_cached = true;
}

// edge-face adjacencies
template <MemSpace_kind MEM>
void
MeshCache<MEM>::cacheEdgeFaces()
{
  assert(framework_mesh_.get());
  if (data_.edge_faces_cached) return;
  framework_mesh_->hasEdgesOrThrow();
  auto lambda = [this](Entity_ID f, View_type<const Entity_ID, MemSpace_kind::HOST>& efaces) {
    this->framework_mesh_->getEdgeFaces(f, efaces);
  };
  data_.edge_faces = asRaggedArray_DualView<Entity_ID>(lambda, nedges_all);
  data_.edge_faces_cached = true;
}

// edge-node adjacencies
template <MemSpace_kind MEM>
void
MeshCache<MEM>::cacheEdgeNodes()
{
  assert(framework_mesh_.get());
  if (data_.edge_nodes_cached) return;
  framework_mesh_->hasEdgesOrThrow();
  auto lambda = [this](Entity_ID f, View_type<const Entity_ID, MemSpace_kind::HOST>& efaces) {
    this->framework_mesh_->getEdgeNodes(f, efaces);
  };
  data_.edge_nodes = asRaggedArray_DualView<Entity_ID>(lambda, nedges_all);
  data_.edge_nodes_cached = true;
}

// edge coordinates
template <MemSpace_kind MEM>
void
MeshCache<MEM>::cacheEdgeCoordinates()
{
  assert(framework_mesh_.get());
  if (data_.edge_coordinates_cached) return;
  framework_mesh_->hasEdgesOrThrow();
  auto lambda = [this](Entity_ID f, cPoint_View& enodes) {
    enodes = this->framework_mesh_->getEdgeCoordinates(f);
  };
  data_.edge_coordinates =
    asRaggedArray_DualView<Amanzi::AmanziGeometry::Point>(lambda, nedges_all);
  data_.edge_coordinates_cached = true;
}

// // node-cell adjacencies
template <MemSpace_kind MEM>
void
MeshCache<MEM>::cacheNodeCells()
{
  assert(framework_mesh_.get());
  if (data_.node_cells_cached) return;
  auto lambda = [this](Entity_ID f, View_type<const Entity_ID, MemSpace_kind::HOST>& ncells) {
    this->framework_mesh_->getNodeCells(f, ncells);
  };
  data_.node_cells = asRaggedArray_DualView<Entity_ID>(lambda, nnodes_all);
  data_.node_cells_cached = true;
}

// // node-face adjacencies
template <MemSpace_kind MEM>
void
MeshCache<MEM>::cacheNodeFaces()
{
  assert(framework_mesh_.get());
  if (data_.node_faces_cached) return;
  auto lambda = [this](Entity_ID f, View_type<const Entity_ID, MemSpace_kind::HOST>& nfaces) {
    this->framework_mesh_->getNodeFaces(f, nfaces);
  };
  data_.node_faces = asRaggedArray_DualView<Entity_ID>(lambda, nnodes_all);
  data_.node_faces_cached = true;
}

// // node-edge adjacencies
template <MemSpace_kind MEM>
void
MeshCache<MEM>::cacheNodeEdges()
{
  assert(framework_mesh_.get());
  if (data_.node_edges_cached) return;
  framework_mesh_->hasEdgesOrThrow();
  auto lambda = [this](Entity_ID f, View_type<const Entity_ID, MemSpace_kind::HOST>& nedges) {
    this->framework_mesh_->getNodeEdges(f, nedges);
  };
  data_.node_edges = asRaggedArray_DualView<Entity_ID>(lambda, nnodes_all);
  data_.node_edges_cached = true;
}


// node coordinates
template <MemSpace_kind MEM>
void
MeshCache<MEM>::cacheNodeCoordinates()
{
  if (data_.node_coordinates_cached) return;
  data_.node_coordinates.resize(nnodes_all);
  for (Entity_ID i = 0; i != nnodes_all; ++i) {
    view<MemSpace_kind::HOST>(data_.node_coordinates)[i] = framework_mesh_->getNodeCoordinate(i);
  }

  Kokkos::deep_copy(data_.node_coordinates.d_view, data_.node_coordinates.h_view);
  data_.node_coordinates_cached = true;
}


// parent entities
template <MemSpace_kind MEM>
void
MeshCache<MEM>::cacheParentEntities()
{
  if (data_.parent_entities_cached) return;
  assert(framework_mesh_.get());

  data_.parent_cells.resize(ncells_all);
  for (Entity_ID i = 0; i != ncells_all; ++i) {
    view<MemSpace_kind::HOST>(data_.parent_cells)[i] =
      framework_mesh_->getEntityParent(AmanziMesh::Entity_kind::CELL, i);
  }
  Kokkos::deep_copy(data_.parent_cells.d_view, data_.parent_cells.h_view);

  data_.parent_faces.resize(nfaces_all);
  for (Entity_ID i = 0; i != nfaces_all; ++i) {
    view<MemSpace_kind::HOST>(data_.parent_faces)[i] =
      framework_mesh_->getEntityParent(AmanziMesh::Entity_kind::FACE, i);
  }
  Kokkos::deep_copy(data_.parent_faces.d_view, data_.parent_faces.h_view);

  data_.parent_nodes.resize(nnodes_all);
  for (Entity_ID i = 0; i != nnodes_all; ++i) {
    view<MemSpace_kind::HOST>(data_.parent_nodes)[i] =
      framework_mesh_->getEntityParent(AmanziMesh::Entity_kind::NODE, i);
  }
  Kokkos::deep_copy(data_.parent_nodes.d_view, data_.parent_nodes.h_view);

  data_.parent_entities_cached = true;
}


// // Note that regions are cached on demand the first time they are requested,
// // but labeled sets must be pre-cached if the framework mesh is to be
// // destroyed.
// void MeshCache<MEM>precacheLabeledSets()
// {
//   for (const auto& rgn : *getGeometricModel()) {
//     if (rgn->get_type() == AmanziGeometry::RegionType::LABELEDSET) {
//       auto rgn_lbl = Teuchos::rcp_dynamic_cast<const AmanziGeometry::RegionLabeledSet>(rgn);
//       AMANZI_ASSERT(rgn_lbl);

//       if (getParentMesh() != Teuchos::null) {
//         AMANZI_ASSERT(false); // not yet implemented lifted sets
//       } else {
//         auto entity_kind = createEntityKind(rgn_lbl->entity_str());
//         Entity_ID_View set_ents;
//         framework_mesh_->getSetEntities(*rgn_lbl, entity_kind, Parallel_kind::ALL, set_ents);
//         auto key = std::make_tuple(rgn->get_name(), entity_kind, Parallel_kind::ALL);
//         sets_[key] = asDualView(set_ents);
//       }
//     }
//   }
// }

// build the MeshMaps object
// void MeshCache<MEM>::buildMaps_(bool natural)
// {
// }


// // common error messaging
// void MeshCache<MEM>throwAccessError_(const std::string& func_name) const
// {
//   Errors::Message msg;
//   msg << "MeshCache<MEM>" << func_name << " cannot compute this quantity -- not cached and framework does not exist.";
//   Exceptions::amanzi_throw(msg);
// }

template <MemSpace_kind MEM>
Entity_ID
MeshCache<MEM>::getEntityParent(const Entity_kind kind, const Entity_ID entid) const
{
  if (data_.parent_entities_cached) {
    switch (kind) {
    case Entity_kind::CELL:
      return view<MEM>(data_.parent_cells)[entid];
      break;
    case Entity_kind::FACE:
      return view<MEM>(data_.parent_faces)[entid];
      break;
    case Entity_kind::EDGE:
      return view<MEM>(data_.parent_edges)[entid];
      break;
    case Entity_kind::NODE:
      return view<MEM>(data_.parent_nodes)[entid];
    default: {
    }
    }
  } else if (framework_mesh_.get()) {
    return framework_mesh_->getEntityParent(kind, entid);
  }
  return -1;
}

template <MemSpace_kind MEM>
typename MeshCache<MEM>::cEntity_ID_View
MeshCache<MEM>::getEntityParents(const Entity_kind kind) const
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
  default: {
  }
  }
  return MeshCache<MEM>::cEntity_ID_View();
}

template <MemSpace_kind MEM>
bool
MeshCache<MEM>::isPointInCell(const AmanziGeometry::Point& p, const Entity_ID cellid) const
{
  if (manifold_dim_ == 3) {
    // 3D Elements with possibly curved faces
    // We have to build a description of the element topology
    // and send it into the polyhedron volume and centroid
    // calculation routine
    int nf;
    std::vector<std::size_t> nfnodes;
    std::vector<AmanziGeometry::Point> cfcoords;

    auto [faces, fdirs] = getCellFacesAndDirections(cellid);

    nf = faces.size();
    nfnodes.resize(nf);

    for (int j = 0; j < nf; j++) {
      auto fcoords = getFaceCoordinates(faces[j]);
      nfnodes[j] = fcoords.size();

      if (fdirs[j] == 1) {
        for (int k = 0; k < nfnodes[j]; k++) cfcoords.push_back(fcoords[k]);
      } else {
        for (int k = nfnodes[j] - 1; k >= 0; k--) cfcoords.push_back(fcoords[k]);
      }
    }

    auto ccoords = getCellCoordinates(cellid);
    return AmanziGeometry::point_in_polyhed(p, ccoords, nf, nfnodes, cfcoords);

  } else if (manifold_dim_ == 2) {
    auto ccoords = getCellCoordinates(cellid);
    return AmanziGeometry::point_in_polygon(p, ccoords);
  }

  return false;
}


template <MemSpace_kind MEM>
void
MeshCache<MEM>::PrintMeshStatistics() const
{
  auto vo_ = Teuchos::rcp(new VerboseObject("Mesh Output", *plist_));
  if (vo_.get() && vo_->getVerbLevel() >= Teuchos::VERB_LOW) {
    int ncells = getNumEntities(AmanziMesh::CELL, AmanziMesh::Parallel_kind::OWNED);
    int nfaces = getNumEntities(AmanziMesh::FACE, AmanziMesh::Parallel_kind::OWNED);
    int nnodes = getNumEntities(AmanziMesh::NODE, AmanziMesh::Parallel_kind::OWNED);
    int nedges(0);
    if (has_edges_) nedges = getNumEntities(AmanziMesh::EDGE, AmanziMesh::Parallel_kind::OWNED);

    int min_out[4], max_out[4], sum_out[4], tmp_in[4] = { ncells, nfaces, nedges, nnodes };
    getComm()->MinAll(tmp_in, min_out, 4);
    getComm()->MaxAll(tmp_in, max_out, 4);
    getComm()->SumAll(tmp_in, sum_out, 4);

    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "cells, tot/min/max: " << sum_out[0] << "/" << min_out[0] << "/" << max_out[0]
               << "\n";
    *vo_->os() << "faces, tot/min/max: " << sum_out[1] << "/" << min_out[1] << "/" << max_out[1]
               << "\n";
    if (has_edges_)
      *vo_->os() << "edges, tot/min/max: " << sum_out[2] << "/" << min_out[2] << "/" << max_out[2]
                 << "\n";
    *vo_->os() << "nodes, tot/min/max: " << sum_out[3] << "/" << min_out[3] << "/" << max_out[3]
               << "\n\n";
  }
}


template <MemSpace_kind MEM>
void
MeshCache<MEM>::recacheGeometry()
{
  // recaches the geometry, as presumably the nodal coordinates have changed.
  if (data_.edge_coordinates_cached) {
    data_.edge_coordinates_cached = false;
    cacheEdgeCoordinates();
  }
  if (data_.edge_geometry_cached) {
    data_.edge_geometry_cached = false;
    cacheEdgeGeometry();
  }
  if (data_.face_coordinates_cached) {
    data_.face_coordinates_cached = false;
    cacheFaceCoordinates();
  }
  if (data_.face_geometry_cached) {
    data_.face_geometry_cached = false;
    cacheFaceGeometry();
  }
  if (data_.cell_coordinates_cached) {
    data_.cell_coordinates_cached = false;
    cacheCellCoordinates();
  }
  if (data_.cell_geometry_cached) {
    data_.cell_geometry_cached = false;
    cacheCellGeometry();
  }
}


template <MemSpace_kind MEM>
void
cacheAll(MeshCache<MEM>& mesh)
{
  mesh.cacheCellNodes();
  mesh.cacheCellCoordinates();
  mesh.cacheFaceCoordinates();
  mesh.cacheNodeCells();
  mesh.cacheNodeFaces();
  if (mesh.hasEdges()) {
    mesh.cacheCellEdges();
    mesh.cacheEdgeCells();
    mesh.cacheNodeEdges();
    mesh.cacheEdgeNodes();
    mesh.cacheEdgeCoordinates();
  }
}


template <MemSpace_kind MEM>
void
cacheDefault(MeshCache<MEM>& mesh)
{
  if (mesh.hasNodes()) { mesh.cacheNodeCoordinates(); }
  mesh.cacheCellFaces();
  mesh.cacheFaceCells();
  if (mesh.hasNodes()) { mesh.cacheFaceNodes(); }
  mesh.cacheCellGeometry();
  mesh.cacheFaceGeometry();
  if (mesh.hasEdges()) {
    mesh.cacheFaceEdges();
    mesh.cacheEdgeFaces();
    mesh.cacheEdgeGeometry();
  }
}


// -----------------------------------------------------------------------------
// Memory space transfer
// -----------------------------------------------------------------------------
template <MemSpace_kind OUT, MemSpace_kind IN>
inline const MeshCache<OUT>
onMemSpace(const MeshCache<IN>& mc_in)
{
  if constexpr (IN == OUT) return mc_in;
  MeshCache<IN>* mc_in_nc = const_cast<MeshCache<IN>*>(&mc_in);
  MeshCache<OUT> mc_out(*mc_in_nc);
  return mc_out;
}


template <MemSpace_kind OUT, MemSpace_kind IN>
inline Teuchos::RCP<const MeshCache<OUT>>
onMemSpace(const Teuchos::RCP<const MeshCache<IN>>& mc_in)
{
  if constexpr (IN == OUT) return mc_in;
  MeshCache<IN>* mc_in_nc = const_cast<MeshCache<IN>*>(mc_in.ptr().get());
  Teuchos::RCP<MeshCache<OUT>> mc_on_out = Teuchos::rcp(new MeshCache<OUT>(*mc_in_nc));
  return (Teuchos::RCP<const MeshCache<OUT>>)mc_on_out;
}


template <MemSpace_kind OUT, MemSpace_kind IN>
inline Teuchos::RCP<MeshCache<OUT>>
onMemSpace(const Teuchos::RCP<MeshCache<IN>>& mc_in)
{
  if constexpr (IN == OUT) return mc_in;
  return Teuchos::rcp(new MeshCache<OUT>(*mc_in));
}


} // namespace AmanziMesh
} // namespace Amanzi
