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

#include "MeshUtils.hh"
#include "MeshFramework.hh"

namespace Amanzi {
namespace AmanziMesh {

template<MemSpace_kind MEM>
MeshCache<MEM>::MeshCache()
  : is_ordered_(false),
    has_edges_(false),
    has_nodes_(true)
{}


template<MemSpace_kind MEM>
void MeshCache<MEM>::setParentMesh(const Teuchos::RCP<const MeshCache>& parent)
{
  if (parent_ != Teuchos::null && parent_ != parent) {
    Errors::Message msg("MeshCache::setParentMesh given conflicting parent mesh.");
    Exceptions::amanzi_throw(msg);
  }
  parent_ = parent;
}

template<MemSpace_kind MEM>
void MeshCache<MEM>::setMeshFramework(const Teuchos::RCP<MeshFramework>& framework_mesh) {
  framework_mesh_ = framework_mesh;
  // always save the algorithms, so we can throw away the data
  algorithms_ = framework_mesh->getAlgorithms();
  has_edges_ = framework_mesh->hasEdges();
  has_nodes_ = framework_mesh->hasNodes();
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
  if(hasNodes()){
    nboundary_nodes_owned = maps_.getNBoundaryNodes(Parallel_kind::OWNED);
    nboundary_nodes_all = maps_.getNBoundaryNodes(Parallel_kind::ALL);
  }
}


template<MemSpace_kind MEM>
bool
MeshCache<MEM>::isValidSetType(const AmanziGeometry::RegionType rtype, const Entity_kind kind) const {
  if (rtype == AmanziGeometry::RegionType::LABELEDSET && framework_mesh_.get()) {
    return framework_mesh_->isValidSetType(rtype, kind);
  }
  return true;
}


template<MemSpace_kind MEM>
int
MeshCache<MEM>::getSetSize(const std::string& region_name,
                      const Entity_kind kind,
                      const Parallel_kind ptype) const {
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
template<MemSpace_kind MEM>
decltype(auto)
MeshCache<MEM>::getSetEntities(const std::string& region_name,
        const Entity_kind kind,
        const Parallel_kind ptype) const
{
  auto key = std::make_tuple(region_name, kind, ptype);
  if (!sets_.count(key)) {
    auto region = getGeometricModel()->FindRegion(region_name);
    if (region == Teuchos::null) {
      Errors::Message msg;
      msg << "Cannot find region of name \"" << region_name << "\" in the geometric model.";
      Exceptions::amanzi_throw(msg);
    }
    MeshCache<MemSpace_kind::HOST> this_on_host(*this);
    sets_[key] =  asDualView(resolveMeshSet(*region, kind, ptype, this_on_host));
    // Error on zero -- some zeros already error internally (at the framework
    // level) but others don't.  This is the highest level we can catch these at.
    int lsize = view<MEM>(sets_.at(key)).size();
    int gsize = 0;
    getComm()->SumAll(&lsize, &gsize, 1);
    if (gsize == 0) {
      Errors::Message msg;
      msg << "AmanziMesh::getSetEntities: Region \"" << region->get_name() << "\" of type \"" 
          << to_string(region->get_type()) << "\" is empty (globally).";
      Exceptions::amanzi_throw(msg);
    }
  }
  return view<MEM>(sets_.at(key)); 
}

template<MemSpace_kind MEM>
decltype(auto)
MeshCache<MEM>::getSetEntitiesAndVolumeFractions(const std::string& region_name,
        const Entity_kind kind,
        const Parallel_kind ptype) const
{
  auto key = std::make_tuple(region_name, kind, ptype);
  if (!set_vol_fracs_.count(key)) {
    auto region = getGeometricModel()->FindRegion(region_name);
    if (region == Teuchos::null) {
      Errors::Message msg;
      msg << "Cannot find region of name \"" << region_name << "\" in the geometric model.";
      Exceptions::amanzi_throw(msg);
    }

    Double_View vol_fracs_list;
    MeshCache<MemSpace_kind::HOST> this_on_host(*this);
    sets_[key] = asDualView(resolveMeshSetVolumeFractions(*region, kind, ptype, vol_fracs_list, this_on_host));
    set_vol_fracs_[key] = asDualView<double>(vol_fracs_list);
  }
  return std::tie(view<MEM>(sets_.at(key)), view<MEM>(set_vol_fracs_.at(key)));
}


template<MemSpace_kind MEM>
Entity_ID
MeshCache<MEM>::getNumEntities(const Entity_kind kind, const Parallel_kind ptype) const
{
  Entity_ID nowned, nall;
  switch(kind) {
    case (Entity_kind::CELL) :
      nowned = ncells_owned; nall = ncells_all;
      break;
    case (Entity_kind::FACE) :
      nowned = nfaces_owned; nall = nfaces_all;
      break;
    case (Entity_kind::EDGE) :
      nowned = nedges_owned; nall = nedges_all;
      break;
    case (Entity_kind::NODE) :
      nowned = nnodes_owned; nall = nnodes_all;
      break;
    case (Entity_kind::BOUNDARY_FACE) :
      nowned = nboundary_faces_owned; nall = nboundary_faces_all;
      break;
    case (Entity_kind::BOUNDARY_NODE) :
      nowned = nboundary_nodes_owned; nall = nboundary_nodes_all;
      break;
    default :
      nowned = -1; nall = -1;
  }

  switch(ptype) {
    case (Parallel_kind::OWNED) :
      return nowned;
      break;
    case (Parallel_kind::ALL) :
      return nall;
      break;
    case Parallel_kind::GHOST :
      return nall - nowned;
      break;
    default :
      return 0;
  }
}

template<MemSpace_kind MEM> 
template<AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION
decltype(auto) MeshCache<MEM>::getCentroid(const Entity_kind kind, const Entity_ID ent) const {
  switch(kind) {
    case (Entity_kind::CELL) :
      return getCellCentroid<AP>(ent); 
    case (Entity_kind::FACE) :
      return getFaceCentroid<AP>(ent);
    case (Entity_kind::EDGE) :
      return getEdgeCentroid<AP>(ent);
    case (Entity_kind::NODE) :
      return getNodeCoordinate<AP>(ent); 
    default :
      assert(false); 
  }
  return AmanziGeometry::Point(); 
}

//=================== 
//    getFace* 
//===================


template<MemSpace_kind MEM>
KOKKOS_INLINE_FUNCTION
decltype(auto) 
MeshCache<MEM>::getFaceEdges(const Entity_ID f) const{
  return RaggedGetter<MEM>::get(
    data_.face_edges_cached, 
    data_.face_edges,
    framework_mesh_,
    [&](const int i) { 
      cEntity_ID_View res; 
      framework_mesh_->getFaceEdges(f,res); 
      return res;
    }, 
    nullptr,
    nullptr, 
    f);
}

template<MemSpace_kind MEM>
KOKKOS_INLINE_FUNCTION
const Entity_ID& 
MeshCache<MEM>::getFaceEdge(const Entity_ID f, const size_type i) const{
  assert(data_.face_edges_cached);
  return data_.face_edges.get<MEM>(f,i);
}

template<MemSpace_kind MEM>
KOKKOS_INLINE_FUNCTION
decltype(auto) // Kokkos::pair<cEntity_ID_View, cEntity_Direction_View>
MeshCache<MEM>::getFaceEdgesAndDirections(const Entity_ID f) const{
  cEntity_ID_View edges; 
  cEntity_Direction_View dirs; 
  framework_mesh_->getFaceEdgesAndDirs(f,edges,&dirs);
  return Kokkos::make_pair(edges, dirs);
}

template<MemSpace_kind MEM>
KOKKOS_INLINE_FUNCTION
void 
MeshCache<MEM>::getFaceEdges(const Entity_ID f, 
                            cEntity_ID_View& edges) const{
  auto [fedges, dirs] = getFaceEdgesAndDirections(f);
  return fedges;
}

template<MemSpace_kind MEM>
KOKKOS_INLINE_FUNCTION
void 
MeshCache<MEM>::getFaceEdgesAndDirs(const Entity_ID f,
                          cEntity_ID_View& edges,
                          cEntity_Direction_View * const dirs) const{
  if constexpr(MEM == MemSpace_kind::DEVICE) {

    if (data_.face_edges_cached) {
      edges = data_.face_edges.getRow<MEM>(f);
      if (dirs) *dirs = data_.face_edge_directions.getRow<MEM>(f);
      return;
    }

  } else {
    if (data_.face_edges_cached) {
      edges = data_.face_edges.getRow<MEM>(f);
      if (dirs) *dirs = data_.face_edge_directions.getRow<MEM>(f);
      return;
    }

    if (framework_mesh_.get()) {
      framework_mesh_->getFaceEdgesAndDirs(f, edges, dirs);
      return;
    }
  }
  assert(false); 
}


template<MemSpace_kind MEM>
template<AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION
size_type
MeshCache<MEM>::getFaceNumCells(const Entity_ID f, const Parallel_kind ptype) const
{
  static_assert(AP != AccessPattern_kind::COMPUTE);
  static_assert(AP != AccessPattern_kind::FRAMEWORK);
  // this is where a generic function would probably help?
  if constexpr(AP == AccessPattern_kind::CACHE) {
    assert(data_.face_cells_cached);
    if (ptype == Parallel_kind::ALL) {
      return data_.face_cells.size<MEM>(f);
    } else {
      int count = 0;
      int n_all = data_.face_cells.size<MEM>(f);
      for (int j=0; j!=n_all; ++j) {
        if (getFaceCell(f,j) < ncells_owned) ++count;
        else break;
      }
      return count;
    }
  } else {
    if (data_.face_cells_cached) return getFaceNumCells<AccessPattern_kind::CACHE>(f);
    return getFaceCells(f).size();
  }
}


template<MemSpace_kind MEM>
KOKKOS_INLINE_FUNCTION
decltype(auto) // cEntity_ID_View
MeshCache<MEM>::getFaceCells(const Entity_ID f, const Parallel_kind ptype) const
{
  cEntity_ID_View fcells; 
  getFaceCells(f, ptype, fcells);
  return fcells;
}


template<MemSpace_kind MEM>
KOKKOS_INLINE_FUNCTION
const Entity_ID&
MeshCache<MEM>::getFaceCell(const Entity_ID f, const size_type i) const
{
  assert(data_.face_cells_cached);
  return data_.face_cells.get<MEM>(f,i);
}


template<MemSpace_kind MEM>
template<AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION
void
MeshCache<MEM>::getFaceCells(const Entity_ID f,
                             const Parallel_kind ptype,
                             cEntity_ID_View & fcells) const
{
  fcells = RaggedGetter<MEM,AP>::get(data_.face_cells_cached,
    data_.face_cells,
    framework_mesh_,
    [&](const int i) { 
      cEntity_ID_View cells; 
      framework_mesh_->getFaceCells(f, ptype, cells); 
      return cells; 
    }, 
    nullptr, 
    nullptr,
    f);
}


// Normal vector of a face
template<MemSpace_kind MEM>
template<AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION
AmanziGeometry::Point MeshCache<MEM>::getFaceNormal(const Entity_ID f) const
{
  return getFaceNormal<AP>(f, -1, nullptr);
}

template<MemSpace_kind MEM>
template<AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION
AmanziGeometry::Point MeshCache<MEM>::getFaceNormal(const Entity_ID f, const Entity_ID c,
        int* orientation) const
{
  AmanziGeometry::Point normal;
  if constexpr (MEM == MemSpace_kind::DEVICE){
    assert(data_.face_geometry_cached); 
  }else {
    if(!data_.face_geometry_cached)
      if (framework_mesh_.get())
        return framework_mesh_->getFaceNormal(f, c, orientation); 
  }

  auto fcells = getFaceCells(f, Parallel_kind::ALL);
  if (orientation) *orientation = 0;

  Entity_ID cc;
  std::size_t i;
  if (c < 0) {
    cc = fcells[0];
    i = 0;
  } else {
    cc = c;
    auto ncells = fcells.size();
    for (i=0; i!=ncells; ++i)
      if (fcells[i] == cc) break;
  }
  normal = data_.face_normals.get<MEM>(f,i);

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
  return normal;

}

template<MemSpace_kind MEM>
template<AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION
decltype(auto) MeshCache<MEM>::getFaceCentroid(const Entity_ID f) const
{
  return Getter<MEM,AP>::get(data_.face_geometry_cached,
    data_.face_centroids,
    framework_mesh_,
    [&](const int i) { return framework_mesh_->getFaceCentroid(i); }, 
    [&](const int i) { return MeshAlgorithms::getFaceCentroid(*this,i);}, 
    f);
}


template<MemSpace_kind MEM>
KOKKOS_INLINE_FUNCTION
decltype(auto) // cEntity_ID_View
MeshCache<MEM>::getFaceNodes(const Entity_ID f) const
{
  cEntity_ID_View fcells; 
  getFaceNodes(f, fcells);
  return fcells;
}


template<MemSpace_kind MEM>
KOKKOS_INLINE_FUNCTION
const Entity_ID&
MeshCache<MEM>::getFaceNode(const Entity_ID f, const size_type i) const
{
  assert(data_.face_nodes_cached);
  return data_.face_nodes.get<MEM>(f,i);
}


template<MemSpace_kind MEM>
KOKKOS_INLINE_FUNCTION
void
MeshCache<MEM>::getFaceNodes(const Entity_ID f,
                             cEntity_ID_View & fcells) const
{
  if constexpr(MEM == MemSpace_kind::DEVICE) {
    if (data_.face_nodes_cached) {
      fcells = data_.face_nodes.getRow<MEM>(f);
      return;
    }

  } else {
    if (data_.face_nodes_cached) {
      fcells = data_.face_nodes.getRow<MEM>(f);
      return;
    }
    if (framework_mesh_.get()) {
      framework_mesh_->getFaceNodes(f, fcells);
      return;
    }
  }
  assert(false); 
}

template<MemSpace_kind MEM>
template<AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION
double MeshCache<MEM>::getFaceArea(const Entity_ID f) const
{
  return Getter<MEM,AP>::get(data_.face_geometry_cached,
    data_.face_areas,
    framework_mesh_,
    [&](const int i) { return framework_mesh_->getFaceArea(i); }, 
    nullptr, 
    f);
}

template<MemSpace_kind MEM> 
template<AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION
decltype(auto) MeshCache<MEM>::getFaceCoordinates(const Entity_ID f) const {
  return RaggedGetter<MEM,AP>::get(data_.face_coordinates_cached,
  data_.face_coordinates,
  framework_mesh_, 
  [&](const int i) { 
    return framework_mesh_->getFaceCoordinates(i);
  }, 
  nullptr, 
  nullptr,
  f);
}

template<MemSpace_kind MEM> 
KOKKOS_INLINE_FUNCTION
std::vector<int> MeshCache<MEM>::getFaceCellEdgeMap(const Entity_ID faceid,
                                 const Entity_ID cellid) const
{
  std::vector<int> map; 
  cEntity_ID_View fedgeids;
  cEntity_Direction_View fedgedirs;

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

template<MemSpace_kind MEM> 
KOKKOS_INLINE_FUNCTION
std::size_t MeshCache<MEM>::getCellMaxFaces() const
{
  std::size_t n(0);
  int ncells = getNumEntities(Entity_kind::CELL, Parallel_kind::OWNED);
  for (int c = 0; c < ncells; ++c) { 
    auto v = getCellNumFaces(c); 
    if(n < v)
      n = v;
  }
  return n;
}

template<MemSpace_kind MEM> 
KOKKOS_INLINE_FUNCTION
std::size_t MeshCache<MEM>::getCellMaxEdges() const
{
  std::size_t n(0);
  if (hasEdges()) {
    int ncells = getNumEntities(Entity_kind::CELL, Parallel_kind::OWNED);
    for (int c = 0; c < ncells; ++c) {
      const auto& edges = getCellEdges(c);
      if(n < edges.size())
        n = edges.size(); 
    }
  }
  return n;
}

// extent
template<MemSpace_kind MEM>
template<AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION
decltype(auto) MeshCache<MEM>::getCellVolume(const Entity_ID c) const
{
  return Getter<MEM,AP>::get(data_.cell_geometry_cached,
    data_.cell_volumes,
    framework_mesh_,
    [&](const int i) { 
      auto v = framework_mesh_->getCellVolume(i); 
      return v; }, 
    nullptr, 
    c);
}

template<MemSpace_kind MEM>
template<AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION
size_type
MeshCache<MEM>::getCellNumFaces(const Entity_ID c) const
{
  static_assert(AP != AccessPattern_kind::COMPUTE);
  static_assert(AP != AccessPattern_kind::FRAMEWORK);
  // this is where a generic function would probably help?
  if constexpr(AP == AccessPattern_kind::CACHE) {
    assert(data_.cell_faces_cached);
    return data_.cell_faces.size<MEM>(c);
  } else {
    if (data_.cell_faces_cached) return getCellNumFaces<AccessPattern_kind::CACHE>(c);
    return getCellFaces(c).size();
  }
}

template<MemSpace_kind MEM>
template<AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION
decltype(auto)
MeshCache<MEM>::getCellFaces(const Entity_ID c) const
{
  cEntity_ID_View cfaces; 
  getCellFaces<AP>(c, cfaces);
  return cfaces; 
}


template<MemSpace_kind MEM>
template<AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION
void
MeshCache<MEM>::getCellFaces(const Entity_ID c,
  cEntity_ID_View& cfaces) const
{
  cfaces = RaggedGetter<MEM,AP>::get(data_.cell_faces_cached,
    data_.cell_faces,
    framework_mesh_, 
    [&](const int i) { 
      cEntity_ID_View cf; 
      framework_mesh_->getCellFaces(i, cf);
      return cf; }, 
    nullptr, 
    nullptr,
    c);
}



template<MemSpace_kind MEM>
KOKKOS_INLINE_FUNCTION
const Entity_ID& MeshCache<MEM>::getCellFace(const Entity_ID c, const size_type i) const
{
  assert(data_.cell_faces_cached);
  return data_.cell_faces.get<MEM>(c,i);
}


template<MemSpace_kind MEM>
KOKKOS_INLINE_FUNCTION
decltype(auto) // Kokkos::pair<cEntity_ID_View, cEntity_Direction_View>
MeshCache<MEM>::getCellFacesAndDirections(const Entity_ID c) const
{
  cEntity_ID_View cfaces; 
  cEntity_Direction_View dirs; 
  getCellFacesAndDirs(c, cfaces, &dirs);
  return Kokkos::pair(cfaces, dirs);
}


template<MemSpace_kind MEM>
KOKKOS_INLINE_FUNCTION
void
MeshCache<MEM>::getCellFacesAndDirs(const Entity_ID c,
                         cEntity_ID_View& faces,
                         cEntity_Direction_View * const dirs) const
{
  if constexpr(MEM == MemSpace_kind::DEVICE) {

    if (data_.cell_faces_cached) {
      faces = data_.cell_faces.getRow<MEM>(c);
      if (dirs) *dirs = data_.cell_face_directions.getRow<MEM>(c);
      return;
    }

  } else {
    if (data_.cell_faces_cached) {
      faces = data_.cell_faces.getRow<MEM>(c);
      if (dirs) *dirs = data_.cell_face_directions.getRow<MEM>(c);
      return;
    }

    if (framework_mesh_.get()) {
      framework_mesh_->getCellFacesAndDirs(c, faces, dirs);
      return;
    }
  }
  assert(false); 
}


template<MemSpace_kind MEM>
KOKKOS_INLINE_FUNCTION
decltype(auto) // Kokkos::pair<cEntity_ID_View, cPoint_View>
MeshCache<MEM>::getCellFacesAndBisectors(const Entity_ID c) const
{
  cEntity_ID_View cfaces;
  cPoint_View bisectors; 
  getCellFacesAndBisectors(c, cfaces, &bisectors);
  return Kokkos::make_pair(cfaces, bisectors);
}


template<MemSpace_kind MEM>
KOKKOS_INLINE_FUNCTION
void
MeshCache<MEM>::getCellFacesAndBisectors(
  const Entity_ID c,
  cEntity_ID_View& faces,
  cPoint_View * const bisectors) const
{
  if constexpr(MEM == MemSpace_kind::DEVICE) {

    if (data_.cell_faces_cached) {
      faces = data_.cell_faces.getRow<MEM>(c);
      if (bisectors) *bisectors = data_.cell_face_bisectors.getRow<MEM>(c);
      return;
    }

  } else {
    if (data_.cell_faces_cached) {
      faces = data_.cell_faces.getRow<MEM>(c);
      if (bisectors) *bisectors = data_.cell_face_bisectors.getRow<MEM>(c);
      return;
    }

    if (framework_mesh_.get()) {
      framework_mesh_->getCellFacesAndBisectors(c, faces, bisectors);
      return;
    }
  }
  assert(false); 
}

template<MemSpace_kind MEM>
template<AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION
decltype(auto) MeshCache<MEM>::getCellCentroid(const Entity_ID c) const
{
  return Getter<MEM,AP>::get(data_.cell_geometry_cached,
    data_.cell_centroids,
    framework_mesh_,
    [&](const int i) { return framework_mesh_->getCellCentroid(i); }, 
    [&](const int i) {
      if constexpr(MEM == MemSpace_kind::DEVICE){
        constexpr auto nnodes = MeshCacheData::static_max_nnodes_; 
        Entity_ID v[nnodes];
        Kokkos::MeshView<Entity_ID*, Kokkos::DefaultExecutionSpace, 
                  Kokkos::MemoryTraits<Kokkos::Unmanaged>> nodes (v, nnodes); 
        getCellNodes(c,nodes);
        AmanziGeometry::Point res;
        for(const auto& n: nodes){
          // !!! res += getnode  doesnt work, it should! 
          res = res + getNodeCoordinate(n); 
        }
        return res/nodes.size(); 
      } else {
        constexpr auto nnodes = MeshCacheData::static_max_nnodes_; 
        cEntity_ID_View nodes;
        getCellNodes(c,nodes);
        AmanziGeometry::Point res;
        for(const auto& n: nodes){
          // !!! res += getnode  doesnt work, it should! 
          res = res + getNodeCoordinate(n); 
        }
        return res/nodes.size(); 
      }
    },
    c);
}

template<MemSpace_kind MEM>
KOKKOS_INLINE_FUNCTION
decltype(auto)
MeshCache<MEM>::getCellNodes(const Entity_ID c) const{
  cEntity_ID_View nodes;
  getCellNodes(c,nodes);
  return nodes;
}

template<MemSpace_kind MEM>
KOKKOS_INLINE_FUNCTION
Entity_ID 
MeshCache<MEM>::getCellNode(const Entity_ID c, const size_type i) const{
  // Compute list and use only one? 
  cEntity_ID_View nodes;
  getCellNodes(c,nodes);  
  return nodes[i];
}

template<MemSpace_kind MEM>
KOKKOS_INLINE_FUNCTION
void
MeshCache<MEM>::getCellNodes(const Entity_ID c, cEntity_ID_View& nodes) const{
  
  nodes = RaggedGetter<MEM,AccessPattern_kind::DEFAULT>::get(data_.cell_nodes_cached,
    data_.cell_nodes,
    framework_mesh_, 
    [&](const int i) { 
      cEntity_ID_View cf; 
      framework_mesh_->getCellNodes(i, cf);
      return cf; }, 
    nullptr, 
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
      auto v = MeshAlgorithms::computeCellNodes(*this,c); 
      nodes = ViewType{v.data(),v.size()}; 
    } else {
      nodes = MeshAlgorithms::computeCellNodes(*this,c); 
    }
  }
  #endif 
}

template<MemSpace_kind MEM> 
template<AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION
decltype(auto) MeshCache<MEM>::getCellCoordinates(const Entity_ID c) const {
  return RaggedGetter<MEM,AP>::get(data_.cell_coordinates_cached,
  data_.cell_coordinates,
  framework_mesh_, 
  [&](const int i) { 
    return framework_mesh_->getCellCoordinates(i);
  }, 
  nullptr, 
  nullptr,
  c);
}


template<MemSpace_kind MEM>
template<AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION
decltype(auto)
MeshCache<MEM>::getCellEdges(const Entity_ID c) const
{
  cEntity_ID_View cedges; 
  getCellEdges<AP>(c, cedges);
  return cedges; 
}


template<MemSpace_kind MEM>
template<AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION
void
MeshCache<MEM>::getCellEdges(const Entity_ID c,
  cEntity_ID_View& cedges) const
{
  cedges = RaggedGetter<MEM,AP>::get(data_.cell_edges_cached,
    data_.cell_edges,
    framework_mesh_, 
    [&](const int i) { 
      cEntity_ID_View ce; 
      framework_mesh_->getCellEdges(i, ce);
      return ce; }, 
    nullptr, 
    nullptr,
    c);
}



template<MemSpace_kind MEM>
KOKKOS_INLINE_FUNCTION
const Entity_ID& MeshCache<MEM>::getCellEdge(const Entity_ID c, const size_type i) const
{
  assert(data_.cell_edges_cached);
  return data_.cell_edges.get<MEM>(c,i);
}

template<MemSpace_kind MEM>
std::size_t MeshCache<MEM>::getCellMaxNodes() const
{
  std::size_t n(0);
  int ncells = getNumEntities(Entity_kind::CELL, Parallel_kind::OWNED);
  for (int c = 0; c < ncells; ++c) {
    auto nodes = getCellNodes(c);
    if(n < nodes.size())
      n = nodes.size(); 
  }
  return n;
}

//=================== 
//    getNode* 
//===================

template<MemSpace_kind MEM> 
template<AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION 
decltype(auto) 
MeshCache<MEM>::getNodeCoordinate(const Entity_ID n) const{
  return Getter<MEM,AP>::get( 
    data_.node_coordinates_cached,
    data_.node_coordinates,
    framework_mesh_,
    [&](const int i) { 
      return framework_mesh_->getNodeCoordinate(i); }, 
    nullptr,  
    n);
}

template<MemSpace_kind MEM> 
template<AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION
decltype(auto) MeshCache<MEM>::getNodeCells(const Entity_ID n,
        const Parallel_kind ptype) const{
  cEntity_ID_View cells; 
  getNodeCells<AP>(n,ptype,cells); 
  return cells; 
}

template<MemSpace_kind MEM> 
template<AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION
void MeshCache<MEM>::getNodeCells(const Entity_ID n,
                  const Parallel_kind ptype,
                  cEntity_ID_View& cells) const{
  cells = RaggedGetter<MEM,AP>::get(data_.node_cells_cached,
    data_.node_cells,
    framework_mesh_, 
    [&](const int i) { 
      cEntity_ID_View cells; 
      framework_mesh_->getNodeCells(i,ptype,cells);
      return cells; 
    }, 
    nullptr,
    nullptr,
    n);
}

template<MemSpace_kind MEM>
template<AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION
decltype(auto) MeshCache<MEM>::getNodeFaces(const Entity_ID n, const Parallel_kind ptype) const
{
  cEntity_ID_View faces; 
  getNodeFaces<AP>(n,ptype,faces); 
  return faces; 
}

template<MemSpace_kind MEM>
template<AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION
void MeshCache<MEM>::getNodeFaces(const Entity_ID n,
                  const Parallel_kind ptype,
                  cEntity_ID_View& faces) const
{
  faces = RaggedGetter<MEM,AP>::get(data_.node_faces_cached,
    data_.node_faces,
    framework_mesh_, 
    [&](const int i) { 
      cEntity_ID_View faces; 
      framework_mesh_->getNodeFaces(i,ptype,faces);
      return faces; 
    }, 
    nullptr, 
    nullptr,
    n);
}

template<MemSpace_kind MEM>
KOKKOS_INLINE_FUNCTION
Cell_kind MeshCache<MEM>::getCellType(const Entity_ID c) const
{
  return MeshAlgorithms::getCellType(*this, c);
}

template<MemSpace_kind MEM>
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

template<MemSpace_kind MEM>
template<AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION
decltype(auto) MeshCache<MEM>::getEdgeNodes(const Entity_ID e) const {
  return RaggedGetter<MEM,AP>::get(data_.edge_nodes_cached,
    data_.edge_nodes,
    framework_mesh_, 
    [&](const int i) { 
      cEntity_ID_View nodes; 
      framework_mesh_->getEdgeNodes(i,nodes);
      return nodes; 
    }, 
    nullptr, 
    nullptr,
    e);
}

template<MemSpace_kind MEM>
template<AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION
decltype(auto) MeshCache<MEM>::getEdgeNode(const Entity_ID e, const size_type i) const{
  // Compute list and use only one? 
  cEntity_ID_View nodes;
  getEdgeNodes(e,nodes);  
  return nodes[i];
}


  //[[deprecated("Prefer to use non-void variant that returns nodes directly")]]
template<MemSpace_kind MEM>
template<AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION
void MeshCache<MEM>::getEdgeNodes(const Entity_ID e, cEntity_ID_View& nodes) const{
  nodes = RaggedGetter<MEM,AP>::get(data_.edge_nodes_cached,
    data_.edge_nodes,
    framework_mesh_, 
    [&](const int i) { 
      cEntity_ID_View nodes; 
      framework_mesh_->getEdgeNodes(i,nodes);
      return nodes; 
    }, 
    nullptr, 
    nullptr,
    e);
}

template<MemSpace_kind MEM>
template<AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION
decltype(auto) MeshCache<MEM>::getEdgeCoordinates(const Entity_ID c) const {
  return RaggedGetter<MEM,AP>::get(data_.edge_coordinates_cached,
  data_.edge_coordinates,
  framework_mesh_, 
  [&](const int i) { 
    return framework_mesh_->getEdgeCoordinates(i);
  }, 
  nullptr, 
  nullptr,
  c);
}

template<MemSpace_kind MEM>
template<AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION
decltype(auto) MeshCache<MEM>::getEdgeCentroid(const Entity_ID c) const {
  return Getter<MEM,AP>::get(data_.edge_geometry_cached,
  data_.edge_centroids,
  framework_mesh_, 
  [&](const int i) { 
    return framework_mesh_->getEdgeCentroid(i);
  }, 
  nullptr,
  c);
}

template<MemSpace_kind MEM>
template<AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION
decltype(auto) MeshCache<MEM>::getEdgeVector(const Entity_ID e) const {
  return Getter<MEM,AP>::get(data_.edge_geometry_cached,
    data_.edge_vectors,
    framework_mesh_, 
    [&](const int i) { 
      return framework_mesh_->getEdgeVector(i);
    }, 
    nullptr,
    e);
}


template<MemSpace_kind MEM>
template<AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION
double MeshCache<MEM>::getEdgeLength(const Entity_ID e) const {
  return Getter<MEM,AP>::get(data_.edge_lengths_cached,
    data_.edge_lengths,
    framework_mesh_, 
    [&](const int i) { 
      return framework_mesh_->getEdgeLength(i);
    }, 
    nullptr,
    e);
}

template<MemSpace_kind MEM>
template<AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION
decltype(auto) MeshCache<MEM>::getEdgeCells(const Entity_ID e, const Parallel_kind ptype) const {
  cEntity_ID_View cells; 
  getEdgeCells<AP>(e,ptype,cells); 
  return cells; 
}

template<MemSpace_kind MEM>
template<AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION
void MeshCache<MEM>::getEdgeCells(const Entity_ID e,
                  const Parallel_kind ptype,
                  cEntity_ID_View& cells) const {
  cells = RaggedGetter<MEM,AP>::get(data_.edge_cells_cached,
  data_.edge_cells,
  framework_mesh_, 
  [&](const int i) { 
    cEntity_ID_View cells; 
    framework_mesh_->getEdgeCells(i,ptype,cells);
    return cells; 
  }, 
  nullptr, 
  nullptr,
  e);
}

template<MemSpace_kind MEM>
template<AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION
decltype(auto) MeshCache<MEM>::getEdgeFaces(const Entity_ID e, const Parallel_kind ptype) const
{
  cEntity_ID_View faces; 
  getEdgeFaces<AP>(e,ptype,faces); 
  return faces; 
}

template<MemSpace_kind MEM>
template<AccessPattern_kind AP>
KOKKOS_INLINE_FUNCTION
void MeshCache<MEM>::getEdgeFaces(const Entity_ID e,
                  const Parallel_kind ptype,
                  cEntity_ID_View& faces) const
{
  faces = RaggedGetter<MEM,AP>::get(data_.edge_faces_cached,
  data_.edge_faces,
  framework_mesh_, 
  [&](const int i) { 
    cEntity_ID_View faces; 
    framework_mesh_->getEdgeFaces(i,ptype,faces);
    return faces; 
  }, 
  nullptr, 
  nullptr,
  e);
}

//
// Build the cache, fine grained control
// =============================================
template<MemSpace_kind MEM>
void MeshCache<MEM>::cacheCellGeometry()
{
  assert(framework_mesh_.get());
  if (data_.cell_geometry_cached) return;
  data_.cell_volumes.resize(ncells_all);
  data_.cell_centroids.resize(ncells_all);
  for (Entity_ID i=0; i!=ncells_all; ++i) {
    // note this must be on host
    std::tie(view<MemSpace_kind::HOST>(data_.cell_volumes)[i],
             view<MemSpace_kind::HOST>(data_.cell_centroids)[i]) =
      framework_mesh_->computeCellGeometry(i);
  }
  Kokkos::deep_copy(data_.cell_volumes.view_device(), data_.cell_volumes.view_host()); 
  Kokkos::deep_copy(data_.cell_centroids.view_device(), data_.cell_centroids.view_host()); 
  data_.cell_geometry_cached = true;
}

// cell-face adjacencies
template<MemSpace_kind MEM>
void MeshCache<MEM>::cacheCellFaces()
{
  assert(framework_mesh_.get());
  if (data_.cell_faces_cached) return;
  int num_cells = framework_mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::ALL);  
  auto lambda1 = [this](Entity_ID c, cEntity_ID_View& cfaces) {this->framework_mesh_->getCellFaces(c, cfaces);};
  data_.cell_faces = asRaggedArray_DualView<Entity_ID>(lambda1, num_cells);
  auto lambda2 = [this](Entity_ID c, cEntity_Direction_View& dirs) { this->framework_mesh_->getCellFaceDirs(c, dirs); };
  data_.cell_face_directions = asRaggedArray_DualView<int>(lambda2, num_cells);
  data_.cell_faces_cached = true;
}

template<MemSpace_kind MEM>
void MeshCache<MEM>::cacheCellCoordinates()
{
  assert(framework_mesh_.get());
  if (data_.cell_coordinates_cached) return;
  int num_cells = framework_mesh_->getNumEntities(
    AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::ALL); 
  auto lambda = [this](Entity_ID c, cPoint_View& ccoords) 
    { ccoords = this->framework_mesh_->getCellCoordinates(c); };
  data_.cell_coordinates = 
    asRaggedArray_DualView<Amanzi::AmanziGeometry::Point>(lambda, num_cells);
  data_.cell_coordinates_cached = true;
}

template<MemSpace_kind MEM>
void MeshCache<MEM>::cacheCellEdges()
{
  assert(framework_mesh_.get());
  if (data_.cell_edges_cached) return;
  framework_mesh_->hasEdgesOrThrow();
  int num_cells = framework_mesh_->getNumEntities(
    AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::ALL); 
  auto lambda = [this](Entity_ID c, cEntity_ID_View& cedges) 
    { this->framework_mesh_->getCellEdges(c,cedges); };
  data_.cell_edges = 
    asRaggedArray_DualView<Entity_ID>(lambda, num_cells);
  data_.cell_edges_cached = true;
}

template<MemSpace_kind MEM>
void MeshCache<MEM>::cacheCellNodes()
{
  assert(framework_mesh_.get());
  if (data_.cell_nodes_cached) return;
  int num_cells = framework_mesh_->getNumEntities(
    AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::ALL); 
  auto lambda = [this](Entity_ID c, cEntity_ID_View& cnodes) 
    { this->framework_mesh_->getCellNodes(c,cnodes); };
  data_.cell_nodes = 
    asRaggedArray_DualView<Entity_ID>(lambda, num_cells);
  data_.cell_nodes_cached = true;
}

template<MemSpace_kind MEM>
void MeshCache<MEM>::cacheFaceGeometry()
{
  assert(framework_mesh_.get());
  if (data_.face_geometry_cached) return;
  data_.face_areas.resize(nfaces_all);
  data_.face_centroids.resize(nfaces_all);

  // slurp down the RaggedArray for normals using a lambda that, as a side
  // effect, captures area and centroid too.
  auto area_view = view<MemSpace_kind::HOST>(data_.face_areas);
  auto centroid_view = view<MemSpace_kind::HOST>(data_.face_centroids);
  auto lambda = [&,this](const Entity_ID& f, cPoint_View& normals) {
    auto area_cent_normal = this->framework_mesh_->computeFaceGeometry(f);
    area_view[f] = std::get<0>(area_cent_normal);
    centroid_view[f] = std::get<1>(area_cent_normal);
    normals = std::get<2>(area_cent_normal);
  };
  data_.face_normals = asRaggedArray_DualView<AmanziGeometry::Point>(lambda, nfaces_all);

  // still must sync areas/centroids
  Kokkos::deep_copy(view<MemSpace_kind::DEVICE>(data_.face_areas),view<MemSpace_kind::HOST>(data_.face_areas)); 
  Kokkos::deep_copy(view<MemSpace_kind::DEVICE>(data_.face_centroids),view<MemSpace_kind::HOST>(data_.face_centroids)); 
  data_.face_geometry_cached = true;

  // cache normal directions -- make this a separate call?  Think about
  // granularity here.
  auto lambda2 = [&,this](const Entity_ID& f, cEntity_Direction_View& dirs) {
    Entity_Direction_View ldirs; 
    // This NEEDS to call the framework or be passed an host mesh to call the function on the host. 
    cEntity_ID_View fcells; 
    framework_mesh_->getFaceCells(f, Parallel_kind::ALL, fcells);
    Kokkos::resize(ldirs,fcells.size());
    for (int i=0; i!=fcells.size(); ++i) {
      this->framework_mesh_->getFaceNormal(f, fcells[i], &ldirs[i]);
    }
    dirs = ldirs; 
  };
  data_.face_normal_directions = asRaggedArray_DualView<int>(lambda2, nfaces_all);

  // cache cell-face-bisectors -- make this a separate call?  Think about
  // granularity here.
  auto lambda3 = [&,this](const Entity_ID& c, cPoint_View& bisectors) {
    cEntity_ID_View cfaces;
    this->framework_mesh_->getCellFacesAndBisectors(c, cfaces, &bisectors);
  };
  data_.cell_face_bisectors = asRaggedArray_DualView<AmanziGeometry::Point>(lambda3, ncells_all);

}


// face-cell adjacencies
template<MemSpace_kind MEM>
void MeshCache<MEM>::cacheFaceCells()
{
  if (data_.face_cells_cached) return;
  auto lambda = [this](Entity_ID f, cEntity_ID_View& fcells) { this->framework_mesh_->getFaceCells(f, Parallel_kind::ALL, fcells); };
  data_.face_cells = asRaggedArray_DualView<Entity_ID>( lambda, nfaces_all);
  data_.face_cells_cached = true;
}

// face-cell adjacencies
template<MemSpace_kind MEM>
void MeshCache<MEM>::cacheFaceEdges()
{
  if (data_.face_edges_cached) return;
  framework_mesh_->hasEdgesOrThrow();
  auto lambda = [this](Entity_ID f, cEntity_ID_View& fedges) 
    { this->framework_mesh_->getFaceEdges(f, fedges);};
  auto lambda_dir = [this](Entity_ID f, cEntity_Direction_View& dirs)
    { cEntity_ID_View l;  this->framework_mesh_->getFaceEdgesAndDirs(f,l,&dirs);};
  data_.face_edges = asRaggedArray_DualView<Entity_ID>( lambda, nfaces_all);
  data_.face_edge_directions = asRaggedArray_DualView<int>( lambda_dir, nfaces_all);
  data_.face_edges_cached = true;
}

// face-node adjacencies
template<MemSpace_kind MEM>
void MeshCache<MEM>::cacheFaceNodes()
{
  if (data_.face_nodes_cached) return;
  auto lambda = [this](Entity_ID f, cEntity_ID_View& fnodes) { this->framework_mesh_->getFaceNodes(f, fnodes); };
  data_.face_nodes = asRaggedArray_DualView<Entity_ID>(lambda, nfaces_all);
  data_.face_nodes_cached = true;
}

template<MemSpace_kind MEM>
void MeshCache<MEM>::cacheFaceCoordinates()
{
  assert(framework_mesh_.get());
  if (data_.face_coordinates_cached) return;
  auto lambda = [this](Entity_ID f, cPoint_View& fcoords) 
    { fcoords = this->framework_mesh_->getFaceCoordinates(f); };
  data_.face_coordinates = asRaggedArray_DualView<Amanzi::AmanziGeometry::Point>(lambda, nfaces_all);

  data_.face_coordinates_cached = true;
}

// edge centroid, length, vector
template<MemSpace_kind MEM>
void MeshCache<MEM>::cacheEdgeGeometry()
{
  assert(framework_mesh_.get());
  if (data_.edge_geometry_cached) return;
  framework_mesh_->hasEdgesOrThrow();
  data_.edge_vectors.resize(nedges_all);
  data_.edge_centroids.resize(nedges_all);

  auto evector = view<MemSpace_kind::HOST>(data_.edge_vectors);
  auto ecentroids = view<MemSpace_kind::HOST>(data_.edge_centroids);
  for(int i = 0 ; i < nedges_all; ++i){
    auto egeometry = this->framework_mesh_->computeEdgeGeometry(i);
    evector[i] = std::get<0>(egeometry);
    ecentroids[i] = std::get<1>(egeometry);
  }
  Kokkos::deep_copy(view<MemSpace_kind::DEVICE>(data_.edge_centroids),view<MemSpace_kind::HOST>(data_.edge_centroids)); 
  Kokkos::deep_copy(view<MemSpace_kind::DEVICE>(data_.edge_vectors),view<MemSpace_kind::HOST>(data_.edge_vectors)); 
  data_.edge_geometry_cached = true;
}

// edge-cell adjacencies
template<MemSpace_kind MEM>
void MeshCache<MEM>::cacheEdgeCells()
{
  assert(framework_mesh_.get());
  if (data_.edge_cells_cached) return;
  framework_mesh_->hasEdgesOrThrow();
  auto lambda = [this](Entity_ID f, cEntity_ID_View& ecells) 
    { this->framework_mesh_->getEdgeCells(f, Parallel_kind::ALL, ecells); };
  data_.edge_cells = asRaggedArray_DualView<Entity_ID>(lambda, nedges_all);
  data_.edge_cells_cached = true;
}

// edge-face adjacencies
template<MemSpace_kind MEM>
void MeshCache<MEM>::cacheEdgeFaces()
{
  assert(framework_mesh_.get());
  if (data_.edge_faces_cached) return;
  framework_mesh_->hasEdgesOrThrow();
  auto lambda = [this](Entity_ID f, cEntity_ID_View& efaces) 
    { this->framework_mesh_->getEdgeFaces(f, Parallel_kind::ALL, efaces); };
  data_.edge_faces = asRaggedArray_DualView<Entity_ID>(lambda, nedges_all);
  data_.edge_faces_cached = true;
}

// edge-node adjacencies
template<MemSpace_kind MEM>
void MeshCache<MEM>::cacheEdgeNodes()
{
  assert(framework_mesh_.get());
  if (data_.edge_nodes_cached) return;
  framework_mesh_->hasEdgesOrThrow();
  auto lambda = [this](Entity_ID f, cEntity_ID_View& efaces) 
    { this->framework_mesh_->getEdgeNodes(f, efaces); };
  data_.edge_nodes = asRaggedArray_DualView<Entity_ID>(lambda, nedges_all);
  data_.edge_nodes_cached = true;
}

// edge coordinates
template<MemSpace_kind MEM>
void MeshCache<MEM>::cacheEdgeCoordinates()
{
  assert(framework_mesh_.get());
  if (data_.edge_coordinates_cached) return;
  framework_mesh_->hasEdgesOrThrow();
  auto lambda = [this](Entity_ID f, cPoint_View& enodes) 
    { enodes = this->framework_mesh_->getEdgeCoordinates(f); };
  data_.edge_coordinates = asRaggedArray_DualView<Amanzi::AmanziGeometry::Point>(lambda, nedges_all);
  data_.edge_coordinates_cached = true;
}

// // node-cell adjacencies
template<MemSpace_kind MEM>
void MeshCache<MEM>::cacheNodeCells()
{
  assert(framework_mesh_.get());
  if (data_.node_cells_cached) return;
  auto lambda = [this](Entity_ID f, cEntity_ID_View& ncells) 
    { this->framework_mesh_->getNodeCells(f,Parallel_kind::ALL, ncells); };
  data_.node_cells = 
    asRaggedArray_DualView<Entity_ID>(lambda, nnodes_all);
  data_.node_cells_cached = true;
}

// // node-face adjacencies
template<MemSpace_kind MEM>
void MeshCache<MEM>::cacheNodeFaces()
{
  assert(framework_mesh_.get());
  if (data_.node_faces_cached) return;
  auto lambda = [this](Entity_ID f, cEntity_ID_View& nfaces) 
    { this->framework_mesh_->getNodeFaces(f,Parallel_kind::ALL, nfaces); };
  data_.node_faces = 
    asRaggedArray_DualView<Entity_ID>(lambda, nnodes_all);
  data_.node_faces_cached = true;
}

// // node-edge adjacencies
template<MemSpace_kind MEM>
void MeshCache<MEM>::cacheNodeEdges()
{
  assert(framework_mesh_.get());
  if (data_.node_edges_cached) return;
  framework_mesh_->hasEdgesOrThrow();
  auto lambda = [this](Entity_ID f, cEntity_ID_View& nedges) 
    { this->framework_mesh_->getNodeEdges(f,Parallel_kind::ALL, nedges); };
  data_.node_edges = 
    asRaggedArray_DualView<Entity_ID>(lambda, nnodes_all);
  data_.node_edges_cached = true;
}


// node coordinates
template<MemSpace_kind MEM>
void MeshCache<MEM>::cacheNodeCoordinates()
{
  if (data_.node_coordinates_cached) return;
  data_.node_coordinates.resize(nnodes_all);
  for (Entity_ID i=0; i!=nnodes_all; ++i) {
    view<MemSpace_kind::HOST>(data_.node_coordinates)[i] = framework_mesh_->getNodeCoordinate(i);
  }

  Kokkos::deep_copy(data_.node_coordinates.d_view, data_.node_coordinates.h_view);
  data_.node_coordinates_cached = true;
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

template<MemSpace_kind MEM> 
void MeshCache<MEM>::setNodeCoordinate(const Entity_ID n, const AmanziGeometry::Point& coord)
{
  auto bf = view<MemSpace_kind::HOST>(data_.node_coordinates)[n]; 
  if (framework_mesh_.get())
    framework_mesh_->setNodeCoordinate(n, coord);
  if (data_.node_coordinates_cached)
    view<MemSpace_kind::HOST>(data_.node_coordinates)[n] = coord;
}

// // common error messaging
// void MeshCache<MEM>throwAccessError_(const std::string& func_name) const
// {
//   Errors::Message msg;
//   msg << "MeshCache<MEM>" << func_name << " cannot compute this quantity -- not cached and framework does not exist.";
//   Exceptions::amanzi_throw(msg);
// }

template<MemSpace_kind MEM> 
Entity_ID
MeshCache<MEM>::getEntityParent(const Entity_kind kind, const Entity_ID entid) const
{
  if (data_.parent_entities_cached) {
    switch(kind) {
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
      default: {}
    }
  } else if (framework_mesh_.get()) {
    return framework_mesh_->getEntityParent(kind, entid);
  }
  return -1;
}

template<MemSpace_kind MEM> 
void
MeshCache<MEM>::PrintMeshStatistics() const
{
  assert(false && "Need verbose object"); 
  #if 0 
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
  #endif 
}

template<MemSpace_kind MEM> 
void MeshCache<MEM>::recacheGeometry()
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
  if(data_.face_coordinates_cached) {
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


namespace MeshAlgorithms {

template<MemSpace_kind MEM>
void cacheAll(MeshCache<MEM>& mesh)
{
  std::cout<<"############# CacheAll"<<std::endl;
  // caches everything, likely just for testing
  mesh.cacheNodeCoordinates();

  mesh.cacheCellFaces();
  mesh.cacheCellNodes();
  mesh.cacheCellCoordinates();
  mesh.cacheFaceCells();
  mesh.cacheFaceNodes();
  mesh.cacheFaceCoordinates();
  mesh.cacheNodeCells();
  mesh.cacheNodeFaces();
  mesh.cacheCellGeometry();
  mesh.cacheFaceGeometry();
  if (mesh.hasEdges()) {
    mesh.cacheCellEdges();
    mesh.cacheEdgeCells();
    mesh.cacheFaceEdges();
    mesh.cacheEdgeFaces();
    mesh.cacheNodeEdges();
    mesh.cacheEdgeNodes();
    mesh.cacheEdgeGeometry();
    mesh.cacheEdgeCoordinates();
   }
}

template<MemSpace_kind MEM>
void cacheDefault(MeshCache<MEM>& mesh)
{
  std::cout<<"############# CacheDefault"<<std::endl;
  // caches what the developers currently think is best
  mesh.cacheNodeCoordinates();
  mesh.cacheCellFaces();
  mesh.cacheFaceCells();
  if (mesh.hasNodes()) {
    mesh.cacheFaceNodes();
  }
  mesh.cacheCellGeometry();
  mesh.cacheFaceGeometry();
  if (mesh.hasEdges()) {
    mesh.cacheFaceEdges();
    mesh.cacheEdgeFaces();
    mesh.cacheEdgeGeometry();
  }
}

} // namespace MeshAlgorithm
} // namespace AmanziMesh
} // namespace Amanzi
