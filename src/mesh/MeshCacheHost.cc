#include "Teuchos_CommHelpers.hpp"

#include "Geometry.hh"
#include "Iterators.hh"
#include "MeshUtils.hh"
#include "MeshFramework.hh"
#include "MeshCacheHost_impl.hh"
#include "MeshAlgorithms.hh"
#include "MeshSets.hh"

namespace Amanzi::AmanziMesh {

MeshCacheHost::MeshCacheHost(const MeshCacheDevice& other): MeshCacheBase(other) {
  parent_ = Teuchos::RCP<const MeshCacheHost>(other.getParentMesh() == Teuchos::null ? Teuchos::null :
          onMemHost(other.getParentMesh()));
}

void
MeshCacheHost::setParentMesh(const Teuchos::RCP<const MeshCacheHost>& parent)
{
  if (parent_ != Teuchos::null && parent_ != parent) {
    Errors::Message msg("MeshCacheHost::setParentMesh given conflicting parent mesh.");
    Exceptions::amanzi_throw(msg);
  } else {
    parent_ = parent;
    cacheParentEntities();
  }
}


void
MeshCacheHost::setMeshFramework(const Teuchos::RCP<MeshFramework>& framework_mesh)
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



bool
MeshCacheHost::isValidSetType(const AmanziGeometry::RegionType rtype, const Entity_kind kind) const
{
  if (rtype == AmanziGeometry::RegionType::LABELEDSET && framework_mesh_.get()) {
    return framework_mesh_->isValidSetType(rtype, kind);
  }
  return true;
}


int
MeshCacheHost::getSetSize(const std::string& region_name,
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
typename MeshCacheHost::cEntity_ID_View
MeshCacheHost::getSetEntities(const std::string& region_name,
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

    sets_[key_all] = asDualView(resolveMeshSet(*region, kind, Parallel_kind::ALL, *this));
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

Kokkos::pair<typename MeshCacheHost::cEntity_ID_View, typename MeshCacheHost::cDouble_View>
MeshCacheHost::getSetEntitiesAndVolumeFractions(const std::string& region_name,
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
    sets_[key] =
      asDualView(resolveMeshSetVolumeFractions(*region, kind, ptype, vol_fracs_list, *this));
    set_vol_fracs_[key] = asDualView<double>(vol_fracs_list);
  }
  return Kokkos::pair(view<MEM>(sets_.at(key)), view<MEM>(set_vol_fracs_.at(key)));
}


Entity_ID
MeshCacheHost::getNumEntities(const Entity_kind kind, const Parallel_kind ptype) const
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

//===================
//    getFace*
//===================

typename MeshCacheHost::cEntity_ID_View
MeshCacheHost::getFaceEdges(const Entity_ID f) const
{
  return Impl::RaggedGetter<MEM>::get(
    data_.face_edges_cached,
    data_.face_edges,
    framework_mesh_,
    [&](const int i) {
      MeshFramework::cEntity_ID_View res;
      framework_mesh_->getFaceEdges(f, res);
      return res;
    },
    nullptr,
    f);
}

const Entity_ID&
MeshCacheHost::getFaceEdge(const Entity_ID f, const size_type i) const
{
  assert(data_.face_edges_cached);
  return data_.face_edges.get<MEM>(f, i);
}

  Kokkos::pair<typename MeshCacheHost::cEntity_ID_View, typename MeshCacheHost::cDirection_View>
  MeshCacheHost::getFaceEdgesAndDirections(const Entity_ID f) const
{
  MeshFramework::cEntity_ID_View edges;
  MeshFramework::cDirection_View dirs;
  framework_mesh_->getFaceEdgesAndDirs(f, edges, &dirs);
  return Kokkos::make_pair(edges, dirs);
}

void MeshCacheHost::getFaceEdges(const Entity_ID f, cEntity_ID_View& edges) const
{
  auto [fedges, dirs] = getFaceEdgesAndDirections(f);
  edges = fedges;
}

void 
MeshCacheHost::getFaceEdgesAndDirs(const Entity_ID f,
                                    cEntity_ID_View& edges,
                                    cDirection_View* const dirs) const
{
    if (data_.face_edges_cached) {
      edges = data_.face_edges.getRowUnmanaged<MEM>(f);
      if (dirs) *dirs = data_.face_edge_directions.getRowUnmanaged<MEM>(f);
      return;
    }

    if (framework_mesh_.get()) {
      framework_mesh_->getFaceEdgesAndDirs(f, edges, dirs);
      return;
    }
    assert(false);
}

const Entity_ID&
MeshCacheHost::getFaceCell(const Entity_ID f, const size_type i) const
{
  assert(data_.face_cells_cached);
  return data_.face_cells.get<MEM>(f, i);
}

typename MeshCacheHost::cEntity_ID_View
MeshCacheHost::getFaceNodes(const Entity_ID f) const
{
  cEntity_ID_View fcells;
  getFaceNodes(f, fcells);
  return fcells;
}

const Entity_ID&
MeshCacheHost::getFaceNode(const Entity_ID f, const size_type i) const
{
  assert(data_.face_nodes_cached);
  return data_.face_nodes.get<MEM>(f, i);
}

void
MeshCacheHost::getFaceNodes(const Entity_ID f, cEntity_ID_View& fcells) const
{
    if (data_.face_nodes_cached) {
      fcells = data_.face_nodes.getRowUnmanaged<MEM>(f);
      return;
    }
    if (framework_mesh_.get()) {
      framework_mesh_->getFaceNodes(f, fcells);
      return;
    }
  assert(false);
}

std::vector<int>
MeshCacheHost::getFaceCellEdgeMap(const Entity_ID faceid, const Entity_ID cellid) const
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

std::size_t
MeshCacheHost::getCellMaxFaces() const
{
  std::size_t n(0);
  int ncells = getNumEntities(Entity_kind::CELL, Parallel_kind::OWNED);
  for (int c = 0; c < ncells; ++c) {
    auto v = getCellNumFaces(c);
    if (n < v) n = v;
  }
  return n;
}

std::size_t
MeshCacheHost::getCellMaxEdges() const
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

const Entity_ID&
MeshCacheHost::getCellFace(const Entity_ID c, const size_type i) const
{
  assert(data_.cell_faces_cached);
  return data_.cell_faces.get<MEM>(c, i);
}

  Kokkos::pair<typename MeshCacheHost::cEntity_ID_View, typename MeshCacheHost::cDirection_View>
  MeshCacheHost::getCellFacesAndDirections(const Entity_ID c) const
{
  cEntity_ID_View cfaces;
  cDirection_View dirs;
  getCellFacesAndDirs(c, cfaces, &dirs);
  return Kokkos::pair(cfaces, dirs);
}

void
MeshCacheHost::getCellFacesAndDirs(const Entity_ID c,
                                    cEntity_ID_View& faces,
                                    cDirection_View* const dirs) const
{
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


  Kokkos::pair<typename MeshCacheHost::cEntity_ID_View, typename MeshCacheHost::cPoint_View>
  MeshCacheHost::getCellFacesAndBisectors(const Entity_ID c) const
{
  cEntity_ID_View cfaces;
  cPoint_View bisectors;
  getCellFacesAndBisectors(c, cfaces, &bisectors);
  return Kokkos::make_pair(cfaces, bisectors);
}

void
MeshCacheHost::getCellFacesAndBisectors(const Entity_ID c,
                                         cEntity_ID_View& faces,
                                         cPoint_View* const bisectors) const
{
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

typename MeshCacheHost::cEntity_ID_View
MeshCacheHost::getCellNodes(const Entity_ID c) const
{
  cEntity_ID_View nodes;
  getCellNodes(c, nodes);
  return nodes;
}

Entity_ID
MeshCacheHost::getCellNode(const Entity_ID c, const size_type i) const
{
  // Compute list and use only one?
  cEntity_ID_View nodes;
  getCellNodes(c, nodes);
  return nodes[i];
}

void
MeshCacheHost::getCellNodes(const Entity_ID c, cEntity_ID_View& nodes) const
{
  nodes = Impl::RaggedGetter<MEM, AccessPattern_kind::DEFAULT>::get(
    data_.cell_nodes_cached,
    data_.cell_nodes,
    framework_mesh_,
    [&](const int i) {
      MeshFramework::cEntity_ID_View cf;
      framework_mesh_->getCellNodes(i, cf);
      return cf;
    },
    nullptr,
    c);
}

const Entity_ID&
MeshCacheHost::getCellEdge(const Entity_ID c, const size_type i) const
{
  assert(data_.cell_edges_cached);
  return data_.cell_edges.get<MEM>(c, i);
}

std::size_t
MeshCacheHost::getCellMaxNodes() const
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

void
MeshCacheHost::setNodeCoordinate(const Entity_ID n, const AmanziGeometry::Point& p)
{
  if (data_.node_coordinates_cached) view<MemSpace_kind::HOST>(data_.node_coordinates)(n) = p;
  if (framework_mesh_.get()) framework_mesh_->setNodeCoordinate(n, p);
}

void
MeshCacheHost::setNodeCoordinates(const cEntity_ID_View& nodes, const cPoint_View& new_coords)
{
  auto bf = view<MEM>(data_.node_coordinates);
  if (data_.node_coordinates_cached) {
    for (int i = 0; i != nodes.size(); ++i) { bf(nodes(i)) = new_coords(i); }
  }

  if (framework_mesh_.get()) {
    MeshFramework::cEntity_ID_View nodes_on_host;
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

Cell_kind
MeshCacheHost::getCellType(const Entity_ID c) const
{
  return Impl::getCellType(*this, c);
}

Parallel_kind
MeshCacheHost::getParallelType(const Entity_kind& kind, const Entity_ID id) const
{
  if (id < getNumEntities(kind, Parallel_kind::OWNED)) {
    return Parallel_kind::OWNED;
  } else if (id < getNumEntities(kind, Parallel_kind::ALL)) {
    return Parallel_kind::GHOST;
  }
  return Parallel_kind::UNKNOWN;
}


Entity_ID
MeshCacheHost::getEntityParent(const Entity_kind kind, const Entity_ID entid) const
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
  } else if (framework_mesh_) {
    return framework_mesh_->getEntityParent(kind, entid);
  }
  return -1;
}

typename MeshCacheHost::cEntity_ID_View
MeshCacheHost::getEntityParents(const Entity_kind kind) const
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
  return MeshCacheHost::cEntity_ID_View();
}

bool
MeshCacheHost::isPointInCell(const AmanziGeometry::Point& p, const Entity_ID cellid) const
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


void
MeshCacheHost::PrintMeshStatistics() const
{
  // auto vo_ = Teuchos::rcp(new VerboseObject("Mesh Output", *plist_));
  // if (vo_.get() && vo_->getVerbLevel() >= Teuchos::VERB_LOW) {
  //   int ncells = getNumEntities(AmanziMesh::CELL, AmanziMesh::Parallel_kind::OWNED);
  //   int nfaces = getNumEntities(AmanziMesh::FACE, AmanziMesh::Parallel_kind::OWNED);
  //   int nnodes = getNumEntities(AmanziMesh::NODE, AmanziMesh::Parallel_kind::OWNED);
  //   int nedges(0);
  //   if (has_edges_) nedges = getNumEntities(AmanziMesh::EDGE, AmanziMesh::Parallel_kind::OWNED);

  //   int min_out[4], max_out[4], sum_out[4], tmp_in[4] = { ncells, nfaces, nedges, nnodes };
  //   getComm()->MinAll(tmp_in, min_out, 4);
  //   getComm()->MaxAll(tmp_in, max_out, 4);
  //   getComm()->SumAll(tmp_in, sum_out, 4);

  //   Teuchos::OSTab tab = vo_->getOSTab();
  //   *vo_->os() << "cells, tot/min/max: " << sum_out[0] << "/" << min_out[0] << "/" << max_out[0]
  //              << "\n";
  //   *vo_->os() << "faces, tot/min/max: " << sum_out[1] << "/" << min_out[1] << "/" << max_out[1]
  //              << "\n";
  //   if (has_edges_)
  //     *vo_->os() << "edges, tot/min/max: " << sum_out[2] << "/" << min_out[2] << "/" << max_out[2]
  //                << "\n";
  //   *vo_->os() << "nodes, tot/min/max: " << sum_out[3] << "/" << min_out[3] << "/" << max_out[3]
  //              << "\n\n";
  // }
}


  
//
// Build the cache, fine grained control
// =============================================

void
MeshCacheHost::cacheCellGeometry()
{
  static_assert(MEM == MemSpace_kind::HOST);
  assert(framework_mesh_.get());
  if (data_.cell_geometry_cached) return;
  auto lambda = [this](Entity_ID c, Double_View cvol, Point_View ccent) {
    std::tie(cvol[c], ccent[c]) = algorithms_->computeCellGeometry(*this, c);
  };
  std::tie(data_.cell_volumes, data_.cell_centroids) = asDualView(lambda, ncells_all);
  data_.cell_geometry_cached = true;
}

// cell-face adjacencies

void
MeshCacheHost::cacheCellFaces()
{
  static_assert(MEM == MemSpace_kind::HOST);
  assert(framework_mesh_.get());
  if (data_.cell_faces_cached) return;
  int num_cells =
    framework_mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::ALL);
  auto lambda1 = [this](Entity_ID c, MeshFramework::cEntity_ID_View& cfaces) {
    this->framework_mesh_->getCellFaces(c, cfaces);
  };
  data_.cell_faces = asRaggedArray_DualView<Entity_ID>(lambda1, num_cells);
  auto lambda2 = [this](Entity_ID c, MeshFramework::cDirection_View& dirs) {
    this->framework_mesh_->getCellFaceDirs(c, dirs);
  };
  data_.cell_face_directions = asRaggedArray_DualView<int>(lambda2, num_cells);
  data_.cell_faces_cached = true;
}


void
MeshCacheHost::cacheCellCoordinates()
{
  assert(framework_mesh_.get());
  if (data_.cell_coordinates_cached) return;
  int num_cells =
    framework_mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::ALL);
  auto lambda = [this](Entity_ID c,
                       MeshFramework::cPoint_View& ccoords) {
    ccoords = this->framework_mesh_->getCellCoordinates(c);
  };
  data_.cell_coordinates = asRaggedArray_DualView<Amanzi::AmanziGeometry::Point>(lambda, num_cells);
  data_.cell_coordinates_cached = true;
}


void
MeshCacheHost::cacheCellEdges()
{
  assert(framework_mesh_.get());
  if (data_.cell_edges_cached) return;
  framework_mesh_->hasEdgesOrThrow();
  int num_cells =
    framework_mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::ALL);
  auto lambda = [this](Entity_ID c, MeshFramework::cEntity_ID_View& cedges) {
    this->framework_mesh_->getCellEdges(c, cedges);
  };
  data_.cell_edges = asRaggedArray_DualView<Entity_ID>(lambda, num_cells);
  data_.cell_edges_cached = true;
}


void
MeshCacheHost::cacheCellNodes()
{
  assert(framework_mesh_.get());
  if (data_.cell_nodes_cached) return;
  int num_cells =
    framework_mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::ALL);
  auto lambda = [this](Entity_ID c, MeshFramework::cEntity_ID_View& cnodes) {
    this->framework_mesh_->getCellNodes(c, cnodes);
  };
  data_.cell_nodes = asRaggedArray_DualView<Entity_ID>(lambda, num_cells);
  data_.cell_nodes_cached = true;
}


void
MeshCacheHost::cacheFaceGeometry()
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
                           MeshFramework::cDirection_View& dirs) {
    Direction_View ldirs;
    // This NEEDS to call the framework or be passed an host mesh to call the function on the host.
    MeshFramework::cEntity_ID_View fcells;
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
    MeshFramework::cEntity_ID_View cfaces;
    algorithms_->computeCellFacesAndBisectors(*this, c, cfaces, &bisectors);
  };
  data_.cell_face_bisectors = asRaggedArray_DualView<AmanziGeometry::Point>(lambda3, ncells_all);

  // Cache the cell global indices for the function getFaceNormal on device
  auto cgi = getEntityGIDs(Entity_kind::CELL, true);
  data_.cell_global_indices.resize(cgi.size()); 
  Kokkos::deep_copy(data_.cell_global_indices.h_view, cgi);
  Kokkos::deep_copy(data_.cell_global_indices.d_view, cgi); 
  data_.cell_global_indices_cached = true; 
}


// face-cell adjacencies

void
MeshCacheHost::cacheFaceCells()
{
  if (data_.face_cells_cached) return;
  auto lambda = [this](Entity_ID f, MeshFramework::cEntity_ID_View& fcells) {
    this->framework_mesh_->getFaceCells(f, fcells);
  };
  data_.face_cells = asRaggedArray_DualView<Entity_ID>(lambda, nfaces_all);
  data_.face_cells_cached = true;
}

  
// face-cell adjacencies

void
MeshCacheHost::cacheFaceEdges()
{
  if (data_.face_edges_cached) return;
  framework_mesh_->hasEdgesOrThrow();
  auto lambda = [this](Entity_ID f, MeshFramework::cEntity_ID_View& fedges) {
    this->framework_mesh_->getFaceEdges(f, fedges);
  };
  auto lambda_dir = [this](Entity_ID f,
                           MeshFramework::cDirection_View& dirs) {
    MeshFramework::cEntity_ID_View l;
    framework_mesh_->getFaceEdgesAndDirs(f, l, &dirs);
  };
  data_.face_edges = asRaggedArray_DualView<Entity_ID>(lambda, nfaces_all);
  data_.face_edge_directions = asRaggedArray_DualView<int>(lambda_dir, nfaces_all);
  data_.face_edges_cached = true;
}

// face-node adjacencies

void
MeshCacheHost::cacheFaceNodes()
{
  if (data_.face_nodes_cached) return;
  auto lambda = [this](Entity_ID f, MeshFramework::cEntity_ID_View& fnodes) {
    this->framework_mesh_->getFaceNodes(f, fnodes);
  };
  data_.face_nodes = asRaggedArray_DualView<Entity_ID>(lambda, nfaces_all);
  data_.face_nodes_cached = true;
}


void
MeshCacheHost::cacheFaceCoordinates()
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

void
MeshCacheHost::cacheEdgeGeometry()
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

void
MeshCacheHost::cacheEdgeCells()
{
  assert(framework_mesh_.get());
  if (data_.edge_cells_cached) return;
  framework_mesh_->hasEdgesOrThrow();
  auto lambda = [this](Entity_ID f, MeshFramework::cEntity_ID_View& ecells) {
    this->framework_mesh_->getEdgeCells(f, ecells);
  };
  data_.edge_cells = asRaggedArray_DualView<Entity_ID>(lambda, nedges_all);
  data_.edge_cells_cached = true;
}

// edge-face adjacencies

void
MeshCacheHost::cacheEdgeFaces()
{
  assert(framework_mesh_.get());
  if (data_.edge_faces_cached) return;
  framework_mesh_->hasEdgesOrThrow();
  auto lambda = [this](Entity_ID f, MeshFramework::cEntity_ID_View& efaces) {
    this->framework_mesh_->getEdgeFaces(f, efaces);
  };
  data_.edge_faces = asRaggedArray_DualView<Entity_ID>(lambda, nedges_all);
  data_.edge_faces_cached = true;
}

  
// edge-node adjacencies
void
MeshCacheHost::cacheEdgeNodes()
{
  assert(framework_mesh_.get());
  if (data_.edge_nodes_cached) return;
  framework_mesh_->hasEdgesOrThrow();
  auto lambda = [this](Entity_ID f, MeshFramework::cEntity_ID_View& efaces) {
    this->framework_mesh_->getEdgeNodes(f, efaces);
  };
  data_.edge_nodes = asRaggedArray_DualView<Entity_ID>(lambda, nedges_all);
  data_.edge_nodes_cached = true;
}

// edge coordinates

void
MeshCacheHost::cacheEdgeCoordinates()
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

void
MeshCacheHost::cacheNodeCells()
{
  assert(framework_mesh_.get());
  if (data_.node_cells_cached) return;
  auto lambda = [this](Entity_ID f, MeshFramework::cEntity_ID_View& ncells) {
    this->framework_mesh_->getNodeCells(f, ncells);
  };
  data_.node_cells = asRaggedArray_DualView<Entity_ID>(lambda, nnodes_all);
  data_.node_cells_cached = true;
}




// // node-face adjacencies

void
MeshCacheHost::cacheNodeFaces()
{
  assert(framework_mesh_.get());
  if (data_.node_faces_cached) return;
  auto lambda = [this](Entity_ID f, MeshFramework::cEntity_ID_View& nfaces) {
    this->framework_mesh_->getNodeFaces(f, nfaces);
  };
  data_.node_faces = asRaggedArray_DualView<Entity_ID>(lambda, nnodes_all);
  data_.node_faces_cached = true;
}

// // node-edge adjacencies
void
MeshCacheHost::cacheNodeEdges()
{
  assert(framework_mesh_.get());
  if (data_.node_edges_cached) return;
  framework_mesh_->hasEdgesOrThrow();
  auto lambda = [this](Entity_ID f, MeshFramework::cEntity_ID_View& nedges) {
    this->framework_mesh_->getNodeEdges(f, nedges);
  };
  data_.node_edges = asRaggedArray_DualView<Entity_ID>(lambda, nnodes_all);
  data_.node_edges_cached = true;
}


// node coordinates
void
MeshCacheHost::cacheNodeCoordinates()
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
void
MeshCacheHost::cacheParentEntities()
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


void
MeshCacheHost::recacheGeometry()
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


}
