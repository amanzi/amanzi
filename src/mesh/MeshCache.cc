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

#include "MeshCache.hh"
#include "MeshFramework.hh"
#include "MeshMaps_impl.hh"
#include "MeshSets.hh"
#include "Mesh_Helpers.hh"
#include "RegionLabeledSet.hh"

namespace Amanzi {
namespace AmanziMesh {

MeshCache::MeshCache() :
  cell_geometry_cached(false),
  cell_faces_cached(false),
  cell_edges_cached(false),
  cell_nodes_cached(false),
  cell_coordinates_cached(false),
  face_geometry_cached(false),
  face_cells_cached(false),
  face_edges_cached(false),
  face_nodes_cached(false),
  face_coordinates_cached(false),
  edge_geometry_cached(false),
  edge_cells_cached(false),
  edge_faces_cached(false),
  edge_nodes_cached(false),
  edge_coordinates_cached(false),
  node_cells_cached(false),
  node_faces_cached(false),
  node_edges_cached(false),
  node_coordinates_cached(false),
  parent_entities_cached(false),
  is_ordered_(false),
  has_edges_(false),
  has_nodes_(true)
{}


void MeshCache::setParentMesh(const Teuchos::RCP<const MeshCache>& parent)
{
  if (parent_ != Teuchos::null && parent_ != parent) {
    Errors::Message msg("MeshCache::setParentMesh given conflicting parent mesh.");
    Exceptions::amanzi_throw(msg);
  }
  parent_ = parent;
}

void MeshCache::setMeshFramework(const Teuchos::RCP<MeshFramework>& framework_mesh) {
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

  ncells_owned = framework_mesh_->getNumEntities(Entity_kind::CELL, Parallel_type::OWNED);
  ncells_all = framework_mesh_->getNumEntities(Entity_kind::CELL, Parallel_type::ALL);

  nfaces_owned = framework_mesh_->getNumEntities(Entity_kind::FACE, Parallel_type::OWNED);
  nfaces_all = framework_mesh_->getNumEntities(Entity_kind::FACE, Parallel_type::ALL);

  if (has_edges_) {
    nedges_owned = framework_mesh_->getNumEntities(Entity_kind::EDGE, Parallel_type::OWNED);
    nedges_all = framework_mesh_->getNumEntities(Entity_kind::EDGE, Parallel_type::ALL);
  }
  if (has_nodes_) {
    nnodes_owned = framework_mesh_->getNumEntities(Entity_kind::NODE, Parallel_type::OWNED);
    nnodes_all = framework_mesh_->getNumEntities(Entity_kind::NODE, Parallel_type::ALL);
  }

  bool natural_ordered_maps = plist_->get<bool>("natural map ordering", false);
  buildMaps(natural_ordered_maps);

  nboundary_faces_owned = getMap(Entity_kind::BOUNDARY_FACE, false).NumMyElements();
  nboundary_faces_all = getMap(Entity_kind::BOUNDARY_FACE, true).NumMyElements();

  if (has_nodes_) {
    nboundary_nodes_owned = getMap(Entity_kind::BOUNDARY_NODE, false).NumMyElements();
    nboundary_nodes_all = getMap(Entity_kind::BOUNDARY_NODE, true).NumMyElements();
  }
}


bool
MeshCache::isValidSetType(const AmanziGeometry::RegionType rtype, const Entity_kind kind) const {
  if (rtype == AmanziGeometry::RegionType::LABELEDSET && framework_mesh_.get()) {
    return framework_mesh_->isValidSetType(rtype, kind);
  }
  return true;
}


int
MeshCache::getSetSize(const std::string& region_name,
                      const Entity_kind kind,
                      const Parallel_type ptype) const {
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
template<MemSpace_type MEM>
cEntity_ID_View
MeshCache<MEM>::getSetEntities(const std::string& region_name,
        const Entity_kind kind,
        const Parallel_type ptype) const
{
  auto key = std::make_tuple(region_name, kind, ptype);
  if (!sets_.count(key)) {
    auto region = getGeometricModel()->FindRegion(region_name);
    if (region == Teuchos::null) {
      Errors::Message msg;
      msg << "Cannot find region of name \"" << region_name << "\" in the geometric model.";
      Exceptions::amanzi_throw(msg);
    }
    MeshCache<MemSpace_type::HOST> this_on_host(*this);
    sets_[key] = asDualView(resolveMeshSet(*region, kind, ptype, this_on_host));

    // Error on zero -- some zeros already error internally (at the framework
    // level) but others don't.  This is the highest level we can catch these at.
    int lsize = sets_[key].size();
    int gsize = 0;
    getComm()->SumAll(&lsize, &gsize, 1);
    if (gsize == 0) {
      Errors::Message msg;
      msg << "Region \"" << region_name << "\" has no entities on this mesh of type " << to_string(kind);
      Exceptions::amanzi_throw(msg);
    }
  }
  return view<MEM>(sets_.at(key));
}

template<MemSpace_type MEM>
cEntity_ID_View
MeshCache<MEM>::getSetEntitiesAndVolumeFractions(const std::string& region_name,
        const Entity_kind kind,
        const Parallel_type ptype,
        Double_View* vol_fracs) const
{
  if (!vol_fracs) return getSetEntities(region_name, kind, ptype);

  auto key = std::make_tuple(region_name, kind, ptype);
  if (!set_vol_fracs_.count(key)) {
    auto region = getGeometricModel()->FindRegion(region_name);
    if (region == Teuchos::null) {
      Errors::Message msg;
      msg << "Cannot find region of name \"" << region_name << "\" in the geometric model.";
      Exceptions::amanzi_throw(msg);
    }

    Double_List vol_facs_list;
    MeshCache<MemSpace_type::HOST> this_on_host(*this);
    sets_[key] = asDualView(resolveMeshSetVolumeFractions(*region, kind, ptype, vol_fracs_list, this_on_host));
    set_vol_fracs_[key] = asDualView(vol_fracs_list);
  } else {
    *vol_fracs = view<MEM>(set_vol_fracs_.at(key));
  }
  return view<MEM>(sets_.at(key));
}


std::size_t
MeshCache::getNumEntities(const Entity_kind kind, const Parallel_type ptype) const
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
    case (Parallel_type::OWNED) :
      return nowned;
      break;
    case (Parallel_type::ALL) :
      return nall;
      break;
    case Parallel_type::GHOST :
      return nall - nowned;
      break;
    default :
      return 0;
  }
}


Cell_type MeshCache::getCellType(const Entity_ID c) const
{
  return MeshAlgorithms::getCellType(*this, c);
}

Parallel_type
MeshCache::getParallelType(const Entity_kind& kind, const Entity_ID id) const
{
  if (id < getNumEntities(kind, Parallel_type::OWNED)) {
    return Parallel_type::OWNED;
  } else if (id < getNumEntities(kind, Parallel_type::ALL)) {
    return Parallel_type::GHOST;
  }
  return Parallel_type::UNKNOWN;
}



//
// Build the cache, fine grained control
// =============================================

// cell centroid, volume
void MeshCache::cacheCellGeometry()
{
  if (cell_geometry_cached) return;
  cell_volumes.resize(ncells_all);
  cell_centroids.resize(ncells_all);
  for (Entity_ID i=0; i!=ncells_all; ++i) {
    std::tie(cell_volumes[i], cell_centroids[i]) = framework_mesh_->computeCellGeometry(i);
  }
  cell_geometry_cached = true;
}

// cell-face adjacencies
void MeshCache::cacheCellFaces()
{
  if (cell_faces_cached) return;
  cell_faces.resize(ncells_all);
  cell_face_directions.resize(ncells_all);
  for (Entity_ID i=0; i!=ncells_all; ++i) {
    framework_mesh_->getCellFacesAndDirs(i, cell_faces[i], &cell_face_directions[i]);
  }
  cell_faces_cached = true;
}

// cell-edge adjacencies
void MeshCache::cacheCellEdges()
{
  if (cell_edges_cached) return;
  framework_mesh_->hasEdgesOrThrow();
  cell_edges.resize(ncells_all);
  for (Entity_ID i=0; i!=ncells_all; ++i) {
    framework_mesh_->getCellEdges(i, cell_edges[i]);
  }
  cell_edges_cached = true;
}

// cell-node adjacencies
void MeshCache::cacheCellNodes()
{
  if (cell_nodes_cached) return;
  cell_nodes.resize(ncells_all);
  for (Entity_ID i=0; i!=ncells_all; ++i) {
    framework_mesh_->getCellNodes(i, cell_nodes[i]);
  }
  cell_nodes_cached = true;
}

// cell coordinates
void MeshCache::cacheCellCoordinates()
{
  if (cell_coordinates_cached) return;
  cell_coordinates.resize(ncells_all);
  for (Entity_ID i=0; i!=ncells_all; ++i) {
    cell_coordinates[i] = framework_mesh_->getCellCoordinates(i);
  }
  cell_coordinates_cached = true;
}


// face centroid, area, normal
void MeshCache::cacheFaceGeometry()
{
  if (face_geometry_cached) return;
  double area;
  face_centroids.resize(nfaces_all);
  face_areas.resize(nfaces_all);
  face_normals.resize(nfaces_all);
  for (Entity_ID i=0; i!=nfaces_all; ++i) {
    std::tie(face_areas[i], face_centroids[i], face_normals[i]) = framework_mesh_->computeFaceGeometry(i);
  }

  // this could get its own call...
  face_normal_directions.resize(nfaces_all);
  AMANZI_ASSERT(face_cells_cached);
  for (Entity_ID f=0; f!=face_cells.size(); ++f) {
    face_normal_directions[f].resize(face_cells[f].size());
    for (Entity_ID i=0; i!=face_cells[f].size(); ++i) {
      face_normals[f][i] = framework_mesh_->getFaceNormal(f, face_cells[f][i], &face_normal_directions[f][i]);
    }
  }

  // this could get its own call...
  cell_face_bisectors.resize(ncells_all);
  Entity_ID_List faces;
  for (Entity_ID i=0; i!=ncells_all; ++i) {
    framework_mesh_->getCellFacesAndBisectors(i, faces, &cell_face_bisectors[i]);
  }
  face_geometry_cached = true;
}

// face-cell adjacencies
void MeshCache::cacheFaceCells()
{
  if (face_cells_cached) return;
  face_cells.resize(nfaces_all);
  for (Entity_ID i=0; i!=nfaces_all; ++i) {
    framework_mesh_->getFaceCells(i, Parallel_type::ALL, face_cells[i]);
  }
  face_cells_cached = true;
}

// face-edge adjacencies
void MeshCache::cacheFaceEdges()
{
  if (face_edges_cached) return;
  framework_mesh_->hasEdgesOrThrow();
  face_edges.resize(nfaces_all);
  face_edge_directions.resize(nfaces_all);
  for (Entity_ID i=0; i!=nfaces_all; ++i) {
    framework_mesh_->getFaceEdgesAndDirs(i, face_edges[i], &face_edge_directions[i]);
  }
  face_edges_cached = true;
}

// face-node adjacencies
void MeshCache::cacheFaceNodes()
{
  if (face_nodes_cached) return;
  face_nodes.resize(nfaces_all);
  for (Entity_ID i=0; i!=nfaces_all; ++i) {
    framework_mesh_->getFaceNodes(i, face_nodes[i]);
  }
  face_nodes_cached = true;
}

// face coordinates
void MeshCache::cacheFaceCoordinates()
{
  if (face_coordinates_cached) return;
  face_coordinates.resize(nfaces_all);
  for (Entity_ID i=0; i!=nfaces_all; ++i) {
    face_coordinates[i] = framework_mesh_->getFaceCoordinates(i);
  }
  face_coordinates_cached = true;
}


// edge centroid, length, vector
void MeshCache::cacheEdgeGeometry()
{
  if (edge_geometry_cached) return;
  framework_mesh_->hasEdgesOrThrow();

  edge_centroids.resize(nedges_all);
  edge_vectors.resize(nedges_all);
  for (Entity_ID i=0; i!=nedges_all; ++i) {
    std::tie(edge_vectors[i], edge_centroids[i]) = framework_mesh_->computeEdgeGeometry(i);
  }
  edge_geometry_cached = true;
}

// edge-cell adjacencies
void MeshCache::cacheEdgeCells()
{
  if (edge_cells_cached) return;
  framework_mesh_->hasEdgesOrThrow();
  edge_cells.resize(nedges_all);
  for (Entity_ID i=0; i!=nedges_all; ++i) {
    framework_mesh_->getEdgeCells(i, Parallel_type::ALL, edge_cells[i]);
  }
  edge_cells_cached = true;
}

// edge-face adjacencies
void MeshCache::cacheEdgeFaces()
{
  if (edge_faces_cached) return;
  framework_mesh_->hasEdgesOrThrow();
  edge_faces.resize(nedges_all);
  for (Entity_ID i=0; i!=nedges_all; ++i) {
    framework_mesh_->getEdgeFaces(i, Parallel_type::ALL, edge_faces[i]);
  }
  edge_faces_cached = true;
}

// edge-node adjacencies
void MeshCache::cacheEdgeNodes()
{
  if (edge_nodes_cached) return;
  framework_mesh_->hasEdgesOrThrow();
  edge_nodes.resize(nedges_all);
  for (Entity_ID i=0; i!=nedges_all; ++i) {
    framework_mesh_->getEdgeNodes(i, edge_nodes[i]);
  }
  edge_nodes_cached = true;
}

// edge coordinates
void MeshCache::cacheEdgeCoordinates()
{
  if (edge_coordinates_cached) return;
  edge_coordinates.resize(nedges_all);
  for (Entity_ID i=0; i!=nedges_all; ++i) {
    edge_coordinates[i] = framework_mesh_->getEdgeCoordinates(i);
  }
  edge_coordinates_cached = true;
}


// node-cell adjacencies
void MeshCache::cacheNodeCells()
{
  if (node_cells_cached) return;
  node_cells.resize(nnodes_all);
  for (Entity_ID i=0; i!=nnodes_all; ++i) {
    framework_mesh_->getNodeCells(i, Parallel_type::ALL, node_cells[i]);
  }
  node_cells_cached = true;
}

// node-face adjacencies
void MeshCache::cacheNodeFaces()
{
  if (node_faces_cached) return;
  node_faces.resize(nnodes_all);
  for (Entity_ID i=0; i!=nnodes_all; ++i) {
    framework_mesh_->getNodeFaces(i, Parallel_type::ALL, node_faces[i]);
  }
  node_faces_cached = true;
}

// node-edge adjacencies
void MeshCache::cacheNodeEdges()
{
  if (node_edges_cached) return;
  framework_mesh_->hasEdgesOrThrow();
  node_edges.resize(nnodes_all);
  for (Entity_ID i=0; i!=nnodes_all; ++i) {
    framework_mesh_->getNodeEdges(i, Parallel_type::ALL, node_edges[i]);
  }
  node_edges_cached = true;
}

// node coordinates
void MeshCache::cacheNodeCoordinates()
{
  if (node_coordinates_cached) return;
  node_coordinates.resize(nnodes_all);
  for (Entity_ID i=0; i!=nnodes_all; ++i) {
    node_coordinates[i] = framework_mesh_->getNodeCoordinate(i);
  }
  node_coordinates_cached = true;
}


// Note that regions are cached on demand the first time they are requested,
// but labeled sets must be pre-cached if the framework mesh is to be
// destroyed.
void MeshCache::precacheLabeledSets()
{
  for (const auto& rgn : *getGeometricModel()) {
    if (rgn->get_type() == AmanziGeometry::RegionType::LABELEDSET) {
      auto rgn_lbl = Teuchos::rcp_dynamic_cast<const AmanziGeometry::RegionLabeledSet>(rgn);
      AMANZI_ASSERT(rgn_lbl);

      if (getParentMesh() != Teuchos::null) {
        AMANZI_ASSERT(false); // not yet implemented lifted sets
      } else {
        auto entity_kind = createEntityKind(rgn_lbl->entity_str());
        Entity_ID_List set_ents;
        framework_mesh_->getSetEntities(*rgn_lbl, entity_kind, Parallel_type::ALL, set_ents);
        auto key = std::make_tuple(rgn->get_name(), entity_kind, Parallel_type::ALL);
        sets_[key] = asDualView(set_ents);
      }
    }
  }
}

// build the MeshMaps object
void MeshCache::buildMaps(bool natural)
{
  maps_.initialize(*framework_mesh_, natural);
}


void MeshCache::setNodeCoordinate(const Entity_ID n, const AmanziGeometry::Point& coord)
{
  if (framework_mesh_.get())
    framework_mesh_->setNodeCoordinate(n, coord);
  if (node_coordinates_cached)
    node_coordinates[n] = coord;
}

// common error messaging
void MeshCache::throwAccessError_(const std::string& func_name) const
{
  Errors::Message msg;
  msg << "MeshCache::" << func_name << " cannot compute this quantity -- not cached and framework does not exist.";
  Exceptions::amanzi_throw(msg);
}


Entity_ID
MeshCache::getEntityParent(const Entity_kind kind, const Entity_ID entid) const
{
  if (parent_entities_cached) {
    switch(kind) {
      case Entity_kind::CELL:
        return parent_cells[entid];
        break;
      case Entity_kind::FACE:
        return parent_faces[entid];
        break;
      case Entity_kind::EDGE:
        return parent_edges[entid];
        break;
      case Entity_kind::NODE:
        return parent_nodes[entid];
      default: {}
    }
  } else if (framework_mesh_.get()) {
    return framework_mesh_->getEntityParent(kind, entid);
  }
  return -1;
}




namespace MeshAlgorithms {

void cacheAll(MeshCache& mesh)
{
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
  mesh.cacheFaceGeometry();
  mesh.cacheCellGeometry();

  if (mesh.hasEdges()) {
    mesh.cacheCellEdges();
    mesh.cacheEdgeCells();
    mesh.cacheFaceEdges();
    mesh.cacheEdgeFaces();
    mesh.cacheNodeEdges();
    mesh.cacheEdgeNodes();
    mesh.cacheEdgeGeometry();
  }
}

void cacheDefault(MeshCache& mesh)
{
  // caches what the developers currently think is best
  mesh.cacheNodeCoordinates();

  mesh.cacheCellFaces();
  mesh.cacheFaceCells();
  mesh.cacheFaceGeometry();
  mesh.cacheCellGeometry();

  // if (mesh.hasEdges()) {
  //   mesh.cacheFaceEdges();
  //   mesh.cacheEdgeFaces();
  //   mesh.cacheEdgeGeometry();
  // }
}


void recacheGeometry(MeshCache& mesh)
{
  // recaches the geometry, as presumably the nodal coordinates have changed.
  if (mesh.edge_coordinates_cached) {
    mesh.edge_coordinates_cached = false;
    mesh.cacheEdgeCoordinates();
  }
  if (mesh.edge_geometry_cached) {
    mesh.edge_geometry_cached = false;
    mesh.cacheEdgeGeometry();
  }
  if(mesh.face_coordinates_cached) {
    mesh.face_coordinates_cached = false;
    mesh.cacheFaceCoordinates();
  }
  if (mesh.face_geometry_cached) {
    mesh.face_geometry_cached = false;
    mesh.cacheFaceGeometry();
  }
  if (mesh.cell_coordinates_cached) {
    mesh.cell_coordinates_cached = false;
    mesh.cacheCellCoordinates();
  }
  if (mesh.cell_geometry_cached) {
    mesh.cell_geometry_cached = false;
    mesh.cacheCellGeometry();
  }
}


} // namespace MeshAlgorithms
} // namespace AmanziMesh
} // namespace Amanzi
