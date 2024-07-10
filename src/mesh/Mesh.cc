/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
           Julien Loiseau (jloiseau@lanl.gov)
           Rao Garimella (rao@lanl.gov)
*/

#include "Teuchos_CommHelpers.hpp"

#include "Geometry.hh"
#include "Iterators.hh"
#include "MeshUtils.hh"
#include "MeshFramework.hh"

#include "MeshInternals.hh"
#include "Mesh.hh"

#include "MeshCache.hh"
#include "MeshAlgorithms.hh"
#include "MeshSets.hh"


namespace Amanzi {
namespace AmanziMesh {

// ----------------------------------
// Constructors
// ----------------------------------
Mesh::Mesh(const Teuchos::RCP<Teuchos::ParameterList>& plist)
  : plist_(plist),
    maps_(Teuchos::rcp(new MeshMaps())),
    sets_(Teuchos::rcp(new MeshSets())),
    set_vol_fracs_(Teuchos::rcp(new MeshSetVolumeFractions()))
{
  if (plist_ == Teuchos::null) {
    plist_ = Teuchos::rcp(new Teuchos::ParameterList());
  }
}


Mesh::Mesh(const Teuchos::RCP<MeshFramework>& framework_mesh,
           const Teuchos::RCP<MeshAlgorithms>& algorithms,
           const Teuchos::RCP<Teuchos::ParameterList>& plist)
  : Mesh(plist)
{
  algorithms_ = algorithms;
  setMeshFramework(framework_mesh);
}


// ---------------------------------------
// Build the cache, fine grained control
// ---------------------------------------
//
// cache cell quantities
//
void
Mesh::cacheCellGeometry()
{
  assert(framework_mesh_.get());
  if (data_.cell_geometry_cached) return;
  auto lambda = [this](Entity_ID c, Double_View cvol, Point_View ccent) {
    std::tie(cvol[c], ccent[c]) = algorithms_->computeCellGeometry(*this, c);
  };
  std::tie(data_.cell_volumes, data_.cell_centroids) = asDualView(lambda, data_.ncells_all);
  data_.cell_geometry_cached = true;
}


void
Mesh::cacheCellFaces()
{
  static_assert(MEM == MemSpace_kind::HOST);
  assert(framework_mesh_.get());
  if (data_.cell_faces_cached) return;
  int num_cells =
    framework_mesh_->getNumEntities(Entity_kind::CELL, Parallel_kind::ALL);
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
Mesh::cacheCellEdges()
{
  assert(framework_mesh_.get());
  if (data_.cell_edges_cached) return;
  framework_mesh_->hasEdgesOrThrow();
  auto lambda = [this](Entity_ID c, MeshFramework::cEntity_ID_View& cedges) {
    this->framework_mesh_->getCellEdges(c, cedges);
  };
  data_.cell_edges = asRaggedArray_DualView<Entity_ID>(lambda, data_.ncells_all);
  data_.cell_edges_cached = true;
}


void
Mesh::cacheCellNodes()
{
  assert(framework_mesh_.get());
  if (data_.cell_nodes_cached) return;
  int num_cells =
    framework_mesh_->getNumEntities(Entity_kind::CELL, Parallel_kind::ALL);
  auto lambda = [this](Entity_ID c, MeshFramework::cEntity_ID_View& cnodes) {
    this->framework_mesh_->getCellNodes(c, cnodes);
  };
  data_.cell_nodes = asRaggedArray_DualView<Entity_ID>(lambda, num_cells);
  data_.cell_nodes_cached = true;
}


void
Mesh::cacheCellCoordinates()
{
  assert(framework_mesh_.get());
  if (data_.cell_coordinates_cached) return;
  auto lambda = [this](Entity_ID c, MeshFramework::cPoint_View& ccoords) {
    ccoords = this->framework_mesh_->getCellCoordinates(c);
  };
  data_.cell_coordinates = asRaggedArray_DualView<Amanzi::AmanziGeometry::Point>(lambda, data_.ncells_all);
  data_.cell_coordinates_cached = true;
}


//
// cache face quantities
//
void
Mesh::cacheFaceGeometry()
{
  if (data_.face_geometry_cached) return;
  assert(framework_mesh_.get());
  data_.face_areas.resize(data_.nfaces_all);
  data_.face_centroids.resize(data_.nfaces_all);

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
  data_.face_normals = asRaggedArray_DualView<AmanziGeometry::Point>(lambda, data_.nfaces_all);

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
  auto lambda2 = [&, this](const Entity_ID& f, MeshFramework::cDirection_View& dirs) {
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
  data_.face_normal_orientations = asRaggedArray_DualView<int>(lambda2, data_.nfaces_all);

  // cache cell-face-bisectors -- make this a separate call?  Think about
  // granularity here.
  auto lambda3 = [&, this](const Entity_ID& c, cPoint_View& bisectors) {
    MeshFramework::cEntity_ID_View cfaces;
    algorithms_->computeCellFacesAndBisectors(*this, c, cfaces, &bisectors);
  };
  data_.cell_face_bisectors = asRaggedArray_DualView<AmanziGeometry::Point>(lambda3, data_.ncells_all);

  // Cache the cell global indices for the function getFaceNormal on device
  auto cgi = getEntityGIDs(Entity_kind::CELL, true);
  data_.cell_global_indices.resize(cgi.size());
  Kokkos::deep_copy(data_.cell_global_indices.h_view, cgi);
  Kokkos::deep_copy(data_.cell_global_indices.d_view, cgi);
  data_.cell_global_indices_cached = true;
}


void
Mesh::cacheFaceCells()
{
  if (data_.face_cells_cached) return;
  auto lambda = [this](Entity_ID f, MeshFramework::cEntity_ID_View& fcells) {
    this->framework_mesh_->getFaceCells(f, fcells);
  };
  data_.face_cells = asRaggedArray_DualView<Entity_ID>(lambda, data_.nfaces_all);
  data_.face_cells_cached = true;
}


void
Mesh::cacheFaceEdges()
{
  if (data_.face_edges_cached) return;
  framework_mesh_->hasEdgesOrThrow();
  auto lambda = [this](Entity_ID f, MeshFramework::cEntity_ID_View& fedges) {
    this->framework_mesh_->getFaceEdges(f, fedges);
  };
  auto lambda_dir = [this](Entity_ID f, MeshFramework::cDirection_View& dirs) {
    MeshFramework::cEntity_ID_View l;
    framework_mesh_->getFaceEdgesAndDirs(f, l, &dirs);
  };
  data_.face_edges = asRaggedArray_DualView<Entity_ID>(lambda, data_.nfaces_all);
  data_.face_edge_directions = asRaggedArray_DualView<int>(lambda_dir, data_.nfaces_all);
  data_.face_edges_cached = true;
}


void
Mesh::cacheFaceNodes()
{
  if (data_.face_nodes_cached) return;
  auto lambda = [this](Entity_ID f, MeshFramework::cEntity_ID_View& fnodes) {
    this->framework_mesh_->getFaceNodes(f, fnodes);
  };
  data_.face_nodes = asRaggedArray_DualView<Entity_ID>(lambda, data_.nfaces_all);
  data_.face_nodes_cached = true;
}


void
Mesh::cacheFaceCoordinates()
{
  assert(framework_mesh_.get());
  if (data_.face_coordinates_cached) return;
  auto lambda = [this](Entity_ID f, cPoint_View& fcoords) {
    fcoords = this->framework_mesh_->getFaceCoordinates(f);
  };
  data_.face_coordinates =
    asRaggedArray_DualView<Amanzi::AmanziGeometry::Point>(lambda, data_.nfaces_all);

  data_.face_coordinates_cached = true;
}


//
// cache edge quantities
//
// edge centroid, length, vector
void
Mesh::cacheEdgeGeometry()
{
  assert(framework_mesh_.get());
  if (data_.edge_geometry_cached) return;
  framework_mesh_->hasEdgesOrThrow();
  auto lambda = [this](Entity_ID e, Point_View ev, Point_View ec) {
    std::tie(ev[e], ec[e]) = algorithms_->computeEdgeGeometry(*this, e);
  };
  std::tie(data_.edge_vectors, data_.edge_centroids) = asDualView(lambda, data_.nedges_all);
  data_.edge_geometry_cached = true;
}


void
Mesh::cacheEdgeCells()
{
  assert(framework_mesh_.get());
  if (data_.edge_cells_cached) return;
  framework_mesh_->hasEdgesOrThrow();
  auto lambda = [this](Entity_ID f, MeshFramework::cEntity_ID_View& ecells) {
    this->framework_mesh_->getEdgeCells(f, ecells);
  };
  data_.edge_cells = asRaggedArray_DualView<Entity_ID>(lambda, data_.nedges_all);
  data_.edge_cells_cached = true;
}


void
Mesh::cacheEdgeFaces()
{
  assert(framework_mesh_.get());
  if (data_.edge_faces_cached) return;
  framework_mesh_->hasEdgesOrThrow();
  auto lambda = [this](Entity_ID f, MeshFramework::cEntity_ID_View& efaces) {
    this->framework_mesh_->getEdgeFaces(f, efaces);
  };
  data_.edge_faces = asRaggedArray_DualView<Entity_ID>(lambda, data_.nedges_all);
  data_.edge_faces_cached = true;
}


void
Mesh::cacheEdgeNodes()
{
  assert(framework_mesh_.get());
  if (data_.edge_nodes_cached) return;
  framework_mesh_->hasEdgesOrThrow();
  auto lambda = [this](Entity_ID f, MeshFramework::cEntity_ID_View& efaces) {
    this->framework_mesh_->getEdgeNodes(f, efaces);
  };
  data_.edge_nodes = asRaggedArray_DualView<Entity_ID>(lambda, data_.nedges_all);
  data_.edge_nodes_cached = true;
}


void
Mesh::cacheEdgeCoordinates()
{
  assert(framework_mesh_.get());
  if (data_.edge_coordinates_cached) return;
  framework_mesh_->hasEdgesOrThrow();
  auto lambda = [this](Entity_ID f, cPoint_View& enodes) {
    enodes = this->framework_mesh_->getEdgeCoordinates(f);
  };
  data_.edge_coordinates =
    asRaggedArray_DualView<Amanzi::AmanziGeometry::Point>(lambda, data_.nedges_all);
  data_.edge_coordinates_cached = true;
}


//
// cache node quantities
//
void
Mesh::cacheNodeCells()
{
  assert(framework_mesh_.get());
  if (data_.node_cells_cached) return;
  auto lambda = [this](Entity_ID f, MeshFramework::cEntity_ID_View& ncells) {
    this->framework_mesh_->getNodeCells(f, ncells);
  };
  data_.node_cells = asRaggedArray_DualView<Entity_ID>(lambda, data_.nnodes_all);
  data_.node_cells_cached = true;
}


void
Mesh::cacheNodeFaces()
{
  assert(framework_mesh_.get());
  if (data_.node_faces_cached) return;
  auto lambda = [this](Entity_ID f, MeshFramework::cEntity_ID_View& nfaces) {
    this->framework_mesh_->getNodeFaces(f, nfaces);
  };
  data_.node_faces = asRaggedArray_DualView<Entity_ID>(lambda, data_.nnodes_all);
  data_.node_faces_cached = true;
}


void
Mesh::cacheNodeEdges()
{
  assert(framework_mesh_.get());
  if (data_.node_edges_cached) return;
  framework_mesh_->hasEdgesOrThrow();
  auto lambda = [this](Entity_ID f, MeshFramework::cEntity_ID_View& nedges) {
    this->framework_mesh_->getNodeEdges(f, nedges);
  };
  data_.node_edges = asRaggedArray_DualView<Entity_ID>(lambda, data_.nnodes_all);
  data_.node_edges_cached = true;
}


void
Mesh::cacheNodeCoordinates()
{
  if (data_.node_coordinates_cached) return;
  data_.node_coordinates.resize(data_.nnodes_all);
  for (Entity_ID i = 0; i != data_.nnodes_all; ++i) {
    view<MemSpace_kind::HOST>(data_.node_coordinates)[i] = framework_mesh_->getNodeCoordinate(i);
  }

  Kokkos::deep_copy(data_.node_coordinates.d_view, data_.node_coordinates.h_view);
  data_.node_coordinates_cached = true;
}


//
// cache parent entities
//
void
Mesh::cacheParentEntities()
{
  if (data_.parent_entities_cached) return;
  assert(framework_mesh_.get());

  data_.parent_cells.resize(data_.ncells_all);
  for (Entity_ID i = 0; i != data_.ncells_all; ++i) {
    view<MemSpace_kind::HOST>(data_.parent_cells)[i] =
      framework_mesh_->getEntityParent(Entity_kind::CELL, i);
  }
  Kokkos::deep_copy(data_.parent_cells.d_view, data_.parent_cells.h_view);

  data_.parent_faces.resize(data_.nfaces_all);
  for (Entity_ID i = 0; i != data_.nfaces_all; ++i) {
    view<MemSpace_kind::HOST>(data_.parent_faces)[i] =
      framework_mesh_->getEntityParent(Entity_kind::FACE, i);
  }
  Kokkos::deep_copy(data_.parent_faces.d_view, data_.parent_faces.h_view);

  data_.parent_nodes.resize(data_.nnodes_all);
  for (Entity_ID i = 0; i != data_.nnodes_all; ++i) {
    view<MemSpace_kind::HOST>(data_.parent_nodes)[i] =
      framework_mesh_->getEntityParent(Entity_kind::NODE, i);
  }
  Kokkos::deep_copy(data_.parent_nodes.d_view, data_.parent_nodes.h_view);

  data_.parent_entities_cached = true;
}


void
Mesh::cacheDefault()
{
  // topology
  cacheCellFaces();
  cacheFaceCells();
  if (hasEdges()) cacheFaceEdges();
  if (hasNodes()) cacheFaceNodes();

  if (hasEdges()) {
    cacheEdgeFaces();
    if (hasNodes()) cacheEdgeNodes();
  }

  if (hasNodes()) {
    cacheNodeFaces();
    if (hasEdges()) cacheNodeEdges();
  }

  // geometry
  if (hasNodes()) {
    cacheNodeCoordinates();
    // if (hasEdges()) cacheEdgeCoordinates();
    // cacheFaceCoordinates();
    // cacheCellCoordinates();
  }

  if (hasEdges()) cacheEdgeGeometry();  // edge centroid, length, vector
  cacheFaceGeometry(); // centroid, area, normal
  cacheCellGeometry();  // cell centroid, volume

  syncCache();
}


void
Mesh::cacheAll()
{
  // topology
  cacheCellFaces();
  if (hasEdges()) cacheCellEdges();
  if (hasNodes()) cacheCellNodes();

  cacheFaceCells();
  if (hasEdges()) cacheFaceEdges();
  if (hasNodes()) cacheFaceNodes();

  if (hasEdges()) {
    cacheEdgeCells();
    cacheEdgeFaces();
    if (hasNodes()) cacheEdgeNodes();
  }

  if (hasNodes()) {
    cacheNodeCells();
    cacheNodeFaces();
    if (hasEdges()) cacheNodeEdges();
  }

  // geometry
  if (hasNodes()) {
    cacheNodeCoordinates();
    if (hasEdges()) cacheEdgeCoordinates();
    cacheFaceCoordinates();
    cacheCellCoordinates();
  }

  if (hasEdges()) cacheEdgeGeometry();  // edge centroid, length, vector
  cacheFaceGeometry(); // centroid, area, normal
  cacheCellGeometry();  // cell centroid, volume

  // sync host to device
  syncCache();
}


void
Mesh::recacheGeometry()
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


void
Mesh::syncCache()
{
  if (cache_ == Teuchos::null) cache_ = Teuchos::rcp(new MeshCache(data_));
  else cache_->data = data_;
}


// Columnar semi-structured meshes
void Mesh::buildColumns()
{
  if (columns == Teuchos::null) columns = Teuchos::rcp(new MeshColumns());
  columns->initialize(*this);
  cache_->columns = *columns;
}

void Mesh::buildColumns(const std::vector<std::string>& regions)
{
  if (columns == Teuchos::null) columns = Teuchos::rcp(new MeshColumns());
  columns->initialize(*this, regions);
  cache_->columns = *columns;
}


// ----------------------------------
// Accessors and Mutators
// ----------------------------------
void
Mesh::setMeshFramework(const Teuchos::RCP<MeshFramework>& framework_mesh)
{
  framework_mesh_ = framework_mesh;
  comm_ = framework_mesh_->getComm();
  gm_ = framework_mesh_->getGeometricModel();

  data_.has_edges_ = framework_mesh->hasEdges();
  data_.has_nodes_ = framework_mesh->hasNodes();
  data_.has_node_faces_ = framework_mesh->hasNodeFaces();
  data_.space_dim_ = framework_mesh_->getSpaceDimension();
  data_.manifold_dim_ = framework_mesh_->getManifoldDimension();
  data_.is_logical_ = framework_mesh_->isLogical();
  data_.is_ordered_ = framework_mesh_->isOrdered();

  data_.ncells_owned = framework_mesh_->getNumEntities(Entity_kind::CELL, Parallel_kind::OWNED);
  data_.ncells_all = framework_mesh_->getNumEntities(Entity_kind::CELL, Parallel_kind::ALL);
  data_.nfaces_owned = framework_mesh_->getNumEntities(Entity_kind::FACE, Parallel_kind::OWNED);
  data_.nfaces_all = framework_mesh_->getNumEntities(Entity_kind::FACE, Parallel_kind::ALL);
  if (data_.has_edges_) {
    data_.nedges_owned = framework_mesh_->getNumEntities(Entity_kind::EDGE, Parallel_kind::OWNED);
    data_.nedges_all = framework_mesh_->getNumEntities(Entity_kind::EDGE, Parallel_kind::ALL);
  }
  if (data_.has_nodes_) {
    data_.nnodes_owned = framework_mesh_->getNumEntities(Entity_kind::NODE, Parallel_kind::OWNED);
    data_.nnodes_all = framework_mesh_->getNumEntities(Entity_kind::NODE, Parallel_kind::ALL);
  }
  bool natural_ordered_maps = plist_->get<bool>("natural map ordering", false);
  maps_->initialize(*framework_mesh_, natural_ordered_maps);

  data_.nboundary_faces_owned = maps_->getNBoundaryFaces(Parallel_kind::OWNED);
  data_.nboundary_faces_all = maps_->getNBoundaryFaces(Parallel_kind::ALL);
  if (data_.has_nodes_) {
    data_.nboundary_nodes_owned = maps_->getNBoundaryNodes(Parallel_kind::OWNED);
    data_.nboundary_nodes_all = maps_->getNBoundaryNodes(Parallel_kind::ALL);
  }

  data_.boundary_faces = maps_->boundary_faces;
  data_.boundary_nodes = maps_->boundary_nodes;

  std::string policy = plist_->get<std::string>("cache policy", "default");
  if (policy == "all") {
    cacheAll();
  } else if (policy == "default") {
    cacheDefault();
  } else {
    Errors::Message msg;
    msg << "MeshCache construction: unknown \"cache policy\" = \"" << policy
        << "\", must be one of \"default\" or \"all\"";
    Exceptions::amanzi_throw(msg);
  }
}


void
Mesh::setParentMesh(const Teuchos::RCP<const Mesh>& parent)
{
  if (parent_ != Teuchos::null && parent_ != parent) {
    Errors::Message msg("Mesh::setParentMesh given conflicting parent mesh.");
    Exceptions::amanzi_throw(msg);
  } else {
    parent_ = parent;
    cacheParentEntities();
    syncCache();
  }
}


const Mesh&
Mesh::getVisMesh() const
{
  if (vis_mesh_.get()) return *vis_mesh_;
  return *this;
}


void
Mesh::printMeshStatistics() const
{
  auto vo_ = Teuchos::rcp(new VerboseObject("Mesh Output", *plist_));
  if (vo_.get() && vo_->getVerbLevel() >= Teuchos::VERB_LOW) {
    int ncells = getNumEntities(Entity_kind::CELL, Parallel_kind::OWNED);
    int nfaces = getNumEntities(Entity_kind::FACE, Parallel_kind::OWNED);
    int nnodes = getNumEntities(Entity_kind::NODE, Parallel_kind::OWNED);
    int nedges(0);
    if (data_.has_edges_) nedges = getNumEntities(Entity_kind::EDGE, Parallel_kind::OWNED);

    int min_out[4], max_out[4], sum_out[4], tmp_in[4] = { ncells, nfaces, nedges, nnodes };
    Teuchos::reduceAll(*getComm(), Teuchos::REDUCE_MIN, 4, tmp_in, min_out);
    Teuchos::reduceAll(*getComm(), Teuchos::REDUCE_MAX, 4, tmp_in, max_out);
    Teuchos::reduceAll(*getComm(), Teuchos::REDUCE_SUM, 4, tmp_in, sum_out);

    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "cells, tot/min/max: " << sum_out[0] << "/" << min_out[0] << "/" << max_out[0]
               << "\n";
    *vo_->os() << "faces, tot/min/max: " << sum_out[1] << "/" << min_out[1] << "/" << max_out[1]
               << "\n";
    if (data_.has_edges_)
      *vo_->os() << "edges, tot/min/max: " << sum_out[2] << "/" << min_out[2] << "/" << max_out[2]
                 << "\n";
    *vo_->os() << "nodes, tot/min/max: " << sum_out[3] << "/" << min_out[3] << "/" << max_out[3]
               << "\n\n";
  }
}


// -------------------
// Map objects
// -------------------
//
// maps define GIDs of each Entity_kind
const Map_ptr_type&
Mesh::getMap(const Entity_kind kind, bool is_ghosted) const
{
  return maps_->getMap(kind, is_ghosted);
}

Entity_GID
Mesh::getEntityGID(const Entity_kind kind, const Entity_ID lid) const
{
  return getMap(kind, true)->getGlobalElement(lid);
}

typename Mesh::cEntity_GID_View
Mesh::getEntityGIDs(const Entity_kind kind, bool ghosted) const
{
  return getMap(kind, ghosted)->getMyGlobalIndices();
}

Entity_ID
Mesh::getEntityLID(const Entity_kind kind, const Entity_GID gid, bool ghosted) const
{
  return getMap(kind, ghosted)->getLocalElement(gid);
}

// importers allow scatter/gather operations
const Import_type&
Mesh::getImporter(const Entity_kind kind) const
{
  return maps_->getImporter(kind);
}

// an importer from FACE-indexed objects to BOUNDARY_FACE-indexed objects
//
// Note this is not the same as getImporter(BOUNDARY_FACE), which
// communicates BOUNDARY_FACE-indexed objects to other BOUNDARY_FACE-indexed
// objects.
const Import_type&
Mesh::getBoundaryFaceImporter() const
{
  return maps_->getBoundaryFaceImporter();
}

// an importer from NODE-indexed objects to BOUNDARY_NODE-indexed objects
const Import_type&
Mesh::getBoundaryNodeImporter() const
{
  return maps_->getBoundaryNodeImporter();
}

// an importer from CELL-indexed objects to BOUNDARY_FACE-indexed objects
const Import_type&
Mesh::getBoundaryFaceInternalCellImporter() const
{
  return maps_->getBoundaryFaceInternalCellImporter();
}


// ----------------
// sets of entities
// ----------------
bool
Mesh::isValidSetType(const AmanziGeometry::RegionType rtype, const Entity_kind kind) const
{
  if (rtype == AmanziGeometry::RegionType::LABELEDSET && framework_mesh_.get()) {
    return framework_mesh_->isValidSetType(rtype, kind);
  }
  return true;
}


int
Mesh::getSetSize(const std::string& region_name,
                          const Entity_kind kind,
                          const Parallel_kind ptype) const
{
  return getSetEntities<MEM>(region_name, kind, ptype).size();
}


// ----------------
// Entity meta-data
// ----------------
Entity_ID
Mesh::getNumEntities(const Entity_kind kind, const Parallel_kind ptype) const
{
  Entity_ID nowned, nall;
  switch (kind) {
  case (Entity_kind::CELL):
    nowned = data_.ncells_owned;
    nall = data_.ncells_all;
    break;
  case (Entity_kind::FACE):
    nowned = data_.nfaces_owned;
    nall = data_.nfaces_all;
    break;
  case (Entity_kind::EDGE):
    nowned = data_.nedges_owned;
    nall = data_.nedges_all;
    break;
  case (Entity_kind::NODE):
    nowned = data_.nnodes_owned;
    nall = data_.nnodes_all;
    break;
  case (Entity_kind::BOUNDARY_FACE):
    nowned = data_.nboundary_faces_owned;
    nall = data_.nboundary_faces_all;
    break;
  case (Entity_kind::BOUNDARY_NODE):
    nowned = data_.nboundary_nodes_owned;
    nall = data_.nboundary_nodes_all;
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


Entity_ID
Mesh::getEntityParent(const Entity_kind kind, const Entity_ID entid) const
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


Cell_kind
Mesh::getCellKind(const Entity_ID c) const
{
  return Impl::computeCellKind(*this, c);
}


Parallel_kind
Mesh::getParallelKind(const Entity_kind& kind, const Entity_ID id) const
{
  if (id < getNumEntities(kind, Parallel_kind::OWNED)) {
    return Parallel_kind::OWNED;
  } else if (id < getNumEntities(kind, Parallel_kind::ALL)) {
    return Parallel_kind::GHOST;
  }
  return Parallel_kind::UNKNOWN;
}


//-----------------------
// Geometry: coordinates
//-----------------------
AmanziGeometry::Point
Mesh::getNodeCoordinate(const Entity_ID n) const
{
  return Impl::Getter<MEM>::get(
    data_.node_coordinates_cached,
    data_.node_coordinates,
    framework_mesh_,
    [&](const int i) { return framework_mesh_->getNodeCoordinate(i); },
    nullptr,
    n);
}


void
Mesh::setNodeCoordinate(const Entity_ID n, const AmanziGeometry::Point& p)
{
  if (data_.node_coordinates_cached) view<MemSpace_kind::HOST>(data_.node_coordinates)(n) = p;
  if (framework_mesh_.get()) framework_mesh_->setNodeCoordinate(n, p);
}


void
Mesh::setNodeCoordinates(const cEntity_ID_View& nodes, const cPoint_View& new_coords)
{
  auto bf = view<MemSpace_kind::HOST>(data_.node_coordinates);
  if (data_.node_coordinates_cached) {
    for (int i = 0; i != nodes.size(); ++i) { bf(nodes(i)) = new_coords(i); }
  }

  if (framework_mesh_.get()) {
    MeshFramework::cEntity_ID_View nodes_on_host;
    View_type<const AmanziGeometry::Point, MemSpace_kind::HOST> coords_on_host;

    for (int i = 0; i != nodes.size(); ++i) {
      framework_mesh_->setNodeCoordinate(nodes(i), new_coords(i));
    }
  }
  recacheGeometry();
}


//-----------------------
// Geometry: other
//-----------------------
bool
Mesh::isPointInCell(const AmanziGeometry::Point& p, const Entity_ID cellid) const
{
  if (data_.manifold_dim_ == 3) {
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

  } else if (data_.manifold_dim_ == 2) {
    auto ccoords = getCellCoordinates(cellid);
    return AmanziGeometry::point_in_polygon(p, ccoords);
  }

  return false;
}


//---------------------
// Downward adjacencies
//---------------------
//
// Cell --> Face
//
size_type
Mesh::getCellNumFaces(const Entity_ID c) const
{
  assert(data_.cell_faces_cached);
  return data_.cell_faces.size<MEM>(c);
}


const Entity_ID&
Mesh::getCellFace(const Entity_ID c, const size_type i) const
{
  assert(data_.cell_faces_cached);
  return data_.cell_faces.get<MEM>(c, i);
}


typename Mesh::cEntity_ID_View
Mesh::getCellFaces(const Entity_ID c) const
{
  return Impl::RaggedGetter<MEM>::get(
    data_.cell_faces_cached,
    data_.cell_faces,
    framework_mesh_,
    [&](const int i) {
      MeshFramework::cEntity_ID_View res;
      framework_mesh_->getCellFaces(c, res);
      return res;
    },
    nullptr,
    c);
}


void
Mesh::getCellFaces(const Entity_ID c, Mesh::cEntity_ID_View& faces) const
{
  faces = getCellFaces(c);
}


Kokkos::pair<Mesh::cEntity_ID_View,
             Mesh::cDirection_View>
Mesh::getCellFacesAndDirections(const Entity_ID c) const
{
  if (data_.cell_faces_cached) {
    return Kokkos::make_pair(data_.cell_faces.getRow<MEM>(c),
                             data_.cell_face_directions.getRow<MEM>(c));
  } else if (framework_mesh_.get()) {
    Kokkos::pair<Mesh::cEntity_ID_View,
                 Mesh::cDirection_View> ret;

    framework_mesh_->getCellFacesAndDirs(c, ret.first, &ret.second);
    return ret;
  }
  throwAccessError_("getCellFacesAndDirections");
  return Kokkos::pair<Mesh::cEntity_ID_View, Mesh::cDirection_View>();
}


void
Mesh::getCellFacesAndDirs(const Entity_ID c,
                          Mesh::cEntity_ID_View& faces,
                          Mesh::cDirection_View* const dirs) const
{
  auto [lfaces, ldirs] = getCellFacesAndDirections(c);
  faces = lfaces;
  if (dirs) (*dirs) = ldirs;
}


Kokkos::pair<Mesh::cEntity_ID_View,
             Mesh::cPoint_View>
Mesh::getCellFacesAndBisectors(const Entity_ID c) const
{
  if (data_.cell_faces_cached) {
    return Kokkos::make_pair(data_.cell_faces.getRow<MEM>(c),
                             data_.cell_face_bisectors.getRow<MEM>(c));
  } else if (algorithms_.get()) {
    Kokkos::pair<Mesh::cEntity_ID_View, Mesh::cPoint_View> ret;
    algorithms_->computeCellFacesAndBisectors(*this, c, ret.first, &ret.second);
    return ret;
  }
  throwAccessError_("getCellFacesAndBisectors");
  return Kokkos::pair<Mesh::cEntity_ID_View, Mesh::cPoint_View>();
}

void
Mesh::getCellFacesAndBisectors(const Entity_ID c,
        Mesh::cEntity_ID_View& faces,
        Mesh::cPoint_View* const bisectors) const
{
  auto [lfaces, lbisectors] = getCellFacesAndBisectors(c);
  faces = lfaces;
  if (bisectors) (*bisectors) = lbisectors;
}


//
// Cell --> Edge
//
size_type
Mesh::getCellNumEdges(const Entity_ID c) const
{
  assert(data_.cell_edges_cached);
  return data_.cell_edges.size<MEM>(c);
}


const Entity_ID&
Mesh::getCellEdge(const Entity_ID c, const size_type i) const
{
  assert(data_.cell_edges_cached);
  return data_.cell_edges.get<MEM>(c, i);
}


//
// Cell --> Node
//
size_type
Mesh::getCellNumNodes(const Entity_ID c) const
{
  assert(data_.cell_nodes_cached);
  return data_.cell_nodes.size<MEM>(c);
}


const Entity_ID&
Mesh::getCellNode(const Entity_ID c, const size_type i) const
{
  assert(data_.cell_nodes_cached);
  return data_.cell_nodes.get<MEM>(c, i);
}


//
// Face --> Edge
//
size_type
Mesh::getFaceNumEdges(const Entity_ID c) const
{
  assert(data_.face_edges_cached);
  return data_.face_edges.size<MEM>(c);
}


const Entity_ID&
Mesh::getFaceEdge(const Entity_ID f, const size_type i) const
{
  assert(data_.face_edges_cached);
  return data_.face_edges.get<MEM>(f, i);
}


typename Mesh::cEntity_ID_View
Mesh::getFaceEdges(const Entity_ID f) const
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


Kokkos::pair<Mesh::cEntity_ID_View,
             Mesh::cDirection_View>
Mesh::getFaceEdgesAndDirections(const Entity_ID f) const
{
  if (data_.face_edges_cached) {
    return Kokkos::make_pair(data_.face_edges.getRow<MEM>(f),
                             data_.face_edge_directions.getRow<MEM>(f));
  } else if (framework_mesh_.get()) {
    Kokkos::pair<Mesh::cEntity_ID_View,
                 Mesh::cDirection_View> ret;

    framework_mesh_->getFaceEdgesAndDirs(f, ret.first, &ret.second);
    return ret;
  }
  throwAccessError_("getFaceEdgesAndDirections");
  return Kokkos::pair<Mesh::cEntity_ID_View, Mesh::cDirection_View>();
}


void
Mesh::getFaceEdgesAndDirs(const Entity_ID f,
        Mesh::cEntity_ID_View& edges,
        Mesh::cDirection_View* const dirs) const
{
  auto [ledges, ldirs] = getFaceEdgesAndDirections(f);
  edges = ledges;
  if (dirs) *dirs = ldirs;
}


//
// Face --> Node
//
size_type
Mesh::getFaceNumNodes(const Entity_ID c) const
{
  assert(data_.face_nodes_cached);
  return data_.face_nodes.size<MEM>(c);
}


const Entity_ID&
Mesh::getFaceNode(const Entity_ID f, const size_type i) const
{
  assert(data_.face_nodes_cached);
  return data_.face_nodes.get<MEM>(f, i);
}


typename Mesh::cEntity_ID_View
Mesh::getFaceNodes(const Entity_ID f) const
{
  return Impl::RaggedGetter<MEM>::get(
    data_.face_nodes_cached,
    data_.face_nodes,
    framework_mesh_,
    [&](const int i) {
      MeshFramework::cEntity_ID_View res;
      framework_mesh_->getFaceNodes(f, res);
      return res;
    },
    nullptr,
    f);
}


//
// Edge --> Node
//
size_type
Mesh::getEdgeNumNodes(const Entity_ID e) const
{
  assert(data_.edge_nodes_cached);
  return data_.edge_nodes.size<MEM>(e);
}


const Entity_ID&
Mesh::getEdgeNode(const Entity_ID e, const size_type i) const
{
  assert(data_.edge_nodes_cached);
  return data_.edge_nodes.get<MEM>(e, i);
}


typename Mesh::cEntity_ID_View
Mesh::getEdgeNodes(const Entity_ID e) const
{
  return Impl::RaggedGetter<MEM>::get(
    data_.edge_nodes_cached,
    data_.edge_nodes,
    framework_mesh_,
    [&](const int i) {
      MeshFramework::cEntity_ID_View res;
      framework_mesh_->getEdgeNodes(e, res);
      return res;
    },
    nullptr,
    e);
}

void
Mesh::getEdgeNodes(const Entity_ID e, typename Mesh::cEntity_ID_View& nodes) const
{
  nodes = getEdgeNodes(e);
}


//
// Max counts
//
std::size_t
Mesh::getCellMaxNodes() const
{
  std::size_t n(0);
  Entity_ID ncells = getNumEntities(Entity_kind::CELL, Parallel_kind::OWNED);
  for (int c = 0; c < ncells; ++c) {
    auto nodes = getCellNodes(c);
    if (n < nodes.size()) n = nodes.size();
  }
  return n;
}


std::size_t
Mesh::getCellMaxFaces() const
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
Mesh::getCellMaxEdges() const
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


//-------------------
// Upward adjacencies
//-------------------
//
// Face --> Cell
//
size_type
Mesh::getFaceNumCells(const Entity_ID c) const
{
  assert(data_.face_cells_cached);
  return data_.face_cells.size<MEM>(c);
}


const Entity_ID&
Mesh::getFaceCell(const Entity_ID c, const size_type i) const
{
  assert(data_.face_cells_cached);
  return data_.face_cells.get<MEM>(c, i);
}


typename Mesh::cEntity_ID_View
Mesh::getFaceCells(const Entity_ID f) const
{
  return Impl::RaggedGetter<MEM>::get(
    data_.face_cells_cached,
    data_.face_cells,
    framework_mesh_,
    [&](const int i) {
      MeshFramework::cEntity_ID_View res;
      framework_mesh_->getFaceCells(f, res);
      return res;
    },
    nullptr,
    f);
}

void
Mesh::getFaceCells(const Entity_ID f, typename Mesh::cEntity_ID_View& cells) const
{
  cells = getFaceCells(f);
}


//
// Edge --> Cell
//
size_type
Mesh::getEdgeNumCells(const Entity_ID c) const
{
  assert(data_.edge_cells_cached);
  return data_.edge_cells.size<MEM>(c);
}


const Entity_ID&
Mesh::getEdgeCell(const Entity_ID c, const size_type i) const
{
  assert(data_.edge_cells_cached);
  return data_.edge_cells.get<MEM>(c, i);
}


//
// Edge --> Face
//
size_type
Mesh::getEdgeNumFaces(const Entity_ID c) const
{
  assert(data_.edge_faces_cached);
  return data_.edge_faces.size<MEM>(c);
}


const Entity_ID&
Mesh::getEdgeFace(const Entity_ID c, const size_type i) const
{
  assert(data_.edge_faces_cached);
  return data_.edge_faces.get<MEM>(c, i);
}


typename Mesh::cEntity_ID_View
Mesh::getEdgeFaces(const Entity_ID f) const
{
  return Impl::RaggedGetter<MEM>::get(
    data_.edge_faces_cached,
    data_.edge_faces,
    framework_mesh_,
    [&](const int i) {
      MeshFramework::cEntity_ID_View res;
      framework_mesh_->getEdgeFaces(f, res);
      return res;
    },
    nullptr,
    f);
}


//
// Node --> Face
//
typename Mesh::cEntity_ID_View
Mesh::getNodeFaces(const Entity_ID f) const
{
  return Impl::RaggedGetter<MEM>::get(
    data_.node_faces_cached,
    data_.node_faces,
    framework_mesh_,
    [&](const int i) {
      MeshFramework::cEntity_ID_View res;
      framework_mesh_->getNodeFaces(f, res);
      return res;
    },
    nullptr,
    f);
}


//
// Node --> Edge
//
typename Mesh::cEntity_ID_View
Mesh::getNodeEdges(const Entity_ID f) const
{
  return Impl::RaggedGetter<MEM>::get(
    data_.node_edges_cached,
    data_.node_edges,
    framework_mesh_,
    [&](const int i) {
      MeshFramework::cEntity_ID_View res;
      framework_mesh_->getNodeEdges(f, res);
      return res;
    },
    nullptr,
    f);
}


// --------------------------------------
// All purpose error message
// --------------------------------------
void
Mesh::throwAccessError_(const std::string& func_name) const
{
  Errors::Message msg;
  msg << "Unable to determine valid strategy to find mesh quantity in \"" << func_name << "\"";
  Exceptions::amanzi_throw(msg);
}

} // namespace AmanziMesh
} // namespace Amanzi
