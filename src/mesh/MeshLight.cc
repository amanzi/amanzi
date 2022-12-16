/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/*
  Mesh

  Konstantin Lipnikov
  Rao Garimella

  The smallest mesh class with basic connectivity and geometry
*/

#include "MeshLight.hh"

#define AMANZI_MESH_CACHE_VARS 1

namespace Amanzi {
namespace AmanziMesh {

// -------------------------------------------------------------------
// Cache data
// -------------------------------------------------------------------
void
MeshLight::BuildCache()
{
  if (!cell2face_info_cached_ && faces_requested_) cache_cell2face_info_();
  if (!cell2edge_info_cached_ && edges_requested_) cache_cell2edge_info_();
}


// -------------------------------------------------------------------
// Downward connectivity: f -> e
// -------------------------------------------------------------------
void
MeshLight::face_get_edges_and_dirs(const Entity_ID faceid,
                                   Entity_ID_List* edgeids,
                                   std::vector<int>* edge_dirs,
                                   bool ordered) const
{
#if AMANZI_MESH_CACHE_VARS != 0
  if (!face2edge_info_cached_) cache_face2edge_info_();

  *edgeids = face_edge_ids_[faceid]; // copy operation

  if (edge_dirs) {
    std::vector<int>& fedgedirs = face_edge_dirs_[faceid];
    *edge_dirs = fedgedirs; // copy operation
  }

#else // Non-cached version
  face_get_edges_and_dirs_internal_(faceid, edgeids, edge_dirs, ordered);
#endif
}


// -------------------------------------------------------------------
// Upward connectivity: f -> c
// -------------------------------------------------------------------
void
MeshLight::face_to_cell_edge_map(const Entity_ID faceid,
                                 const Entity_ID cellid,
                                 std::vector<int>* map) const
{
#if AMANZI_MESH_CACHE_VARS != 0
  if (!face2edge_info_cached_) cache_face2edge_info_();
  if (!cell2edge_info_cached_) cache_cell2edge_info_();

  map->resize(face_edge_ids_[faceid].size());
  for (int f = 0; f < face_edge_ids_[faceid].size(); ++f) {
    Entity_ID fedge = face_edge_ids_[faceid][f];

    for (int c = 0; c < cell_edge_ids_[cellid].size(); ++c) {
      if (fedge == cell_edge_ids_[cellid][c]) {
        (*map)[f] = c;
        break;
      }
    }
  }

#else // non-cached version
  Entity_ID_List fedgeids, cedgeids;
  std::vector<int> fedgedirs;

  face_get_edges_and_dirs(faceid, &fedgeids, &fedgedirs, true);
  cell_get_edges(cellid, &cedgeids);

  map->resize(fedgeids.size(), -1);
  for (int f = 0; f < fedgeids.size(); ++f) {
    Entity_ID fedge = fedgeids[f];

    for (int c = 0; c < cedgeids.size(); ++c) {
      if (fedge == cedgeids[c]) {
        (*map)[f] = c;
        break;
      }
    }
  }
#endif
}


// -------------------------------------------------------------------
// Centroids: cells
// -------------------------------------------------------------------
AmanziGeometry::Point
MeshLight::cell_centroid(const Entity_ID cellid, bool recompute) const
{
  if (!cell_geometry_precomputed_) {
    compute_cell_geometric_quantities_();
    return cell_centroids_[cellid];
  } else {
    if (recompute) {
      double volume;
      AmanziGeometry::Point centroid(space_dim_);
      compute_cell_geometry_(cellid, &volume, &centroid);
      return centroid;
    } else
      return cell_centroids_[cellid];
  }
}


// -------------------------------------------------------------------
// Cache: c -> f
// -------------------------------------------------------------------
void
MeshLight::cache_cell2face_info_() const
{
  int ncells = num_entities(CELL, Parallel_type::ALL);
  cell_face_ids_.resize(ncells);
  cell_face_dirs_.resize(ncells);

  int nfaces = num_entities(FACE, Parallel_type::ALL);
  face_cell_ids_.resize(nfaces);
  face_cell_ptype_.resize(nfaces);

  for (int c = 0; c < ncells; c++) {
    cell_get_faces_and_dirs_internal_(c, &(cell_face_ids_[c]), &(cell_face_dirs_[c]), false);

    int nf = cell_face_ids_[c].size();
    for (int jf = 0; jf < nf; jf++) {
      Entity_ID f = cell_face_ids_[c][jf];
      int dir = cell_face_dirs_[c][jf];
      face_cell_ids_[f].push_back(dir > 0 ? c : ~c); // store 1s complement of c if dir is -ve
      face_cell_ptype_[f].push_back(entity_get_ptype(CELL, c));
    }
  }

  cell2face_info_cached_ = true;
  face2cell_info_cached_ = true;
  faces_requested_ = true;
}


// -------------------------------------------------------------------
// Cache: c -> e
// -------------------------------------------------------------------
void
MeshLight::cache_cell2edge_info_() const
{
  int ncells = num_entities(CELL, Parallel_type::ALL);
  cell_edge_ids_.resize(ncells);

  if (space_dim_ == 2) {
    std::vector<int> dirs;
    for (int c = 0; c < ncells; c++)
      cell_get_faces_and_dirs_internal_(c, &(cell_edge_ids_[c]), &dirs);
  } else
    for (int c = 0; c < ncells; c++) cell_get_edges_internal_(c, &(cell_edge_ids_[c]));

  cell2edge_info_cached_ = true;
}


// -------------------------------------------------------------------
// Cache: f -> c
// The results are cached the first time it is called and then return
// the cached results subsequently.
// -------------------------------------------------------------------
void
MeshLight::face_get_cells(const Entity_ID faceid,
                          const Parallel_type ptype,
                          Entity_ID_List* cellids) const
{
#if AMANZI_MESH_CACHE_VARS != 0
  if (!face2cell_info_cached_) cache_cell2face_info_();

  cellids->reserve(2);
  cellids->clear();
  int n = face_cell_ptype_[faceid].size();

  for (int i = 0; i < n; i++) {
    Parallel_type cell_ptype = face_cell_ptype_[faceid][i];
    if (cell_ptype == Parallel_type::PTYPE_UNKNOWN) continue;

    int c = face_cell_ids_[faceid][i];
    if (std::signbit(c)) c = ~c; // strip dir info by taking 1s complement

    if (ptype == Parallel_type::ALL || ptype == cell_ptype) cellids->emplace_back(c);
  }

#else // Non-cached version
  Entity_ID_List fcells;
  face_get_cells_internal_(faceid, Parallel_type::ALL, &fcells);

  cellids->clear();
  int n = face_cell_ptype_[faceid].size();

  switch (ptype) {
  case Parallel_type::ALL:
    for (int i = 0; i < n; i++)
      if (entity_get_ptype(CELL, fcells[i]) != PTYPE_UNKNOWN) cellids->push_back(fcells[i]);
    break;
  case Parallel_type::OWNED:
    for (int i = 0; i < n; i++)
      if (entity_get_ptype(CELL, fcells[i]) == Parallel_type::OWNED) cellids->push_back(fcells[i]);
    break;
  case Parallel_type::GHOST:
    for (int i = 0; i < n; i++)
      if (entity_get_ptype(CELL, fcells[i]) == Parallel_type::GHOST) cellids->push_back(fcells[i]);
    break;
  default:
    break;
  }
#endif
}


// -------------------------------------------------------------------
// Cache: f -> e
// -------------------------------------------------------------------
void
MeshLight::cache_face2edge_info_() const
{
  int nfaces = num_entities(FACE, Parallel_type::ALL);
  face_edge_ids_.resize(nfaces);
  face_edge_dirs_.resize(nfaces);

  for (int f = 0; f < nfaces; f++) {
    Entity_ID_List fedgeids;
    std::vector<int> fedgedirs;

    face_get_edges_and_dirs_internal_(f, &(face_edge_ids_[f]), &(face_edge_dirs_[f]), true);
  }

  face2edge_info_cached_ = true;
  faces_requested_ = true;
  edges_requested_ = true;
}


// -------------------------------------------------------------------
// # (c -> f)
// -------------------------------------------------------------------
unsigned int
MeshLight::cell_get_num_faces(const Entity_ID cellid) const
{
#if AMANZI_MESH_CACHE_VARS != 0
  if (!cell2face_info_cached_) cache_cell2face_info_();
  return cell_face_ids_[cellid].size();

#else // Non-cached version
  Entity_ID_List cfaceids;
  std::vector<int> cfacedirs;

  cell_get_faces_and_dirs_internal_(cellid, &cfaceids, &cfacedirs, false);
  return cfaceids.size();
#endif
}


// -------------------------------------------------------------------
// cache: cell geometry
// -------------------------------------------------------------------
int
MeshLight::compute_cell_geometric_quantities_() const
{
  int ncells = num_entities(CELL, Parallel_type::ALL);

  cell_volumes_.resize(ncells);
  cell_centroids_.resize(ncells);
  for (int i = 0; i < ncells; i++) {
    double volume;
    AmanziGeometry::Point centroid(space_dim_);

    compute_cell_geometry_(i, &volume, &centroid);

    cell_volumes_[i] = volume;
    cell_centroids_[i] = centroid;
  }

  cell_geometry_precomputed_ = true;

  return 1;
}


// -------------------------------------------------------------------
// cache: cell geometry
// -------------------------------------------------------------------
double
MeshLight::cell_volume(const Entity_ID c, const bool recompute) const
{
  if (!cell_geometry_precomputed_) {
    compute_cell_geometric_quantities_();
    return cell_volumes_[c];
  } else {
    if (recompute) {
      double volume;
      AmanziGeometry::Point centroid(space_dim_);
      compute_cell_geometry_(c, &volume, &centroid);
      return volume;
    } else {
      return cell_volumes_[c];
    }
  }
}


// -------------------------------------------------------------------
// cache: face geometry
// The face centroid is computed as the area weighted average of the
// centroids of the triangles from a symmetric triangular
// decomposition of the face. Each triangular facet is formed by the
// connecting the face center (average of face nodes) to the two
// nodes of an edge of the face
// -------------------------------------------------------------------
AmanziGeometry::Point
MeshLight::face_centroid(const Entity_ID f, const bool recompute) const
{
  AMANZI_ASSERT(faces_requested_);

  if (!face_geometry_precomputed_) {
    compute_face_geometric_quantities_();
    return face_centroids_[f];
  } else {
    if (recompute) {
      double area;
      AmanziGeometry::Point centroid(space_dim_);
      std::vector<AmanziGeometry::Point> normals;
      compute_face_geometry_(f, &area, &centroid, &normals);
      return centroid;
    } else {
      return face_centroids_[f];
    }
  }
}


// -------------------------------------------------------------------
// cache: face normal
// -------------------------------------------------------------------
AmanziGeometry::Point
MeshLight::face_normal(const Entity_ID faceid,
                       const bool recompute,
                       const Entity_ID cellid,
                       int* orientation) const
{
  AMANZI_ASSERT(faces_requested_);

  std::vector<AmanziGeometry::Point>* fnormals = nullptr;
  std::vector<AmanziGeometry::Point> fnormals_new;

  if (!face_geometry_precomputed_) {
    compute_face_geometric_quantities_();
    fnormals = &(face_normals_[faceid]);
  } else {
    if (recompute) {
      double area;
      AmanziGeometry::Point centroid(space_dim_);
      compute_face_geometry_(faceid, &area, &centroid, &fnormals_new);
      fnormals = &fnormals_new;
    } else {
      fnormals = &(face_normals_[faceid]);
    }
  }

  AMANZI_ASSERT(fnormals->size() > 0);

  if (cellid == -1) {
    // Return the natural normal. This is the normal with respect to
    // the first cell, appropriately adjusted according to whether the
    // face is pointing into the cell (-ve cell id) or out

    int c = face_cell_ids_[faceid][0];
    return std::signbit(c) ? -(*fnormals)[0] : (*fnormals)[0];
  } else {
    // Find the index of 'cellid' in list of cells connected to face

    int dir;
    int irefcell;
    int nfc = face_cell_ids_[faceid].size();
    for (irefcell = 0; irefcell < nfc; irefcell++) {
      int c = face_cell_ids_[faceid][irefcell];
      if (c == cellid || ~c == cellid) {
        dir = std::signbit(c) ? -1 : 1;
        break;
      }
    }
    AMANZI_ASSERT(irefcell < nfc);
    if (orientation) *orientation = dir; // if orientation was requested
    return (*fnormals)[irefcell];
  }
}


// -------------------------------------------------------------------
// cache: face area
// -------------------------------------------------------------------
double
MeshLight::face_area(const Entity_ID faceid, const bool recompute) const
{
  AMANZI_ASSERT(faces_requested_);

  if (!face_geometry_precomputed_) {
    compute_face_geometric_quantities_();
    return face_areas_[faceid];
  } else {
    if (recompute) {
      double area;
      AmanziGeometry::Point centroid(space_dim_);
      std::vector<AmanziGeometry::Point> normals;
      compute_face_geometry_(faceid, &area, &centroid, &normals);
      return area;
    } else {
      return face_areas_[faceid];
    }
  }
}


// -------------------------------------------------------------------
// edge geometry: centroid
// -------------------------------------------------------------------
AmanziGeometry::Point
MeshLight::edge_centroid(const Entity_ID e) const
{
  Entity_ID p0, p1;
  AmanziGeometry::Point xyz0, xyz1;

  edge_get_nodes(e, &p0, &p1);
  node_get_coordinates(p0, &xyz0);
  node_get_coordinates(p1, &xyz1);
  return (xyz0 + xyz1) / 2;
}


// -------------------------------------------------------------------
// edge geometry: vector
// -------------------------------------------------------------------
AmanziGeometry::Point
MeshLight::edge_vector(const Entity_ID edgeid,
                       const bool recompute,
                       const Entity_ID pointid,
                       int* orientation) const
{
  AMANZI_ASSERT(edges_requested_);

  AmanziGeometry::Point evector(space_dim_);
  AmanziGeometry::Point& evector_ref = evector; // to avoid extra copying

  if (!edge_geometry_precomputed_) compute_edge_geometric_quantities_();

  if (recompute) {
    double length;
    compute_edge_geometry_(edgeid, &length, &evector);
    // evector_ref already points to evector
  } else {
    evector_ref = edge_vectors_[edgeid];
  }

  if (orientation) *orientation = 1;

  if (pointid == -1)
    return evector_ref;
  else {
    Entity_ID p0, p1;
    edge_get_nodes(edgeid, &p0, &p1);

    if (pointid == p0)
      return evector_ref;
    else {
      if (orientation) *orientation = -1;
      return -evector_ref;
    }
  }
}


// -------------------------------------------------------------------
// edge geometry: length
// -------------------------------------------------------------------
double
MeshLight::edge_length(const Entity_ID e, const bool recompute) const
{
  AMANZI_ASSERT(edges_requested_);

  if (!edge_geometry_precomputed_) {
    compute_edge_geometric_quantities_();
    return edge_lengths_[e];
  } else {
    if (recompute) {
      double length;
      AmanziGeometry::Point vector(space_dim_);
      compute_edge_geometry_(e, &length, &vector);
      return length;
    } else {
      return edge_lengths_[e];
    }
  }
}


// -------------------------------------------------------------------
// cache: face geometry
// -------------------------------------------------------------------
int
MeshLight::compute_face_geometric_quantities_() const
{
  if (space_dimension() == 3 && manifold_dimension() == 2) {
    // need cell centroids to compute normals
    if (!cell_geometry_precomputed_) compute_cell_geometric_quantities_();
  }

  int nfaces = num_entities(FACE, Parallel_type::ALL);

  face_areas_.resize(nfaces);
  face_centroids_.resize(nfaces);
  face_normals_.resize(nfaces);

  for (int i = 0; i < nfaces; i++) {
    double area;
    AmanziGeometry::Point centroid(space_dim_);
    std::vector<AmanziGeometry::Point> normals;

    // normal0 and normal1 are outward normals of the face with
    // respect to the cell0 and cell1 of the face. The natural normal
    // of the face points out of cell0 and into cell1. If one of these
    // cells do not exist, then the normal is the null vector.
    compute_face_geometry_(i, &area, &centroid, &normals);

    face_areas_[i] = area;
    face_centroids_[i] = centroid;
    face_normals_[i] = normals;
  }

  face_geometry_precomputed_ = true;

  return 1;
}


// -------------------------------------------------------------------
// cache: edge geometry
// -------------------------------------------------------------------
int
MeshLight::compute_edge_geometric_quantities_() const
{
  int nedges = num_entities(EDGE, Parallel_type::ALL);

  edge_vectors_.resize(nedges);
  edge_lengths_.resize(nedges);

  for (int i = 0; i < nedges; i++) {
    double length;
    AmanziGeometry::Point evector(space_dim_);

    compute_edge_geometry_(i, &length, &evector);

    edge_lengths_[i] = length;
    edge_vectors_[i] = evector;
  }

  edge_geometry_precomputed_ = true;

  return 1;
}


// -------------------------------------------------------------------
// NOTE: this is a on-processor routine
// -------------------------------------------------------------------
unsigned int
MeshLight::cell_get_max_faces() const
{
  unsigned int n(0);
  int ncells = num_entities(CELL, Parallel_type::OWNED);
  for (int c = 0; c < ncells; ++c) { n = std::max(n, cell_get_num_faces(c)); }
  return n;
}


// -------------------------------------------------------------------
// NOTE: this is a on-processor routine
// -------------------------------------------------------------------
unsigned int
MeshLight::cell_get_max_edges() const
{
  unsigned int n(0);
  if (edges_requested_) {
    int ncells = num_entities(CELL, Parallel_type::OWNED);
    for (int c = 0; c < ncells; ++c) {
      const auto& edges = cell_get_edges(c);
      n = std::max(n, (unsigned int)edges.size());
    }
  }
  return n;
}


// -------------------------------------------------------------------
// NOTE: this is a on-processor routine
// -------------------------------------------------------------------
unsigned int
MeshLight::cell_get_max_nodes() const
{
  unsigned int n(0);
  int ncells = num_entities(CELL, Parallel_type::OWNED);
  for (int c = 0; c < ncells; ++c) {
    AmanziMesh::Entity_ID_List nodes;
    cell_get_nodes(c, &nodes);
    n = std::max(n, (unsigned int)nodes.size());
  }
  return n;
}

} // namespace AmanziMesh
} // namespace Amanzi
