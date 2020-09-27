/*
  Mesh

  Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Konstantin Lipnikov
  Rao Garimella

  The smallest mesh class with basic connectivity and geometry
*/

#include "MeshLight.hh"

#define AMANZI_MESH_CACHE_VARS 1

namespace Amanzi {
namespace AmanziMesh {

// -------------------------------------------------------------------
// Downward connectivity: c -> f
// -------------------------------------------------------------------
void MeshLight::cell_get_faces_and_dirs(
    const Entity_ID cellid,
    Entity_ID_List *faceids,
    std::vector<int> *face_dirs,
    const bool ordered) const
{
#if AMANZI_MESH_CACHE_VARS != 0
  if (!cell2face_info_cached_) cache_cell_face_info_();

  if (ordered)
    cell_get_faces_and_dirs_internal_(cellid, faceids, face_dirs, ordered);
  else {
    Entity_ID_List &cfaceids = cell_face_ids_[cellid];
    *faceids = cfaceids; // copy operation

    if (face_dirs) {
      std::vector<int> &cfacedirs = cell_face_dirs_[cellid];
      *face_dirs = cfacedirs; // copy operation
    }
  }
#else // Non-cached version
  cell_get_faces_and_dirs_internal_(cellid, faceids, face_dirs, ordered);
#endif
}


// -------------------------------------------------------------------
// Downward connectivity: c -> e
// -------------------------------------------------------------------
void MeshLight::cell_get_edges(const Entity_ID cellid,
                               Entity_ID_List *edgeids) const
{
#if AMANZI_MESH_CACHE_VARS != 0
  if (!cell2edge_info_cached_) cache_cell2edge_info_();

  *edgeids = cell_edge_ids_[cellid];  // copy operation

#else
  cell_get_edges_internal_(cellid, edgeids);
#endif
}


// -------------------------------------------------------------------
// Cache: c -> f
// -------------------------------------------------------------------
void MeshLight::cache_cell_face_info_() const
{
  int ncells = num_entities(CELL, Parallel_type::ALL);
  cell_face_ids_.resize(ncells);
  cell_face_dirs_.resize(ncells);

  int nfaces = num_entities(FACE, Parallel_type::ALL);
  face_cell_ids_.resize(nfaces);
  face_cell_ptype_.resize(nfaces);

  for (int c = 0; c < ncells; c++) {
    cell_get_faces_and_dirs_internal_(c, &(cell_face_ids_[c]),
                                      &(cell_face_dirs_[c]), false);

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
void MeshLight::cache_cell2edge_info_() const
{
  int ncells = num_entities(CELL,Parallel_type::ALL);
  cell_edge_ids_.resize(ncells);

  if (space_dim_ == 2) {
    cell_2D_edge_dirs_.resize(ncells);
    for (int c = 0; c < ncells; c++)
      cell_2D_get_edges_and_dirs_internal_(c, &(cell_edge_ids_[c]),
              &(cell_2D_edge_dirs_[c]));
  }
  else
    for (int c = 0; c < ncells; c++)
      cell_get_edges_internal_(c, &(cell_edge_ids_[c]));

  cell2edge_info_cached_ = true;
}


// -------------------------------------------------------------------
// # (c -> f)
// -------------------------------------------------------------------
unsigned int MeshLight::cell_get_num_faces(const Entity_ID cellid) const
{
#if AMANZI_MESH_CACHE_VARS != 0
  if (!cell2face_info_cached_) cache_cell_face_info_();
  return cell_face_ids_[cellid].size();

#else  // Non-cached version
  Entity_ID_List cfaceids;
  std::vector<int> cfacedirs;

  cell_get_faces_and_dirs_internal_(cellid, &cfaceids, &cfacedirs, false);
  return cfaceids.size();
#endif
}


// -------------------------------------------------------------------
// NOTE: this is a on-processor routine
// -------------------------------------------------------------------
unsigned int MeshLight::cell_get_max_faces() const
{
  unsigned int n(0);
  int ncells = num_entities(CELL, Parallel_type::OWNED);
  for (int c = 0; c < ncells; ++c) {
    n = std::max(n, cell_get_num_faces(c));
  }
  return n;
}


// -------------------------------------------------------------------
// NOTE: this is a on-processor routine
// -------------------------------------------------------------------
unsigned int MeshLight::cell_get_max_edges() const
{
  unsigned int n(0);
  if (edges_requested_) {
    int ncells = num_entities(CELL, Parallel_type::OWNED);
    for (int c = 0; c < ncells; ++c) {
      AmanziMesh::Entity_ID_List edges;
      cell_get_edges(c, &edges);
      n = std::max(n, (unsigned int) edges.size());
    }
  }
  return n;
}


// -------------------------------------------------------------------
// NOTE: this is a on-processor routine
// -------------------------------------------------------------------
unsigned int MeshLight::cell_get_max_nodes() const
{
  unsigned int n(0);
  int ncells = num_entities(CELL, Parallel_type::OWNED);
  for (int c = 0; c < ncells; ++c) {
    AmanziMesh::Entity_ID_List nodes;
    cell_get_nodes(c, &nodes);
    n = std::max(n, (unsigned int) nodes.size());
  }
  return n;
}

}  // namespace AmanziMesh
}  // namespace Amanzi

