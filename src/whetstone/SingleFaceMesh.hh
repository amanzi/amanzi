/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Single face mesh has one less topological dimension than the 
  base mesh used in the constructor. Like in the base mesh, edge
  and face are equivalent.
*/

#ifndef AMANZI_WHETSTONE_SINGLE_FACE_MESH_HH_
#define AMANZI_WHETSTONE_SINGLE_FACE_MESH_HH_

#include "Point.hh"

#include "MeshLight.hh"
#include "SurfaceCoordinateSystem.hh"

namespace Amanzi {
namespace WhetStone {

class SingleFaceMesh : public AmanziMesh::MeshLight {
 public:
  SingleFaceMesh(const Teuchos::RCP<const AmanziMesh::MeshLight>& mesh, int f,
                 const SurfaceCoordinateSystem& coordsys) {
    BuildCache_(mesh, f, coordsys); 
  }

  SingleFaceMesh(const Teuchos::RCP<const AmanziMesh::MeshLight>& mesh, int f) {
    const auto& xf = mesh->face_centroid(f);
    const auto& normal = mesh->face_normal(f);
    SurfaceCoordinateSystem coordsys(xf, normal);
    BuildCache_(mesh, f, coordsys);
   }

  virtual ~SingleFaceMesh() {};

  // ---------------------
  // Downward connectivity
  // ---------------------
  virtual void cell_get_nodes(const AmanziMesh::Entity_ID c,
                              AmanziMesh::Entity_ID_List *nodes) const override {
    AMANZI_ASSERT(c == 0);
    *nodes = cell_node_ids_;
  }

  virtual void face_get_nodes(const AmanziMesh::Entity_ID f,
                              AmanziMesh::Entity_ID_List *nodes) const override {
    AMANZI_ASSERT(f < nnodes_);
    *nodes = face_node_ids_[f];
  }

  virtual void edge_get_nodes(const AmanziMesh::Entity_ID e,
                              AmanziMesh::Entity_ID* n0, AmanziMesh::Entity_ID* n1) const override {
    AMANZI_ASSERT(e < nnodes_);
    *n0 = face_node_ids_[e][0]; 
    *n1 = face_node_ids_[e][1]; 
  }

  // -------------------
  // Upward connectivity
  // -------------------
  virtual void node_get_cells(
          const AmanziMesh::Entity_ID v,
          const AmanziMesh::Parallel_type ptype,
          AmanziMesh::Entity_ID_List *cells) const override { AMANZI_ASSERT(v < nnodes_); cells->resize(1, 0); }

  virtual void node_get_faces(
          const AmanziMesh::Entity_ID v,
          const AmanziMesh::Parallel_type ptype,
          AmanziMesh::Entity_ID_List *faces) const override {
    AMANZI_ASSERT(v < nnodes_);
    faces->resize(2);
    (*faces)[0] = v; 
    (*faces)[1] = (v + 1) % nnodes_;
  }

  // --------
  // Geometry
  // --------
  virtual void node_get_coordinates(
          const AmanziMesh::Entity_ID v, AmanziGeometry::Point* xp) const override {
    AMANZI_ASSERT(v < nnodes_); *xp = cell_coords_[v];
  }

  virtual void face_get_coordinates(
          const AmanziMesh::Entity_ID f,
          std::vector<AmanziGeometry::Point> *fcoords) const override {
    AMANZI_ASSERT(f < nnodes_);
    fcoords->resize(2);
    const auto& nodes = face_node_ids_[f];
    (*fcoords)[0] = cell_coords_[nodes[0]];
    (*fcoords)[1] = cell_coords_[nodes[1]];
  }

  virtual void cell_get_coordinates(
          const AmanziMesh::Entity_ID c,
          std::vector<AmanziGeometry::Point> *ccoords) const override { *ccoords = cell_coords_; }

  virtual unsigned int num_entities(
          const AmanziMesh::Entity_kind kind,
          const AmanziMesh::Parallel_type ptype) const override {
    return (kind == AmanziMesh::CELL) ? 1 : nnodes_;
  }

  virtual AmanziMesh::Parallel_type entity_get_ptype(
          const AmanziMesh::Entity_kind kind, const AmanziMesh::Entity_ID ent) const override {
    return AmanziMesh::Parallel_type::OWNED;
  }

 protected:
  // is not used and will be removed when light mesh becomes cache-only 
  virtual void cell_get_faces_and_dirs_internal_(
          const AmanziMesh::Entity_ID c,
          AmanziMesh::Entity_ID_List* faces,
          std::vector<int>* fdirs,
          const bool ordered = false) const override {};

  virtual void cell_get_edges_internal_(
          const AmanziMesh::Entity_ID c,
          AmanziMesh::Entity_ID_List* edges) const override {};

  virtual void face_get_edges_and_dirs_internal_(
          const AmanziMesh::Entity_ID f,
          AmanziMesh::Entity_ID_List* edges,
          std::vector<int>* edirs,
          const bool ordered = true) const override {};

  virtual void face_get_cells_internal_(
          const AmanziMesh::Entity_ID f,
          const AmanziMesh::Parallel_type ptype,
          AmanziMesh::Entity_ID_List* cells) const override {};

  // geometries
  virtual int compute_cell_geometry_(
          const AmanziMesh::Entity_ID c,
          double* volume,
          AmanziGeometry::Point* xc) const override { return 0; }

  virtual int compute_face_geometry_(
          const AmanziMesh::Entity_ID f,
          double* area,
          AmanziGeometry::Point* xf,
          std::vector<AmanziGeometry::Point> *normals) const override { return 0; }

  virtual int compute_edge_geometry_(
          const AmanziMesh::Entity_ID e,
          double *length,
          AmanziGeometry::Point* tau) const override { return 0; }

 private:
  void BuildCache_(
    const Teuchos::RCP<const AmanziMesh::MeshLight>& mesh, int f,
    const SurfaceCoordinateSystem& coordsys);

 private:
  int nnodes_;
  AmanziMesh::Entity_ID_List cell_node_ids_;
  std::vector<AmanziGeometry::Point> cell_coords_;

  std::vector<AmanziMesh::Entity_ID_List> face_node_ids_;
};


/* ******************************************************************
* Construct cache data
****************************************************************** */
inline
void SingleFaceMesh::BuildCache_(
    const Teuchos::RCP<const AmanziMesh::MeshLight>& mesh, int f,
    const SurfaceCoordinateSystem& coordsys)
{
  int d = mesh->space_dimension() - 1;
  set_space_dimension(d);
  set_manifold_dimension(d);

  AmanziMesh::Entity_ID_List fedges, fnodes, enodes, cells;
  mesh->face_get_nodes(f, &fnodes);
  nnodes_ = fnodes.size();

  std::vector<int> fdirs(nnodes_, 1);
  if (mesh->valid_edges()) {
    mesh->face_get_edges_and_dirs(f, &fedges, &fdirs);
  }

  // single surface cell
  fedges.resize(nnodes_);
  for (int i = 0; i < nnodes_; ++i) fedges[i] = i;

  cell_face_ids_.resize(1);
  cell_face_ids_[0] = fedges;

  cell_face_dirs_.resize(1);
  cell_face_dirs_[0] = fdirs;

  AmanziGeometry::Point xyz(d);
  const auto& xf = mesh->face_centroid(f);
  cell_centroids_.resize(1, coordsys.Project(xf, true));

  double volume = mesh->face_area(f);
  cell_volumes_.resize(1, volume);

  // cell nodes
  cell_node_ids_.resize(nnodes_);
  cell_coords_.resize(nnodes_);

  for (int i = 0; i < nnodes_; ++i) {
    cell_node_ids_[i] = i;
    mesh->node_get_coordinates(fnodes[i], &xyz);
    cell_coords_[i] = coordsys.Project(xyz, true);
  }

  // cell faces
  enodes.resize(2);
  face_node_ids_.resize(nnodes_);
  face_areas_.resize(nnodes_);
  face_centroids_.resize(nnodes_);
  face_normals_.resize(nnodes_);

  for (int i = 0; i < nnodes_; ++i) {
    int n0(i), n1((i + 1) % nnodes_);

    enodes[0] = (fdirs[i] > 0) ? n0 : n1;
    enodes[1] = (fdirs[i] > 0) ? n1 : n0;
    face_node_ids_[i] = enodes;

    auto tau = (cell_coords_[n1] - cell_coords_[n0]) * fdirs[i];
    face_areas_[i] = norm(tau);

    face_normals_[i].resize(1, AmanziGeometry::Point(tau[1], -tau[0]));
    face_centroids_[i] = (cell_coords_[n1] + cell_coords_[n0]) / 2;
  }

  // cell edges
  edge_vectors_.resize(nnodes_);
  edge_lengths_ = face_areas_;

  for (int i = 0; i < nnodes_; ++i) {
    int n0(i), n1((i + 1) % nnodes_);
    edge_vectors_[i] = cell_coords_[n1] - cell_coords_[n0];
  }

  // backward connectivity
  cells.resize(1, 0);
  std::vector<int> edirs(1, 1);
  face_cell_ids_.resize(nnodes_, cells);

  face_edge_ids_.resize(nnodes_);
  face_edge_dirs_.resize(nnodes_, edirs);

  for (int i = 0; i < nnodes_; ++i) {
    cells.resize(1, i);
    face_edge_ids_[i] = cells;
  }

  // cache is done (this probably is not needed)
  faces_requested_ = true;
  edges_requested_ = true;

  cell2face_info_cached_ = true;
  cell2edge_info_cached_ = true;
  face2cell_info_cached_ = true;
  face2edge_info_cached_ = true;

  cell_geometry_precomputed_ = true;
  face_geometry_precomputed_ = true;
  edge_geometry_precomputed_ = true;
}

} // namespace WhetStone
} // namespace Amanzi

#endif
