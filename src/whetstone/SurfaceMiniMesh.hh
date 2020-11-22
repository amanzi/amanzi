/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Surface mini mesh has one less topological dimension than the 
  base mesh used in the constructor.
*/

#ifndef AMANZI_WHETSTONE_SURFACE_MINI_MESH_HH_
#define AMANZI_WHETSTONE_SURFACE_MINI_MESH_HH_

#include "Point.hh"

#include "SurfaceCoordinateSystem.hh"

namespace Amanzi {
namespace WhetStone {

class SurfaceMiniMesh {
 public:
  SurfaceMiniMesh(const Teuchos::RCP<const AmanziMesh::MeshLight>& mesh,
                  std::shared_ptr<const SurfaceCoordinateSystem> coordsys)
    : mesh_(mesh),
      coordsys_(coordsys),
      d_(mesh->space_dimension() - 1) {};
  ~SurfaceMiniMesh() {};

  // base class functionality
  int space_dimension() const { return d_; }

  // -- downward connectivity
  void cell_get_faces(Entity_ID c, Entity_ID_List *faces) const {
    std::vector<int> dirs;
    mesh_->face_get_edges_and_dirs(c, faces, &dirs);
  }

  void cell_get_faces_and_dirs(Entity_ID c, Entity_ID_List *faces, std::vector<int> *dirs) const {
    mesh_->face_get_edges_and_dirs(c, faces, dirs);
  }

  void cell_get_edges(Entity_ID c, Entity_ID_List *edges) const {
    std::vector<int> dirs;
    mesh_->face_get_edges_and_dirs(c, edges, &dirs);
  }

  void cell_get_nodes(Entity_ID c, Entity_ID_List *nodes) const {
    mesh_->face_get_nodes(c, nodes);
  }

  void face_get_edges_and_dirs(const Entity_ID f, Entity_ID_List *edges, std::vector<int> *dirs) const {
    Entity_ID v0, v1;
    mesh_->edge_get_nodes(f, &v0, &v1);
    edges->clear();
    edges->push_back(v0);
    edges->push_back(v1);
    dirs->resize(2, 1);
  }

  void face_get_nodes(Entity_ID f, Entity_ID_List *nodes) const {
    Entity_ID v0, v1;
    mesh_->edge_get_nodes(f, &v0, &v1);
    nodes->clear();
    nodes->push_back(v0);
    nodes->push_back(v1);
  }

  void edge_get_nodes(const Entity_ID e, Entity_ID* n0, Entity_ID* n1) const {
    *n1 = *n0 = e;
  }

  // -- geometric objects
  AmanziGeometry::Point cell_centroid(Entity_ID c) const {
    const auto& xf = mesh_->face_centroid(c);
    return coordsys_->Project(xf, true);
  }

  double cell_volume(Entity_ID c) const {
    return mesh_->face_area(c);
  }

  // calculus of vectors is the same in 2D and 3D; however, since
  // many loops use space dimensionality, it is safer to project
  // vector on the local coordinate system
  AmanziGeometry::Point face_centroid(Entity_ID f) const {
    const auto& xe = mesh_->edge_centroid(f);
    return coordsys_->Project(xe, true);
  }

  AmanziGeometry::Point face_normal(int f) const {
    const auto& tau = mesh_->edge_vector(f);
    auto normal = tau ^ coordsys_->normal_unit();
    return coordsys_->Project(normal, false);
  }

  double face_area(Entity_ID f) const {
    return mesh_->edge_length(f);
  }

  AmanziGeometry::Point edge_centroid(Entity_ID e) const {
    AmanziGeometry::Point xyz(d_);
    mesh_->node_get_coordinates(e, &xyz);
    return coordsys_->Project(xyz, true);
  }

  AmanziGeometry::Point edge_vector(Entity_ID e) const {
    AmanziGeometry::Point xyz(d_);
    return coordsys_->Project(xyz, false);
  }

  double edge_length(Entity_ID e) const { return 0.0; }

  void node_get_coordinates(Entity_ID v, AmanziGeometry::Point* xyz) const {
    AmanziGeometry::Point tmp(d_);
    mesh_->node_get_coordinates(v, &tmp);
    *xyz = coordsys_->Project(tmp, true);
  }

 private:
  Teuchos::RCP<const AmanziMesh::MeshLight> mesh_;
  std::shared_ptr<const SurfaceCoordinateSystem> coordsys_;
  int d_;
};

} // namespace WhetStone
} // namespace Amanzi

#endif
