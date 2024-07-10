#pragma once

#include "RaggedArray_DualView.hh"
#include "MeshDefs.hh"

namespace Amanzi {
namespace AmanziMesh {

struct MeshCacheData {
  // flags
  bool cell_geometry_cached = false;
  bool cell_faces_cached = false;
  bool cell_edges_cached = false;
  bool cell_nodes_cached = false;
  bool cell_coordinates_cached = false;

  bool face_geometry_cached = false;
  bool face_cells_cached = false;
  bool face_edges_cached = false;
  bool face_nodes_cached = false;
  bool face_coordinates_cached = false;

  bool edge_geometry_cached = false;
  bool edge_cells_cached = false;
  bool edge_faces_cached = false;
  bool edge_nodes_cached = false;
  bool edge_coordinates_cached = false;
  bool edge_lengths_cached = false;

  bool node_cells_cached = false;
  bool node_faces_cached = false;
  bool node_edges_cached = false;
  bool node_coordinates_cached = false;

  bool parent_entities_cached = false;
  bool cell_cellbelow_cached = false;

  bool cell_global_indices_cached = false;

  // geometry
  Point_DualView node_coordinates;
  Point_DualView cell_centroids;
  Point_DualView face_centroids;
  Point_DualView edge_centroids;
  RaggedArray_DualView<AmanziGeometry::Point> cell_coordinates;
  RaggedArray_DualView<AmanziGeometry::Point> face_coordinates;
  RaggedArray_DualView<AmanziGeometry::Point> edge_coordinates;

  Double_DualView cell_volumes;
  Double_DualView face_areas;
  Double_DualView edge_lengths;

  RaggedArray_DualView<AmanziGeometry::Point> face_normals;
  RaggedArray_DualView<Direction_type> face_normal_orientations;
  Point_DualView edge_vectors;

  // downward adjacencies
  RaggedArray_DualView<Entity_ID> cell_faces;
  RaggedArray_DualView<Direction_type> cell_face_directions;
  RaggedArray_DualView<AmanziGeometry::Point> cell_face_bisectors;

  RaggedArray_DualView<Entity_ID> cell_edges;
  RaggedArray_DualView<Entity_ID> cell_nodes;
  RaggedArray_DualView<Entity_ID> face_edges;
  RaggedArray_DualView<Direction_type> face_edge_directions;
  RaggedArray_DualView<Entity_ID> face_nodes;
  RaggedArray_DualView<Entity_ID> edge_nodes;
  Entity_ID_DualView cell_cellbelow;

  // upward adjacencies
  RaggedArray_DualView<Entity_ID> face_cells;
  RaggedArray_DualView<Entity_ID> edge_cells;
  RaggedArray_DualView<Entity_ID> edge_faces;
  RaggedArray_DualView<Entity_ID> node_cells;
  RaggedArray_DualView<Entity_ID> node_faces;
  RaggedArray_DualView<Entity_ID> node_edges;

  // parent entities
  Entity_ID_DualView parent_nodes;
  Entity_ID_DualView parent_edges;
  Entity_ID_DualView parent_faces;
  Entity_ID_DualView parent_cells;

  // Need GIDs for referencing faces
  Entity_GID_DualView cell_global_indices;

  // Need boundary entities for converting to and from standard entities.
  Entity_ID_DualView boundary_faces;
  Entity_ID_DualView boundary_nodes;

  // sizes
  Entity_ID ncells_owned, ncells_all;
  Entity_ID nfaces_owned, nfaces_all;
  Entity_ID nedges_owned, nedges_all;
  Entity_ID nnodes_owned, nnodes_all;
  Entity_ID nboundary_faces_owned, nboundary_faces_all;
  Entity_ID nboundary_nodes_owned, nboundary_nodes_all;

  int space_dim_;
  int manifold_dim_;
  bool is_ordered_;
  bool is_logical_;
  bool has_edges_;
  bool has_nodes_;
  bool has_node_faces_;
  bool is_sfm_;

  MeshCacheData()
    : node_coordinates("node_coordinates", 0),
      cell_centroids("cell_centroids", 0),
      face_centroids("face_centroids", 0),
      edge_centroids("edge_centroids", 0),
      cell_volumes("cell_volumes", 0),
      face_areas("face_areas", 0),
      edge_lengths("edge_lengths", 0),
      cell_cellbelow("cell_cellbelow", 0),
      ncells_owned(-1),
      ncells_all(-1),
      nfaces_owned(-1),
      nfaces_all(-1),
      nedges_owned(-1),
      nedges_all(-1),
      nnodes_owned(-1),
      nnodes_all(-1),
      nboundary_faces_owned(-1),
      nboundary_faces_all(-1),
      nboundary_nodes_owned(-1),
      nboundary_nodes_all(-1),
      space_dim_(-1),
      manifold_dim_(-1),
      is_ordered_(false),
      is_logical_(false),
      has_edges_(false),
      has_nodes_(true),
      has_node_faces_(true),
      is_sfm_(false)
  {}
};


} // namespace AmanziMesh
} // namespace Amanzi
