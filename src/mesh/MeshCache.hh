/*
  Copyright 2010-201x held jointly by LANL, ORNL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
           Julien Loiseau (jloiseau@lanl.gov)
*/
//! Caches mesh information for fast repeated access.

#pragma once

#include <string>

#include "AmanziComm.hh"
#include "GeometricModel.hh"
#include "MeshDefs.hh"
#include "MeshMaps_decl.hh"
#include "MeshSets.hh"
#include "MeshSets_Helpers.hh"

namespace Amanzi {
namespace AmanziMesh {

struct MeshCache {

  // ----------------------
  // Accessors and Mutators
  // ----------------------
  Comm_ptr_type get_comm() const { return comm_; }
  void set_comm(const Comm_ptr_type& comm) { comm_ = comm; }

  Teuchos::RCP<const AmanziGeometry::GeometricModel> get_geometric_model() const { return gm_; }
  void set_geometric_model(const Teuchos::RCP<const AmanziGeometry::GeometricModel>& gm) { gm_ = gm; }

  // space dimension describes the dimension of coordinates in space
  std::size_t get_space_dimension() const { return space_dim_; }
  void set_space_dimension(unsigned int dim) { space_dim_ = dim; }

  // manifold dimension describes the dimensionality of the corresponding R^n
  // manifold onto which this mesh can be projected.
  std::size_t get_manifold_dimension() const { return manifold_dim_; }
  void set_manifold_dimension(const unsigned int dim) { manifold_dim_ = dim; }

  // Some meshes are subsets of or derived from a parent mesh.
  // Usually this is null, but some meshes may provide it.
  Teuchos::RCP<const MeshCache> get_parent() const { return parent_; }
  void set_parent(const Teuchos::RCP<const MeshCache>& parent) {
    parent_ = parent; }

  // Some meshes have a corresponding mesh that is better for visualization.
  const MeshCache& get_vis_mesh() const {
    if (vis_mesh_.get()) return *vis_mesh_;
    return *this;
  }
  void set_vis_mesh(const Teuchos::RCP<const MeshCache>& vis_mesh) {
    vis_mesh_ = vis_mesh; }

  // map access
  const Entity_ID_View& get_boundary_faces() const {
    return maps_.get_boundary_faces();
  }
  const Entity_ID_View& get_boundary_nodes() const {
    return maps_.get_boundary_nodes();
  }
  const Map_type& get_map(const Entity_kind kind, bool is_ghosted) const {
    return maps_.get_map(kind, is_ghosted);
  }
  const Import_type& get_importer(const Entity_kind kind) const {
    return maps_.get_importer(kind);
  }
  const Import_type& get_boundary_face_importer() const {
    return maps_.get_boundary_face_importer();
  }
  const Import_type& get_boundary_node_importer() const {
    return maps_.get_boundary_node_importer();
  }

  // set access
  inline
  const Entity_ID_View& getSetEntities(const std::string& region_name,
          const Entity_kind kind,
          const Parallel_type ptype) const {
    auto key = std::make_tuple(region_name, kind, ptype);
    if (sets_.count(key)) {
      return sets_.at(key);
    } else {
      auto region = get_geometric_model()->FindRegion(region_name);
      sets_[key] = resolveMeshSet(*region, kind, ptype, *this);
      return sets_.at(key);
    }
  }

  // ----------------
  // Entity meta-data
  // ----------------
  std::size_t getNumEntities(const Entity_kind kind, const Parallel_type ptype) const {
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

  // corresponding entity in the parent mesh
  //
  // Note the kind refers to the kind in _this_ mesh -- for some lifted meshes,
  // this may not be the same as the entity kind in the parent mesh.  That
  // logic is left to the user of this class -- we simply store the IDs.
  Entity_ID getEntityParent(const Entity_kind kind, const Entity_ID entid) const {
    AMANZI_ASSERT(parents_initialized);
    switch(kind) {
      case Entity_kind::CELL:
        return parent_cells_[entid];
        break;
      case Entity_kind::FACE:
        return parent_faces_[entid];
        break;
      case Entity_kind::EDGE:
        return parent_edges_[entid];
        break;
      case Entity_kind::NODE:
        return parent_nodes_[entid];
      default: {}
    }
    return -1;
  }


  //---------------------
  // Geometry
  //---------------------
  // locations
  AmanziGeometry::Point getNodeCoordinate(const Entity_ID n) const {
    return node_coordinates[n];
  }
  AmanziGeometry::Point getCellCentroid(const Entity_ID c) const {
    return cell_centroids[c];
  }
  AmanziGeometry::Point getFaceCentroid(const Entity_ID f) const {
    return face_centroids[f];
  }
  AmanziGeometry::Point getEdgeCentroid(const Entity_ID e) const {
    return edge_centroids[e];
  }
  AmanziGeometry::Point getCentroid(const Entity_kind kind, const Entity_ID ent) {
    switch(kind) {
      case (Entity_kind::CELL) :
        return getCellCentroid(ent);
        break;
      case (Entity_kind::FACE) :
        return getFaceCentroid(ent);
        break;
      case (Entity_kind::EDGE) :
        return getEdgeCentroid(ent);
        break;
      case (Entity_kind::NODE) :
        return getNodeCoordinate(ent);
        break;
      default :
        Errors::Message msg("Invalid argument kind to getCentroid");
        Exceptions::amanzi_throw(msg);
    }
    return AmanziGeometry::Point();
  }

  // extent
  double getCellVolume(const Entity_ID c) const {
    return cell_volumes[c];
  }
  double getFaceArea(const Entity_ID f) const {
    return face_areas[f];
  }
  double getEdgeLength(const Entity_ID e) const {
    return edge_lengths[e];
  }

  // Normal vector of a face
  inline
  AmanziGeometry::Point getFaceNormal(const Entity_ID f) const {
    return face_normals[f][0];
  }

  // Normal vector and natural direction of a face, outward with respect to a
  // cell.
  //
  // The vector is normalized and then weighted by the area of the face.
  //
  // The orientation is 1 if the outward normal is the same direction as the
  // natural normal, -1 if in opposite directions, and 0 if there is no natural
  // normal.
  std::pair<AmanziGeometry::Point, int> getFaceNormal(const Entity_ID f,
          const Entity_ID c) const;

  // Vector describing the edge, where the length is the edge length.
  //
  // Orientation is the natural orientation, e.g. that it points from node 0 to
  // node 1 with respect to edge_node adjacency information.
  inline
  AmanziGeometry::Point getEdgeVector(const Entity_ID e) const {
    return edge_vectors[e];
  }

  //---------------------
  // Downward adjacencies
  //---------------------
  // Get faces of a cell
  //
  // On a distributed mesh, this will return all the faces of the
  // cell, OWNED or GHOST. If the framework supports it, the faces will be
  // returned in a standard order according to Exodus II convention
  // for standard cells; in all other situations (not supported or
  // non-standard cells), the list of faces will be in arbitrary order
  //
  // In 3D manifolds and 2D manifolds of 2d space, the natural direction of a
  // face is defined by the right-hand-rule of the node ordering.  In 2D
  // manifolds of 3D space, there is no natural direction.
  //
  // Use of bisectors instead of coordinate geometry enables use of logical
  // meshes which may not have geometric coordinates.
  //
  // New interfaces that return by const reference should be preferred rather
  // than those that return void and expect the return value as an argument --
  // this new-style interface works better with Kokkos and should be more
  // efficient in all cases.
  std::size_t getCellNumFaces(const Entity_ID c) const {
    return cell_faces[c].size();
  }
  const Entity_ID_View& getCellFaces(const Entity_ID c) const {
    return cell_faces[c];
  }
  const Entity_Direction_View& getCellFaceDirections(const Entity_ID c) const {
    return cell_face_directions[c];
  }
  const Point_View& getCellFaceBisectors(const Entity_ID c) const {
    return cell_face_bisectors[c];
  }

  [[deprecated("Prefer to use non-void variant that returns faces directly")]]
  void getCellFaces(const Entity_ID c,
                    Entity_ID_View& faces) const {
    getCellFacesAndDirs(c, faces, nullptr);
  }

  [[deprecated("Prefer to use non-void variant that returns faces directly")]]
  void getCellFacesAndDirs(const Entity_ID c,
                           Entity_ID_View& faces,
                           Entity_Direction_View * const dirs) const {
    if (dirs) *dirs = getCellFaceDirections(c);
    faces = getCellFaces(c);
  }

  [[deprecated("Prefer to use non-void variant that returns faces directly")]]
  void getCellFacesAndBisectors(
          const Entity_ID cellid,
          Entity_ID_View& faceids,
          Point_List * const bisectors) const;

  // Get edges of a cell.
  std::size_t getCellNumEdges(const Entity_ID c) const {
    return cell_edges[c].size();
  }
  const Entity_ID_View& getCellEdges(const Entity_ID c) const {
    return cell_edges[c];
  }
  [[deprecated("Prefer to use non-void variant that returns edges directly")]]
  void getCellEdges(const Entity_ID c, Entity_ID_View& edges) const {
    edges = getCellEdges(c);
  }

  // Get nodes of a cell.
  std::size_t getCellNumNodes(const Entity_ID c) const {
    return cell_nodes[c].size();
  }
  const Entity_ID_View& getCellNodes(const Entity_ID c) const {
    return cell_nodes[c];
  }
  [[deprecated("Prefer to use non-void variant that returns nodes directly")]]
  void getCellNodes(const Entity_ID c, Entity_ID_View& nodes) const {
    nodes = getCellNodes(c);
  }

  // Get edges of a face and directions in which the face uses the edges.
  //
  // In 3D, edge direction is 1 when it is oriented counter clockwise
  // with respect to the face natural normal.
  //
  // On a distributed mesh, this will return all the edges of the
  // face, OWNED or GHOST. If the framework supports it, the edges will be
  // returned in a ccw order around the face as it is naturally defined.
  //
  // IMPORTANT NOTE IN 2D: In meshes where the cells are two
  // dimensional, faces and edges are identical. For such cells, this
  // operator will return a single edge and a direction of 1. However,
  // this direction cannot be relied upon to compute, say, a contour
  // integral around the 2D cell.
  std::size_t getFaceNumEdges(const Entity_ID f) const {
    return face_edges[f].size();
  }
  const Entity_ID_View& getFaceEdges(const Entity_ID f) const {
    return face_edges[f];
  }
  const Entity_Direction_View& getFaceEdgeDirections(const Entity_ID f) const {
    return face_edge_directions[f];
  }
  [[deprecated("Prefer to use non-void variant that returns edges directly")]]
  void getFaceEdgesAndDirs(const Entity_ID f,
                           Entity_ID_View& edges,
                           Entity_Direction_View * const dirs=nullptr) const {
    edges = getFaceEdges(f);
    if (dirs) *dirs = getFaceEdgeDirections(f);
  }

  // Get nodes of face
  //
  // In 3D, the nodes of the face are returned in ccw order consistent
  // with the face normal.
  std::size_t getFaceNumNodes(const Entity_ID f) const {
    return face_nodes[f].size();
  }
  const Entity_ID_View& getFaceNodes(const Entity_ID f) const {
    return face_nodes[f];
  }
  [[deprecated("Prefer to use non-void variant that returns nodes directly")]]
  void getFaceNodes(const Entity_ID f, Entity_ID_View& nodes) const {
    nodes = getFaceNodes(f);
  }

  // Get nodes of edge
  std::size_t getEdgeNumNodes(const Entity_ID e) const {
    return edge_nodes[e].size();
  }
  const Entity_ID_View& getEdgeNodes(const Entity_ID e) const {
    return edge_nodes[e];
  }
  [[deprecated("Prefer to use non-void variant that returns nodes directly")]]
  void getEdgeNodes(const Entity_ID e, Entity_ID_View& nodes) const;

  //-------------------
  // Upward adjacencies
  //-------------------
  // The cells are returned in no particular order. Also, the order of cells
  // is not guaranteed to be the same for corresponding faces on different
  // processors
  std::size_t getFaceNumCells(const Entity_ID f) const {
    return face_cells[f].size();
  }
  const Entity_ID_View& getFaceCells(const Entity_ID f) const {
    return face_cells[f];
  }
  void getFaceCells(const Entity_ID f,
                    const Parallel_type ptype,
                    Entity_ID_View& cells) const {
    if (ptype == Parallel_type::ALL) {
      cells = getFaceCells(f);
    } else if (ptype == Parallel_type::OWNED) {
      auto all_cells = getFaceCells(f);
      cells.resize(0);
      cells.reserve(all_cells.size());
      for (const auto& c : all_cells) {
        if (c < ncells_owned) cells.push_back(c);
      }
    } else if (ptype == Parallel_type::GHOST) {
      auto all_cells = getFaceCells(f);
      cells.resize(0);
      cells.reserve(all_cells.size());
      for (const auto& c : all_cells) {
        if (c >= ncells_owned) cells.push_back(c);
      }
    } else {
      cells.resize(0);
    }
  }

  // Cells of a given Parallel_type connected to an edge
  //
  // The order of cells is not guaranteed to be the same for corresponding
  // edges on different processors
  void getEdgeCells(const Entity_ID edgeid,
                    const Parallel_type ptype,
                    Entity_ID_View& cellids) const {
    Errors::Message msg("MeshCache::getEdgeCells not implemented");
    Exceptions::amanzi_throw(msg);
  }

  // Faces of type 'ptype' connected to an edge
  // NOTE: The order of faces is not guaranteed to be the same for
  // corresponding edges on different processors
  void getEdgeFaces(const Entity_ID edgeid,
                            const Parallel_type ptype,
                    Entity_ID_View& faceids) const {
    Errors::Message msg("MeshCache::getEdgeFaces not implemented");
    Exceptions::amanzi_throw(msg);
  }

  // Cells of type 'ptype' connected to a node
  // NOTE: The order of cells is not guaranteed to be the same for
  // corresponding nodes on different processors
  void getNodeCells(const Entity_ID nodeid,
                    const Parallel_type ptype,
                    Entity_ID_View& cellids) const {
    Errors::Message msg("MeshCache::getNodeCells not implemented");
    Exceptions::amanzi_throw(msg);
  }

  // Faces of type parallel 'ptype' connected to a node
  // NOTE: The order of faces is not guarnateed to be the same for
  // corresponding nodes on different processors
  void getNodeFaces(const Entity_ID nodeid,
                    const Parallel_type ptype,
                    Entity_ID_View& faceids) const {
    Errors::Message msg("MeshCache::getNodeFaces not implemented");
    Exceptions::amanzi_throw(msg);
  }

  // Edges of type 'ptype' connected to a node
  //
  // The order of edges is not guaranteed to be the same for corresponding
  // node on different processors
  void getNodeEdges(const Entity_ID nodeid,
                    const Parallel_type ptype,
                    Entity_ID_View& edgeids) const {
    Errors::Message msg("MeshCache::getNodeEdges not implemented");
    Exceptions::amanzi_throw(msg);
  }


  // all of this is public to allow direct, fast access to data.
  // sizes
  Entity_ID ncells_owned, ncells_all;
  Entity_ID nfaces_owned, nfaces_all;
  Entity_ID nedges_owned, nedges_all;
  Entity_ID nnodes_owned, nnodes_all;
  Entity_ID nboundary_faces_owned, nboundary_faces_all;
  Entity_ID nboundary_nodes_owned, nboundary_nodes_all;

  // geometry
  Point_View node_coordinates;
  Point_View cell_centroids;
  Point_View face_centroids;
  Point_View edge_centroids;

  Double_View cell_volumes;
  Double_View face_areas;
  Double_View edge_lengths;

  CSR<AmanziGeometry::Point> face_normals;
  Point_View edge_vectors;

  // downward adjacencies
  CSR<Entity_ID> cell_faces;
  CSR<int> cell_face_directions;
  CSR<AmanziGeometry::Point> cell_face_bisectors;

  CSR<Entity_ID> cell_edges;
  CSR<Entity_ID> cell_nodes;
  CSR<Entity_ID> face_edges;
  CSR<int> face_edge_directions;
  CSR<Entity_ID> face_nodes;
  CSR<Entity_ID> edge_nodes;

  // upward adjacencies
  CSR<Entity_ID> face_cells;
  CSR<Entity_ID> edge_cells;
  CSR<Entity_ID> edge_faces;
  CSR<Entity_ID> node_cells;
  CSR<Entity_ID> node_faces;
  CSR<Entity_ID> node_edges;

  // parent entities
  Entity_ID_View parent_nodes_;
  Entity_ID_View parent_edges_;
  Entity_ID_View parent_faces_;
  Entity_ID_View parent_cells_;

  // initialized flags
  bool nodes_initialized;
  bool edges_initialized;
  bool faces_initialized;
  bool cells_initialized;
  bool parents_initialized;

  // standard things
  Comm_ptr_type comm_;
  Teuchos::RCP<const AmanziGeometry::GeometricModel> gm_;
  int space_dim_;
  int manifold_dim_;

  // related meshes
  Teuchos::RCP<const MeshCache> parent_;
  Teuchos::RCP<const MeshCache> vis_mesh_;

  // helper classes
  MeshMaps maps_;
  mutable MeshSets sets_;

};




} // namespace AmanziMesh
} // namespace Amanzi
