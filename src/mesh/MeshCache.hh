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

/*

Design discussion:
------------------

Goals:

1. Allow O(1) access to all adjacencies, coordinates, and geometric quantities
   that is as fast as array access, if required.
2. Allow the underlying MeshFramework to be deleted.
3. Avoid caching everything.
4. Allow things not cached to still be calculated, even if the framework has
   been deleted, from existing cached data.

Effectively this set of requirements means there are two or three ways of
getting each piece of information:

1. the MeshCache's cache (e.g. "fast" access)
2. directly from the Framework (e.g. "framework" access)
3. recomputing from some mix of "cache" and "framework" info (e.g. "slow" access)

Note that the third is not always possible -- for instance, if cell-face
adjacencies are not cached and the framework is deleted, there is no way to
recompute these adjacencies.  In this case, the third access method will not
exist.

"Fast" access will have a few interfaces to allow a method-based access (the
standard mesh interface), but MeshCache will also expose all underlying data as
public, allowing direct access.  This "breaking encapsulation" is particularly
crucial for novel architectures, where we might want to perform a kernel launch
over a view of the data.

In general, all things are returned by const reference.  This can be broken in
a few cases, as a const container does not imply that the contained things are
const.  This should not be a common problem though.  When this class migrates
to Kokkos::Views, these will return by value instead of by const reference,
matching the Kokkos View semantics.


Therefore, a standard user pattern might look like PDE_OperatorFV, which needs
face-cell (bidirectional) adjacencies, face normals, cell volumes, and
cell-face bisectors, but doesn't really need coordinates, edges, etc.  This
would result in only using "fast" access patterns.

.. code::

   auto& m = *S->GetMesh(...);
   m.cacheCellGeometry();
   m.cacheCellFaces()
   m.cacheFaceCells()
   m.cacheFaceGeometry();
   m.DestroyFramework();  // destroys Framework after needs are cached

   ...

   // ONE WAY: get the view
   for (Entity_ID f=0; f!=m.nfaces_owned; ++f) {
     const Entity_ID_View& fcells = m.getFaceCells(f);
     for (int i=0; i!=fcells.size(); ++i) {
       ... with fcells[i] ...
     }
   }

   // OR ANOTHER: direct access
   for (Entity_ID f=0; f!=m.nfaces_owned; ++f) {
     for (int i=0; i!=m.getFaceNumCells(f); ++i) {
       ... with m.face_cells[f][i] ..
     }
   }


Alternatively, for instance, PDE_OperatorMFD may need cell coordinates to
calculate mass matrices at the start, but, as long as the mesh isn't deforming,
may not need to do this often.  Therefore, the user might NOT want to cache
cell coordinates, because this uses a lot of memory.  Then the user could do
one of two things -- either don't destroy the Framework mesh, or (probably
preferably) re-construct the cell coordinate lists each time they are accessed
from the cell-node adjacencies and the node coordinates.  This would be the
"slow" access pattern, which recomputes as needed.

.. code::

   auto& m = *S->GetMesh(...);
   m.cacheCellGeometry();
   // m.cacheCellNodes() // Don't cache the cell-to-node adjacencies.
   m.cacheFaceNodes() // Cell-to-node adjacencies will be reconstructed from
   m.cacheCellFaces() // cell-to-face and face-to-node caches.
   m.cacheFaceCells()
   m.cacheFaceGeometry();
   m.cacheNodeCoordinates();
   m.DestroyFramework();  // destroys Framework after needs are cached

   // This would be valid if we HAD called m.cacheCellCoordinates()
   for (Entity_ID c=0; c!=m.ncells_owned; ++c) {
     auto ccoords = m.getCellCoordinates(c); // SEG FAULT, not cached!
     ...
   }

   // This would be valid if we had NOT called m.DestroyFramework()
   for (Entity_ID c=0; c!=m.ncells_owned; ++c) {
     auto ccoords = m.getCellCoordinates<MeshCache::AccessPattern::FRAMEWORK>(c);
     // SEG FAULT, framework mesh destroyed
     ...
   }

   // This will do the slow thing of iterating over all faces in cells, nodes
   // in faces, creating a set of nodes, and creating the Point_View of those
   // node coordinates.
   for (Entity_ID c=0; c!=m.ncells_owned; ++c) {
     auto ccords = m.getCellCoordinates<MeshCache::AccessPattern::RECOMPUTE>(c);
     ...
   }


These snippets hint at the underlying mechanism -- most methods accept a
template parameter for the AccessPattern, which is an enum which may be:

1. AccessPattern::CACHE ("fast" access, no error checking, this is the default!)
2. AccessPattern::RECOMPUTE ("slow" access but will work if at all possible)
3. AccessPattern::FRAMEWORK (go back to the framework mesh to get this quantity)

Note also that users may query whether something was cached or not via, e.g.,
cell_coordinates_cached boolean.

*/


#pragma once

#include <string>

#include "AmanziComm.hh"
#include "GeometricModel.hh"
#include "MeshDefs.hh"
#include "MeshSets.hh"
#include "MeshMaps_decl.hh"

namespace Amanzi {
namespace AmanziMesh {

class MeshFramework;

struct MeshCache {

  MeshCache();

  //
  // The current assumption is that all MeshCaches are made from a framework
  // mesh.
  //
  // This may change.
  //
  MeshCache(const Teuchos::RCP<MeshFramework>& framework_mesh, bool natural=false)
    :  MeshCache() {
    set_framework_mesh(framework_mesh);
    InitializeFramework(natural);
  }
  void DestroyFramework() {
    framework_mesh_ = Teuchos::null;
  }
  void InitializeFramework(bool natural=false);
  const MeshFramework& get_framework_mesh() const { return *framework_mesh_; }
  MeshFramework& get_framework_mesh() { return *framework_mesh_; }
  void set_framework_mesh(const Teuchos::RCP<MeshFramework>& framework_mesh) {
    framework_mesh_ = framework_mesh; }

  //
  // Build the cache, fine grained control
  // =============================================

  // cell centroid, volume
  void cacheCellGeometry();
  bool cell_geometry_cached;
  // cell-face adjacencies
  void cacheCellFaces();
  bool cell_faces_cached;
  // cell-edge adjacencies
  void cacheCellEdges();
  bool cell_edges_cached;
  // cell-node adjacencies
  void cacheCellNodes();
  bool cell_nodes_cached;
  // cell coordinates
  void cacheCellCoordinates();
  bool cell_coordinates_cached;

  // face centroid, area, normals
  void cacheFaceGeometry();
  bool face_geometry_cached;
  // face-cell adjacencies
  void cacheFaceCells();
  bool face_cells_cached;
  // face-edge adjacencies
  void cacheFaceEdges();
  bool face_edges_cached;
  // face-node adjacencies
  void cacheFaceNodes();
  bool face_nodes_cached;
  // face coordinates
  void cacheFaceCoordinates();
  bool face_coordinates_cached;

  // edge centroid, length, vector
  void cacheEdgeGeometry();
  bool edge_geometry_cached;
  // edge-cell adjacencies
  void cacheEdgeCells();
  bool edge_cells_cached;
  // edge-face adjacencies
  void cacheEdgeFaces();
  bool edge_faces_cached;
  // edge-node adjacencies
  void cacheEdgeNodes();
  bool edge_nodes_cached;
  // edge coordinates
  void cacheEdgeCoordinates();
  bool edge_coordinates_cached;

  // node-cell adjacencies
  void cacheNodeCells();
  bool node_cells_cached;
  // node-face adjacencies
  void cacheNodeFaces();
  bool node_faces_cached;
  // node-edge adjacencies
  void cacheNodeEdges();
  bool node_edges_cached;
  // node coordinates
  void cacheNodeCoordinates();
  bool node_coordinates_cached;

  // Note that regions are cached on demand the first time they are requested,
  // but labeled sets must be pre-cached if the framework mesh is to be
  // destroyed.
  void precacheLabeledSets();

  // build the MeshMaps object
  //
  // this is called in InitializeFramework and may not need to be public?
  void buildMaps(bool natural=false);

  //
  // Baseline mesh functionality
  // =============================================

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

  // mesh properties
  bool is_ordered() const { return is_ordered_; }
  bool has_edges() const { return has_edges_; }

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
  const Entity_ID_View& getSetEntities(const std::string& region_name,
          const Entity_kind kind,
          const Parallel_type ptype) const;

  // ----------------
  // Entity meta-data
  // ----------------
  std::size_t getNumEntities(const Entity_kind kind, const Parallel_type ptype) const;

  // corresponding entity in the parent mesh
  //
  // Note the kind refers to the kind in _this_ mesh -- for some lifted meshes,
  // this may not be the same as the entity kind in the parent mesh.  That
  // logic is left to the user of this class -- we simply store the IDs.
  inline Entity_ID getEntityParent(const Entity_kind kind, const Entity_ID entid) const {
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
    return -1;
  }

  Cell_type getCellType(const Entity_ID c) const;

  //---------------------
  // Geometry
  //---------------------
  // locations
  const AmanziGeometry::Point& getNodeCoordinate(const Entity_ID n) const {
    return node_coordinates[n];
  }
  void setNodeCoordinate(const Entity_ID n, const AmanziGeometry::Point& coord);

  const Point_View& getEdgeCoordinates(const Entity_ID n) const {
    return edge_coordinates[n];
  }
  const Point_View& getFaceCoordinates(const Entity_ID n) const {
    return face_coordinates[n];
  }
  const Point_View& getCellCoordinates(const Entity_ID n) const {
    return cell_coordinates[n];
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
    return AmanziGeometry::norm(face_normals[f][0]);
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
  inline
  AmanziGeometry::Point getFaceNormal(const Entity_ID f,
          const Entity_ID c, int* orientation=nullptr) const {
    int i;
    for (i=0; i!=face_normals[f].size(); ++i)
      if (face_cells[f][i] == c) break;

    if (orientation) *orientation = face_normal_directions[f][i];
    return face_normals[f][i];
  }

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
  inline
  std::size_t getCellNumFaces(const Entity_ID c) const {
    return cell_faces[c].size();
  }
  inline
  const Entity_ID_View& getCellFaces(const Entity_ID c) const {
    return cell_faces[c];
  }
  inline
  const Entity_Direction_View& getCellFaceDirections(const Entity_ID c) const {
    return cell_face_directions[c];
  }
  inline
  const Point_View& getCellFaceBisectors(const Entity_ID c) const {
    return cell_face_bisectors[c];
  }

  // NOTE: all deprecated pragmas should go back in after finished refactoring
  // [[deprecated("Prefer to use non-void variant that returns faces directly")]]
  inline
  void getCellFaces(const Entity_ID c,
                    Entity_ID_View& faces) const {
    getCellFacesAndDirs(c, faces, nullptr);
  }

  //[[deprecated("Prefer to use non-void variant that returns faces directly")]]
  inline
  void getCellFacesAndDirs(const Entity_ID c,
                           Entity_ID_View& faces,
                           Entity_Direction_View * const dirs) const {
    if (dirs) *dirs = getCellFaceDirections(c);
    faces = getCellFaces(c);
  }

  //[[deprecated("Prefer to use non-void variant that returns faces directly")]]
  void getCellFacesAndBisectors(
          const Entity_ID cellid,
          Entity_ID_View& faceids,
          Point_List * const bisectors) const;

  // Get edges of a cell.
  inline
  std::size_t getCellNumEdges(const Entity_ID c) const {
    return cell_edges[c].size();
  }
  inline
  const Entity_ID_View& getCellEdges(const Entity_ID c) const {
    return cell_edges[c];
  }
  //[[deprecated("Prefer to use non-void variant that returns edges directly")]]
  inline
  void getCellEdges(const Entity_ID c, Entity_ID_View& edges) const {
    edges = getCellEdges(c);
  }

  // Get nodes of a cell.
  inline
  std::size_t getCellNumNodes(const Entity_ID c) const {
    return cell_nodes[c].size();
  }
  inline
  const Entity_ID_View& getCellNodes(const Entity_ID c) const {
    return cell_nodes[c];
  }
  //[[deprecated("Prefer to use non-void variant that returns nodes directly")]]
  inline
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
  inline
  std::size_t getFaceNumEdges(const Entity_ID f) const {
    return face_edges[f].size();
  }
  inline
  const Entity_ID_View& getFaceEdges(const Entity_ID f) const {
    return face_edges[f];
  }
  inline
  const Entity_Direction_View& getFaceEdgeDirections(const Entity_ID f) const {
    return face_edge_directions[f];
  }
  //[[deprecated("Prefer to use non-void variant that returns edges directly")]]
  inline
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
  inline
  std::size_t getFaceNumNodes(const Entity_ID f) const {
    return face_nodes[f].size();
  }
  inline
  const Entity_ID_View& getFaceNodes(const Entity_ID f) const {
    return face_nodes[f];
  }
  //[[deprecated("Prefer to use non-void variant that returns nodes directly")]]
  inline
  void getFaceNodes(const Entity_ID f, Entity_ID_View& nodes) const {
    nodes = getFaceNodes(f);
  }

  // Get nodes of edge
  inline
  std::size_t getEdgeNumNodes(const Entity_ID e) const {
    return edge_nodes[e].size();
  }
  inline
  const Entity_ID_View& getEdgeNodes(const Entity_ID e) const {
    return edge_nodes[e];
  }
  //[[deprecated("Prefer to use non-void variant that returns nodes directly")]]
  inline
  void getEdgeNodes(const Entity_ID e, Entity_ID_View& nodes) const {
    nodes = getEdgeNodes(e);
  }

  //-------------------
  // Upward adjacencies
  //-------------------
  // The cells are returned in no particular order. Also, the order of cells
  // is not guaranteed to be the same for corresponding faces on different
  // processors
  inline
  std::size_t getFaceNumCells(const Entity_ID f) const {
    return face_cells[f].size();
  }
  inline
  const Entity_ID_View& getFaceCells(const Entity_ID f) const {
    return face_cells[f];
  }
  inline
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
  inline
  std::size_t getEdgeNumCells(const Entity_ID e) const {
    return edge_cells[e].size();
  }
  inline
  const Entity_ID_View& getEdgeCells(const Entity_ID e) const {
    return edge_cells[e];
  }
  inline
  void getEdgeCells(const Entity_ID e,
                    const Parallel_type ptype,
                    Entity_ID_View& cells) const {
    if (ptype == Parallel_type::ALL) {
      cells = getEdgeCells(e);
    } else if (ptype == Parallel_type::OWNED) {
      auto all_cells = getEdgeCells(e);
      cells.resize(0);
      cells.reserve(all_cells.size());
      for (const auto& c : all_cells) {
        if (c < ncells_owned) cells.push_back(c);
      }
    } else if (ptype == Parallel_type::GHOST) {
      auto all_cells = getEdgeCells(e);
      cells.resize(0);
      cells.reserve(all_cells.size());
      for (const auto& c : all_cells) {
        if (c >= ncells_owned) cells.push_back(c);
      }
    } else {
      cells.resize(0);
    }
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
  CSR<AmanziGeometry::Point> cell_coordinates;
  CSR<AmanziGeometry::Point> face_coordinates;
  CSR<AmanziGeometry::Point> edge_coordinates;

  Double_View cell_volumes;
  Double_View edge_lengths;

  CSR<AmanziGeometry::Point> face_normals;
  CSR<int> face_normal_directions;
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
  Entity_ID_View parent_nodes;
  Entity_ID_View parent_edges;
  Entity_ID_View parent_faces;
  Entity_ID_View parent_cells;

  // initialized flags
  bool nodes_initialized;
  bool edges_initialized;
  bool faces_initialized;
  bool cells_initialized;
  bool parents_initialized;

 private:

  // standard things
  Comm_ptr_type comm_;
  Teuchos::RCP<const AmanziGeometry::GeometricModel> gm_;
  int space_dim_;
  int manifold_dim_;
  bool is_ordered_;
  bool has_edges_;

  // related meshes
  Teuchos::RCP<MeshFramework> framework_mesh_;
  Teuchos::RCP<const MeshCache> parent_;
  Teuchos::RCP<const MeshCache> vis_mesh_;

  // helper classes
  MeshMaps maps_;
  mutable MeshSets sets_;

};


namespace Utils {

void cacheAll(MeshCache& mesh);
void recacheGeometry(MeshCache& mesh);

} // namespace Testing
} // namespace AmanziMesh
} // namespace Amanzi