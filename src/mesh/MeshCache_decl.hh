/*
  Copyright 2010-202x held jointly by participating institutions.
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
2. recomputing using a mesh algorithm and the public interface of MeshCache.
3. directly from the MeshFramework (e.g. "framework" access), which may itself
   call a mesh algorithm and its own interface.

Note that the second is not always possible -- for instance, node coordinates
cannot be recomputed; nor can basic adjacency information, as there is no
requirement here that a minimal representation of the framework has been cached
before it is destroyed.  But many of these _can_ be recomputed; for instance,
one can determine cell coordinates from node coordinates and cell-node
adjacency, which itself can be computed from cell-face and face-node adjacency
information.  Therefore, it is likely faster to recompute with MeshCache, using
other cached values, than to ask the framework, which likely will use the same
algorithm but using uncached values.

So the default API typically tries to use a cached value first (if it is
available), then falls back on either the framework (for fundamental things
that probably shouldn't or can't be computed) or the algorithm (for everything
else).

The framework is split into two classes: MeshFramework, which supplies all
topological information, nodal coordinates, and a fully featured API, and
MeshAlgorithms, which is a struct of algorithmic helper functions.
This split is done so that one can destroy the MeshFramework but still use the
original algorithms.

This class is also designed to be used on either a MemSpace_kind::HOST or
MemSpace_kind::DEVICE.  To do this, the cache itself is all dual views, and
these are stored in a struct that can be shared across all MeshCaches of the
same mesh.  Then, the DEVICE template parameter is used to get the "right" view
on any given device.  For now, all will be HOST.


As an example, a standard user pattern might look like PDE_OperatorFV, which
needs face-cell (bidirectional) adjacencies, face normals, cell volumes, and
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
     auto fcells = m.getFaceCells(f);
     for (int i=0; i!=fcells.size(); ++i) {
       ... with fcells[i] ...
     }
   }

   // OR ANOTHER: direct access
   for (Entity_ID f=0; f!=m.nfaces_owned; ++f) {
     auto ncells = m.getFaceNumCells<CACHED>(f)
     for (size_type i=0; i!=ncells; ++i) {
       auto&& c = m.getFaceCell(f, i);
     }
   }

Note the use of `auto&& c = ...` which allows the user to write code that
doesn't care whether the returned thing is returned by value or reference.
This is more important for Point objects, for instance.


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

   // The default will figure out a way to make it work
   //
   // This will do the slow thing of iterating over all faces in cells, nodes
   // in faces, creating a set of nodes, and creating the Point_View of those
   // node coordinates.
   for (Entity_ID c=0; c!=m.ncells_owned; ++c) {
     auto ccoords = m.getCellCoordinates(c);
     ...
   }

   // This would be an out of bounds error, because we did not cache
   for (Entity_ID c=0; c!=m.ncells_owned; ++c) {
     auto ccoords = m.getCellCoordinates<CACHED>(c);
     // out of bounds error
     ...
   }

Note that the MeshCache stores pointers to many other classes that are useful
for layering in the various capabilities that make this a fully-featured mesh.

1. MeshFramework -- the virtual "framework mesh" which may be slow, but
   provides all topological and geometric functionality.

2. MeshAlgorithms -- this is a virtual struct containing a few virtual
   methods for implementing algorithms on either this or the MeshFramework
   object.  This little class is crucial to ensure that, whether you use the
   algorithm on the framework or on the cache, you do the same algorithm.  It
   also allows alternative frameworks to overwrite these algorithms by writing
   derived classes with mesh-specific algorithms.  This would have been
   implemented as virtual methods in the MeshFramework class, but we want to be
   able to delete the framework mesh, and will still need these algorithms
   after deleting that class if not everything is cached.

3. MeshSets -- this is a dictionary used to store cached lists of entity IDs
   that exist in a given region on a given mesh.  Helpers and this using
   declaration are defined in MeshSets.hh

4. MeshMaps -- a helper class that implements providing maps for use in data
   structures.

5. MeshColumns -- a helper class that implements columnar lists of cells/faces
   to define subgrid, 1D processes in the vertical dimension.  This class is
   only created if buildColumns() is called, and has various preconditions on
   this mesh.

*/

#pragma once

#include <string>

#include "AmanziComm.hh"
#include "AmanziMap.hh"
#include "GeometricModel.hh"

#include "MeshDefs.hh"
#include "MeshSets.hh"
#include "MeshColumns.hh"

namespace Amanzi {
namespace AmanziMesh {

class MeshFramework;
struct MeshAlgorithms;

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
};


struct MeshCacheBase {
 protected:
  MeshCacheBase()
    : ncells_owned(-1),
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


  MeshCacheBase(const Teuchos::RCP<Teuchos::ParameterList>& plist) : MeshCacheBase()
  {
    if (plist == Teuchos::null) {
      plist_ = Teuchos::rcp(new Teuchos::ParameterList());
    } else {
      plist_ = plist;
    }
  }

  MeshCacheBase(const Teuchos::RCP<MeshFramework>& framework_mesh,
                const Teuchos::RCP<MeshAlgorithms>& framework_algorithms,
                const Teuchos::RCP<Teuchos::ParameterList>& plist)
    : MeshCacheBase(plist)
  {
    if (framework_algorithms == Teuchos::null) {
      Errors::Message msg("MeshAlgorithms is Teuchos::null.");
      Exceptions::amanzi_throw(msg);
    }
    algorithms_ = framework_algorithms;
  }

  MeshCacheBase(const MeshCacheBase& other) = default;

 public:
  // sizes
  Entity_ID ncells_owned, ncells_all;
  Entity_ID nfaces_owned, nfaces_all;
  Entity_ID nedges_owned, nedges_all;
  Entity_ID nnodes_owned, nnodes_all;
  Entity_ID nboundary_faces_owned, nboundary_faces_all;
  Entity_ID nboundary_nodes_owned, nboundary_nodes_all;

  //
  // Baseline mesh functionality
  // =============================================
  // mesh properties
  bool isOrdered() const { return is_ordered_; }
  bool isLogical() const { return is_logical_; }
  bool isSFM() const { return is_sfm_; } // single face mesh -- special case

  bool hasNodes() const { return has_nodes_; }
  bool hasEdges() const { return has_edges_; }
  bool hasNodeFaces() const { return has_node_faces_; }

  void hasEdgesOrThrow() const
  {
    if (!hasEdges()) {
      Errors::Message msg("MeshFramework does not include edges.");
      Exceptions::amanzi_throw(msg);
    }
  }

  // ----------------------
  // Accessors and Mutators
  // ----------------------
  Teuchos::RCP<const MeshFramework> getMeshFramework() const { return framework_mesh_; }
  Teuchos::RCP<MeshFramework> getMeshFramework() { return framework_mesh_; }
  void destroyFramework() { framework_mesh_ = Teuchos::null; }
  Teuchos::RCP<Teuchos::ParameterList> getParameterList() const { return plist_; }

  // algorithms class, a partner to the framework, that defines how to compute things
  Teuchos::RCP<const MeshAlgorithms> getAlgorithms() const { return algorithms_; }

  Comm_ptr_type getComm() const { return comm_; }
  void setComm(const Comm_ptr_type& comm) { comm_ = comm; }

  Teuchos::RCP<const AmanziGeometry::GeometricModel> getGeometricModel() const { return gm_; }
  void setGeometricModel(const Teuchos::RCP<const AmanziGeometry::GeometricModel>& gm) { gm_ = gm; }

  // space dimension describes the dimension of coordinates in space
  std::size_t getSpaceDimension() const { return space_dim_; }
  void setSpaceDimension(unsigned int dim) { space_dim_ = dim; }

  // manifold dimension describes the dimensionality of the corresponding R^n
  // manifold onto which this mesh can be projected.
  std::size_t getManifoldDimension() const { return manifold_dim_; }
  void setManifoldDimension(const unsigned int dim) { manifold_dim_ = dim; }

 protected:
  // standard things
  Comm_ptr_type comm_;
  Teuchos::RCP<Teuchos::ParameterList> plist_;
  Teuchos::RCP<const AmanziGeometry::GeometricModel> gm_;
  int space_dim_;
  int manifold_dim_;
  bool is_ordered_;
  bool is_logical_;
  bool has_edges_, has_nodes_;
  bool has_node_faces_;
  bool is_sfm_;

  // related meshes
  Teuchos::RCP<MeshFramework> framework_mesh_;
  Teuchos::RCP<const MeshAlgorithms> algorithms_;

  // helper classes
  MeshCacheData data_;
  MeshMaps maps_;
  mutable MeshSets sets_;
  mutable MeshSetVolumeFractions set_vol_fracs_;
};


template <MemSpace_kind MEM>
struct MeshCache : public MeshCacheBase {
  // view types
  static const MemSpace_kind MEM_kind = MEM;
  using Entity_ID_View = View_type<Entity_ID, MEM>;
  using cEntity_ID_View = View_type<const Entity_ID, MEM>;
  using Entity_GID_View = View_type<Entity_GID, MEM>;
  using cEntity_GID_View = View_type<const Entity_GID, MEM>;
  using Direction_View = View_type<Direction_type, MEM>;
  using cDirection_View = View_type<const Direction_type, MEM>;
  using Point_View = View_type<AmanziGeometry::Point, MEM>;
  using cPoint_View = View_type<const AmanziGeometry::Point, MEM>;
  using Double_View = View_type<double, MEM>;
  using cDouble_View = View_type<const double, MEM>;

  MeshCache() : MeshCacheBase() {}

  //
  // To be used by non-framework meshes
  //
  MeshCache(const Teuchos::RCP<Teuchos::ParameterList>& plist) : MeshCacheBase(plist) {}

  //
  // Standard constructor, used by the factory
  //
  MeshCache(const Teuchos::RCP<MeshFramework>& framework_mesh,
            const Teuchos::RCP<MeshAlgorithms>& algorithms,
            const Teuchos::RCP<Teuchos::ParameterList>& plist)
    : MeshCacheBase(framework_mesh, algorithms, plist)
  {
    setMeshFramework(framework_mesh);
  }

  // copy constructor
  MeshCache(const MeshCache& other) = default;

  //
  // Memory-transfer constructor -- used for getting HOST-view meshes from
  // DEVICE-view meshes, and visa versa.
  //
  template <MemSpace_kind MEM_OTHER>
  MeshCache(MeshCache<MEM_OTHER>& other);

  void setMeshFramework(const Teuchos::RCP<MeshFramework>& framework_mesh);

  Entity_GID getEntityGID(const Entity_kind kind, const Entity_ID lid) const
  {
    return getMap(kind, false)->getGlobalElement(lid);
  }
  cEntity_GID_View getEntityGIDs(const Entity_kind kind, const Entity_ID lid) const
  {
    return getMap(kind, false)->getMyGlobalIndices();
  }

  Teuchos::RCP<const MeshAlgorithms> getAlgorithms() const { return algorithms_; }

  //
  // Build the cache, fine grained control
  // =============================================

  // cell centroid, volume
  void cacheCellGeometry();
  void cacheCellFaces();
  void cacheCellEdges();
  void cacheCellNodes();
  void cacheCellCoordinates();

  // // face centroid, area, normals
  void cacheFaceGeometry();
  void cacheFaceCells();
  void cacheFaceEdges();
  void cacheFaceNodes();
  void cacheFaceCoordinates();

  // // edge centroid, length, vector
  void cacheEdgeGeometry();
  void cacheEdgeCells();
  void cacheEdgeFaces();
  void cacheEdgeNodes();
  void cacheEdgeCoordinates();

  // // node-cell adjacencies
  void cacheNodeCells();
  void cacheNodeFaces();
  void cacheNodeEdges();
  void cacheNodeCoordinates();

  // Parent entities may need to be cached too
  void cacheParentEntities();

  // // Note that regions are cached on demand the first time they are requested,
  // // but labeled sets must be pre-cached if the framework mesh is to be
  // // destroyed.
  // void precacheLabeledSets();

  //
  // Baseline mesh functionality
  // =============================================
  // ----------------------
  // Accessors and Mutators
  // ----------------------
  // Some meshes are subsets of or derived from a parent mesh.
  // Usually this is null, but some meshes may provide it.
  Teuchos::RCP<const MeshCache<MEM>> getParentMesh() const { return parent_; }
  void setParentMesh(const Teuchos::RCP<const MeshCache<MEM>>& parent);

  // Some meshes have a corresponding mesh that is better for visualization.
  const MeshCache<MEM>& getVisMesh() const
  {
    if (vis_mesh_.get()) return *vis_mesh_;
    return *this;
  }
  Teuchos::RCP<const MeshCache<MEM>> getVisMeshPtr() const { return vis_mesh_; }
  void setVisMesh(const Teuchos::RCP<const MeshCache<MEM>>& vis_mesh) { vis_mesh_ = vis_mesh; }


  // -------------------
  // Access map objects
  // -------------------
  //
  // a list of ALL face LIDs that are on the boundary
  cEntity_ID_View getBoundaryFaces() const { return maps_.getBoundaryFaces<MEM>(); }

  // a list of ALL node LIDs that are on the boundary
  cEntity_ID_View getBoundaryNodes() const { return maps_.getBoundaryNodes<MEM>(); }

  // maps define GIDs of each Entity_kind
  const Map_ptr_type& getMap(const Entity_kind kind, bool is_ghosted) const
  {
    return maps_.getMap(kind, is_ghosted);
  }

  // importers allow scatter/gather operations
  const Import_type& getImporter(const Entity_kind kind) const { return maps_.getImporter(kind); }

  // an importer from FACE-indexed objects to BOUNDARY_FACE-indexed objects
  //
  // Note this is not the same as getImporter(BOUNDARY_FACE), which
  // communicates BOUNDARY_FACE-indexed objects to other BOUNDARY_FACE-indexed
  // objects.
  const Import_type& getBoundaryFaceImporter() const { return maps_.getBoundaryFaceImporter(); }

  // an importer from NODE-indexed objects to BOUNDARY_NODE-indexed objects
  const Import_type& getBoundaryNodeImporter() const { return maps_.getBoundaryNodeImporter(); }

  // an importer from CELL-indexed objects to BOUNDARY_FACE-indexed objects
  const Import_type& getBoundaryFaceInternalCellImporter() const
  {
    return maps_.getBoundaryFaceInternalCellImporter();
  }

  // ----------------
  // sets of entities
  // ----------------
  // NOTE: no test for this -- unclear exactly the semantic of what is and is not valid.  FIXME!
  bool isValidSetName(const std::string& name, const Entity_kind kind) const { return true; }
  bool isValidSetType(const AmanziGeometry::RegionType rtype, const Entity_kind kind) const;

  int getSetSize(const std::string& region_name,
                 const Entity_kind kind,
                 const Parallel_kind ptype) const;

  cEntity_ID_View getSetEntities(const std::string& region_name,
                                 const Entity_kind kind,
                                 const Parallel_kind ptype) const;

  Kokkos::pair<cEntity_ID_View, cDouble_View>
  getSetEntitiesAndVolumeFractions(const std::string& region_name,
                                   const Entity_kind kind,
                                   const Parallel_kind ptype) const;

  // ----------------
  // Entity meta-data
  // ----------------
  Entity_ID getNumEntities(const Entity_kind kind, const Parallel_kind ptype) const;

  // // corresponding entity in the parent mesh
  // //
  // // Note the kind refers to the kind in _this_ mesh -- for some lifted meshes,
  // // this may not be the same as the entity kind in the parent mesh.  That
  // // logic is left to the user of this class -- we simply store the IDs.
  KOKKOS_INLINE_FUNCTION
  Entity_ID getEntityParent(const Entity_kind kind, const Entity_ID entid) const;

  cEntity_ID_View getEntityParents(const Entity_kind kind) const;

  KOKKOS_INLINE_FUNCTION
  Cell_kind getCellType(const Entity_ID c) const;

  KOKKOS_INLINE_FUNCTION
  Parallel_kind getParallelType(const Entity_kind& kind, const Entity_ID id) const;
  //---------------------
  // Geometry
  //---------------------
  // node locations
  template <AccessPattern_kind = AccessPattern_kind::DEFAULT>
  KOKKOS_INLINE_FUNCTION AmanziGeometry::Point getNodeCoordinate(const Entity_ID n) const;

  KOKKOS_INLINE_FUNCTION
  void setNodeCoordinate(const Entity_ID n, const AmanziGeometry::Point& coord);
  void setNodeCoordinates(const cEntity_ID_View& nodes, const cPoint_View& new_coords);

  // coordinate views
  template <AccessPattern_kind = AccessPattern_kind::DEFAULT>
  KOKKOS_INLINE_FUNCTION cPoint_View getEdgeCoordinates(const Entity_ID e) const;

  template <AccessPattern_kind = AccessPattern_kind::DEFAULT>
  KOKKOS_INLINE_FUNCTION cPoint_View getFaceCoordinates(const Entity_ID f) const;

  template <AccessPattern_kind = AccessPattern_kind::DEFAULT>
  KOKKOS_INLINE_FUNCTION cPoint_View getCellCoordinates(const Entity_ID c) const;

  // cell centroids
  template <AccessPattern_kind AP = AccessPattern_kind::DEFAULT>
  KOKKOS_INLINE_FUNCTION AmanziGeometry::Point getCellCentroid(const Entity_ID c) const;

  // face centroids
  template <AccessPattern_kind = AccessPattern_kind::DEFAULT>
  KOKKOS_INLINE_FUNCTION AmanziGeometry::Point getFaceCentroid(const Entity_ID f) const;

  template <AccessPattern_kind = AccessPattern_kind::DEFAULT>
  KOKKOS_INLINE_FUNCTION AmanziGeometry::Point getEdgeCentroid(const Entity_ID e) const;

  // at run-time, calls one of get*Centroid
  template <AccessPattern_kind = AccessPattern_kind::DEFAULT>
  KOKKOS_INLINE_FUNCTION AmanziGeometry::Point
  getCentroid(const Entity_kind kind, const Entity_ID ent) const;

  // at compile-time, calls one of get*Centroid
  template <Entity_kind, AccessPattern_kind = AccessPattern_kind::DEFAULT>
  KOKKOS_INLINE_FUNCTION AmanziGeometry::Point getCentroid(const Entity_ID ent) const;

  // extent
  template <AccessPattern_kind AP = AccessPattern_kind::DEFAULT>
  KOKKOS_INLINE_FUNCTION // double
    double
    getCellVolume(const Entity_ID c) const;

  template <AccessPattern_kind AP = AccessPattern_kind::DEFAULT>
  KOKKOS_INLINE_FUNCTION double getFaceArea(const Entity_ID f) const;

  template <AccessPattern_kind = AccessPattern_kind::DEFAULT>
  KOKKOS_INLINE_FUNCTION double getEdgeLength(const Entity_ID e) const;

  // at run-time, calls one of getCellVolume, getFaceArea, or getEdgeLength
  template <Entity_kind EK, AccessPattern_kind AP = AccessPattern_kind::DEFAULT>
  double getExtent(const Entity_ID e) const;

  // at compile-time, calls one of getCellVolume, getFaceArea, or getEdgeLength
  template <AccessPattern_kind AP = AccessPattern_kind::DEFAULT>
  double getExtent(const Entity_kind kind, const Entity_ID e) const;

  // // Normal vector of a face
  template <AccessPattern_kind = AccessPattern_kind::DEFAULT>
  KOKKOS_INLINE_FUNCTION AmanziGeometry::Point getFaceNormal(const Entity_ID f) const;

  // Normal vector and natural direction of a face, outward with respect to a
  // cell.
  //
  // The vector is normalized and then weighted by the area of the face.
  //
  // The orientation is 1 if the outward normal is the same direction as the
  // natural normal, -1 if in opposite directions, and 0 if there is no natural
  // normal.
  template <AccessPattern_kind = AccessPattern_kind::DEFAULT>
  KOKKOS_INLINE_FUNCTION AmanziGeometry::Point
  getFaceNormal(const Entity_ID f, const Entity_ID c, int* orientation = nullptr) const;

  // Vector describing the edge, where the length is the edge length.
  //
  // Orientation is the natural orientation, e.g. that it points from node 0 to
  // node 1 with respect to edge_node adjacency information.
  template <AccessPattern_kind = AccessPattern_kind::DEFAULT>
  KOKKOS_INLINE_FUNCTION AmanziGeometry::Point getEdgeVector(const Entity_ID e) const;

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
  template <AccessPattern_kind AP = AccessPattern_kind::DEFAULT>
  KOKKOS_INLINE_FUNCTION size_type getCellNumFaces(const Entity_ID c) const;

  // note, no AccessPattern_kind -- as this creates a view there is no need
  template <AccessPattern_kind AP = AccessPattern_kind::DEFAULT>
  KOKKOS_INLINE_FUNCTION cEntity_ID_View getCellFaces(const Entity_ID c) const;

  // note, no AccessPattern_kind -- only works on cached
  KOKKOS_INLINE_FUNCTION
  const Entity_ID& getCellFace(const Entity_ID c, const size_type i) const;

  KOKKOS_INLINE_FUNCTION
  Kokkos::pair<cEntity_ID_View, cDirection_View> getCellFacesAndDirections(const Entity_ID c) const;

  KOKKOS_INLINE_FUNCTION
  Kokkos::pair<cEntity_ID_View, cPoint_View> getCellFacesAndBisectors(const Entity_ID c) const;

  template <AccessPattern_kind AP = AccessPattern_kind::DEFAULT>
  KOKKOS_INLINE_FUNCTION void getCellFaces(const Entity_ID c, cEntity_ID_View& faces) const;

  KOKKOS_INLINE_FUNCTION
  void
  getCellFacesAndDirs(const Entity_ID c, cEntity_ID_View& faces, cDirection_View* const dirs) const;

  KOKKOS_INLINE_FUNCTION
  void getCellFacesAndBisectors(const Entity_ID c,
                                cEntity_ID_View& faces,
                                cPoint_View* const bisectors) const;

  // //
  // // Downward adjacency -- edges of a cell
  // //
  template <AccessPattern_kind = AccessPattern_kind::DEFAULT>
  KOKKOS_INLINE_FUNCTION size_type getCellNumEdges(const Entity_ID c) const;

  template <AccessPattern_kind AP = AccessPattern_kind::DEFAULT>
  KOKKOS_INLINE_FUNCTION cEntity_ID_View getCellEdges(const Entity_ID c) const;

  KOKKOS_INLINE_FUNCTION
  const Entity_ID& getCellEdge(const Entity_ID c, const size_type i) const;

  //[[deprecated("Prefer to use non-void variant that returns edges directly")]]
  template <AccessPattern_kind AP = AccessPattern_kind::DEFAULT>
  KOKKOS_INLINE_FUNCTION void getCellEdges(const Entity_ID c, cEntity_ID_View& edges) const;

  // Get nodes of a cell.
  template <AccessPattern_kind = AccessPattern_kind::DEFAULT>
  KOKKOS_INLINE_FUNCTION size_type getCellNumNodes(const Entity_ID c) const;

  KOKKOS_INLINE_FUNCTION
  cEntity_ID_View getCellNodes(const Entity_ID c) const;

  KOKKOS_INLINE_FUNCTION
  Entity_ID getCellNode(const Entity_ID c, const size_type i) const;

  //[[deprecated("Prefer to use non-void variant that returns nodes directly")]]
  KOKKOS_INLINE_FUNCTION
  void getCellNodes(const Entity_ID c, cEntity_ID_View& nodes) const;

  KOKKOS_INLINE_FUNCTION
  Entity_ID getCellCellBelow(const Entity_ID cellid) const;

  KOKKOS_INLINE_FUNCTION
  std::size_t getCellMaxNodes() const;

  // // Get edges of a face and directions in which the face uses the edges.
  // //
  // // In 3D, edge direction is 1 when it is oriented counter clockwise
  // // with respect to the face natural normal.
  // //
  // // On a distributed mesh, this will return all the edges of the
  // // face, OWNED or GHOST. If the framework supports it, the edges will be
  // // returned in a ccw order around the face as it is naturally defined.
  // //
  // // IMPORTANT NOTE IN 2D: In meshes where the cells are two
  // // dimensional, faces and edges are identical. For such cells, this
  // // operator will return a single edge and a direction of 1. However,
  // // this direction cannot be relied upon to compute, say, a contour
  // // integral around the 2D cell.
  // template<AccessPattern_kind = AccessPattern_kind::DEFAULT>
  // KOKKOS_INLINE_FUNCTION
  // size_type getFaceNumEdges(const Entity_ID f) const;

  KOKKOS_INLINE_FUNCTION
  cEntity_ID_View getFaceEdges(const Entity_ID f) const;

  KOKKOS_INLINE_FUNCTION
  const Entity_ID& getFaceEdge(const Entity_ID f, const size_type i) const;

  KOKKOS_INLINE_FUNCTION
  Kokkos::pair<cEntity_ID_View, cDirection_View> getFaceEdgesAndDirections(const Entity_ID f) const;

  //[[deprecated("Prefer to use non-void variant that returns edges directly")]]
  KOKKOS_INLINE_FUNCTION
  void getFaceEdges(const Entity_ID f, cEntity_ID_View& fedges) const;

  //[[deprecated("Prefer to use non-void variant that returns edges directly")]]
  KOKKOS_INLINE_FUNCTION
  void getFaceEdgesAndDirs(const Entity_ID f,
                           cEntity_ID_View& edges,
                           cDirection_View* const dirs = nullptr) const;

  KOKKOS_INLINE_FUNCTION
  std::vector<int> getFaceCellEdgeMap(const Entity_ID faceid, const Entity_ID cellid) const;

  // // Get nodes of face
  // //
  // // In 3D, the nodes of the face are returned in ccw order consistent
  // // with the face normal.
  // template<AccessPattern_kind = AccessPattern_kind::DEFAULT>
  // KOKKOS_INLINE_FUNCTION
  // size_type getFaceNumNodes(const Entity_ID f) const;

  KOKKOS_INLINE_FUNCTION
  cEntity_ID_View getFaceNodes(const Entity_ID f) const;

  KOKKOS_INLINE_FUNCTION
  const Entity_ID& getFaceNode(const Entity_ID f, const size_type i) const;

  KOKKOS_INLINE_FUNCTION
  void getFaceNodes(const Entity_ID f, cEntity_ID_View& nodes) const;

  // //
  // // NOT CURRENTLY IMPLEMENTED, here to satisfy the interface
  // //
  KOKKOS_INLINE_FUNCTION
  cPoint_View getFaceHOCoordinates(const Entity_ID f) const { return cPoint_View(); }

  // // Get nodes of edge
  // template<AccessPattern_kind = AccessPattern_kind::DEFAULT>
  // KOKKOS_INLINE_FUNCTION
  // size_type getEdgeNumNodes(const Entity_ID e) const;

  template <AccessPattern_kind AP = AccessPattern_kind::DEFAULT>
  KOKKOS_INLINE_FUNCTION cEntity_ID_View getEdgeNodes(const Entity_ID e) const;

  template <AccessPattern_kind AP = AccessPattern_kind::DEFAULT>
  KOKKOS_INLINE_FUNCTION Entity_ID getEdgeNode(const Entity_ID e, const size_type i) const;

  // //
  // // NOT CURRENTLY IMPLEMENTED, here to satisfy the interface
  // //
  // // NOTE: user code is incorrect too, and requires this to return
  // // non-const.  When or if this gets implemented, it needs to return const so
  // // that user code does not change the entries in the View!
  KOKKOS_INLINE_FUNCTION
  Point_View getEdgeHOCoordinates(const Entity_ID e) const { return Point_View(); }

  //[[deprecated("Prefer to use non-void variant that returns nodes directly")]]
  template <AccessPattern_kind AP = AccessPattern_kind::DEFAULT>
  KOKKOS_INLINE_FUNCTION void getEdgeNodes(const Entity_ID e, cEntity_ID_View& nodes) const;

  //-------------------
  // Upward adjacencies
  //-------------------
  // The cells are returned in no particular order. Also, the order of cells
  // is not guaranteed to be the same for corresponding faces on different
  // processors

  // template<AccessPattern_kind = AccessPattern_kind::DEFAULT>
  // KOKKOS_INLINE_FUNCTION
  // size_type getFaceNumCells(const Entity_ID f, const Parallel_kind ptype) const;

  template <AccessPattern_kind AP = AccessPattern_kind::DEFAULT>
  KOKKOS_INLINE_FUNCTION size_type getFaceNumCells(const Entity_ID f,
                                                   const Parallel_kind ptype) const;

  KOKKOS_INLINE_FUNCTION
  cEntity_ID_View getFaceCells(const Entity_ID f) const;

  KOKKOS_INLINE_FUNCTION
  const Entity_ID& getFaceCell(const Entity_ID f, const size_type i) const;

  template <AccessPattern_kind AP = AccessPattern_kind::DEFAULT>
  KOKKOS_INLINE_FUNCTION void getFaceCells(const Entity_ID f, cEntity_ID_View& cells) const;

  KOKKOS_INLINE_FUNCTION
  std::size_t getCellMaxFaces() const;

  KOKKOS_INLINE_FUNCTION
  std::size_t getCellMaxEdges() const;

  // Cells of a given Parallel_kind connected to an edge
  //
  // The order of cells is not guaranteed to be the same for corresponding
  // edges on different processors
  template <AccessPattern_kind AP = AccessPattern_kind::DEFAULT>
  KOKKOS_INLINE_FUNCTION cEntity_ID_View getEdgeCells(const Entity_ID e) const;


  //[[deprecated("Prefer to use non-void variant that returns edges directly")]]
  template <AccessPattern_kind AP = AccessPattern_kind::DEFAULT>
  KOKKOS_INLINE_FUNCTION void getEdgeCells(const Entity_ID e, cEntity_ID_View& cells) const;

  // Faces of type 'ptype' connected to an edge
  // NOTE: The order of faces is not guaranteed to be the same for
  // corresponding edges on different processors
  template <AccessPattern_kind AP = AccessPattern_kind::DEFAULT>
  KOKKOS_INLINE_FUNCTION cEntity_ID_View getEdgeFaces(const Entity_ID e) const;

  //[[deprecated("Prefer to use non-void variant that returns edges directly")]]
  template <AccessPattern_kind AP = AccessPattern_kind::DEFAULT>
  KOKKOS_INLINE_FUNCTION void getEdgeFaces(const Entity_ID edgeid, cEntity_ID_View& faces) const;

  // Cells of type 'ptype' connected to a node
  // NOTE: The order of cells is not guaranteed to be the same for
  // corresponding nodes on different processors
  template <AccessPattern_kind AP = AccessPattern_kind::DEFAULT>
  KOKKOS_INLINE_FUNCTION cEntity_ID_View
  getNodeCells(const Entity_ID n, const Parallel_kind ptype = Parallel_kind::OWNED) const;

  //[[deprecated("Prefer to use non-void variant that returns edges directly")]]
  template <AccessPattern_kind AP = AccessPattern_kind::DEFAULT>
  KOKKOS_INLINE_FUNCTION void
  getNodeCells(const Entity_ID n, const Parallel_kind ptype, cEntity_ID_View& cells) const;

  // Faces of type parallel 'ptype' connected to a node
  // NOTE: The order of faces is not guarnateed to be the same for
  // corresponding nodes on different processors
  template <AccessPattern_kind AP = AccessPattern_kind::DEFAULT>
  KOKKOS_INLINE_FUNCTION cEntity_ID_View getNodeFaces(const Entity_ID n) const;

  //[[deprecated("Prefer to use non-void variant that returns edges directly")]]
  template <AccessPattern_kind AP = AccessPattern_kind::DEFAULT>
  KOKKOS_INLINE_FUNCTION void getNodeFaces(const Entity_ID n, cEntity_ID_View& faces) const;

  void PrintMeshStatistics() const;

  bool isPointInCell(const AmanziGeometry::Point& p, const Entity_ID cellid) const;

  // // Edges of type 'ptype' connected to a node
  // //
  // // The order of edges is not guaranteed to be the same for corresponding
  // // node on different processors
  // KOKKOS_INLINE_FUNCTION
  // View_type<const Entity_ID,MEM> getNodeEdges(const Entity_ID n,
  //         const Parallel_kind ptype) const {
  //   Errors::Message msg("MeshCache::getNodeEdges not implemented");
  //   Exceptions::amanzi_throw(msg);
  // }


  // //[[deprecated("Prefer to use non-void variant that returns edges directly")]]
  // KOKKOS_INLINE_FUNCTION
  // void getNodeEdges(const Entity_ID n,
  //                   const Parallel_kind ptype,
  //                   cEntity_ID_View& edgeids) const {
  //   Errors::Message msg("MeshCache::getNodeEdges not implemented");
  //   Exceptions::amanzi_throw(msg);
  // }

  // column structure
  MeshColumns columns;

  inline void buildColumns() { columns.initialize(*this); }
  inline void buildColumns(const std::vector<std::string>& regions)
  {
    columns.initialize(*this, regions);
  }

  void recacheGeometry();

 protected:
  // common error messaging
  void throwAccessError_(const std::string& func_name) const;

  // these are kept here, not in base, because they are device-specific
  Teuchos::RCP<const MeshCache> parent_;
  Teuchos::RCP<const MeshCache> vis_mesh_;
};


// -----------------------------------------------------------------------------
// Caching algorithms
// -----------------------------------------------------------------------------
template <MemSpace_kind MEM>
void
cacheDefault(MeshCache<MEM>& mesh);

template <MemSpace_kind MEM>
void
cacheAll(MeshCache<MEM>& mesh);


// -----------------------------------------------------------------------------
// Memory space transfer
// -----------------------------------------------------------------------------
template <MemSpace_kind OUT, MemSpace_kind IN>
const MeshCache<OUT>
onMemSpace(const MeshCache<IN>& mc_in);

template <MemSpace_kind OUT, MemSpace_kind IN>
Teuchos::RCP<const MeshCache<OUT>>
onMemSpace(const Teuchos::RCP<const MeshCache<IN>>& mc_in);

template <MemSpace_kind OUT, MemSpace_kind IN>
Teuchos::RCP<MeshCache<OUT>>
onMemSpace(const Teuchos::RCP<MeshCache<IN>>& mc_in);


// -----------------------------------------------------------------------------
// Default Mesh types
// -----------------------------------------------------------------------------
using Mesh = MeshCache<MemSpace_kind::HOST>;
using MeshHost = MeshCache<MemSpace_kind::HOST>;

// may remove these eventually, but putting these in the AmanziMesh namespace
// for backwards compatibility now.
using Entity_ID_View = Mesh::Entity_ID_View;
using cEntity_ID_View = Mesh::cEntity_ID_View;
using Entity_GID_View = Mesh::Entity_GID_View;
using cEntity_GID_View = Mesh::cEntity_GID_View;
using Direction_View = Mesh::Direction_View;
using cDirection_View = Mesh::cDirection_View;
using Point_View = Mesh::Point_View;
using cPoint_View = Mesh::cPoint_View;
using Double_View = Mesh::Double_View;
using cDouble_View = Mesh::cDouble_View;


} // namespace AmanziMesh
} // namespace Amanzi
