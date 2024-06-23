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
   been deleted, from existing cached data_.

Effectively this set of requirements means there are two or three ways of
getting each piece of information:

1. the MeshCacheHost's cache (e.g. "fast" access)
2. recomputing using a mesh algorithm and the public interface of MeshCacheHost.
3. directly from the MeshFramework (e.g. "framework" access), which may itself
   call a mesh algorithm and its own interface.

Note that the second is not always possible -- for instance, node coordinates
cannot be recomputed; nor can basic adjacency information, as there is no
requirement here that a minimal representation of the framework has been cached
before it is destroyed.  But many of these _can_ be recomputed; for instance,
one can determine cell coordinates from node coordinates and cell-node
adjacency, which itself can be computed from cell-face and face-node adjacency
information.  Therefore, it is likely faster to recompute with MeshCacheHost, using
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
these are stored in a struct that can be shared across all MeshCacheHosts of the
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

Note that the MeshCacheHost stores pointers to many other classes that are useful
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

#include "GeometryDefs.hh"
#include "MeshDefs.hh"
#include "MeshCacheData.hh"

namespace Amanzi {

// -- forward declarations
namespace AmanziGeometry { class GeometricModel; }

namespace AmanziMesh {

class MeshFramework;
struct MeshAlgorithms;
class MeshMaps;
class MeshColumns;
struct MeshCache;

template <class Mesh_type>
class MeshAudit_Sets;
// -- end forward declarations

class Mesh {
 public:
  // view types
  static const MemSpace_kind MEM = MemSpace_kind::HOST;
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

  //
  // To be used by non-framework meshes
  //
  explicit Mesh(const Teuchos::RCP<Teuchos::ParameterList>& plist);

  //
  // Standard constructor, used by the factory
  //
  Mesh(const Teuchos::RCP<MeshFramework>& framework_mesh,
       const Teuchos::RCP<MeshAlgorithms>& algorithms,
       const Teuchos::RCP<Teuchos::ParameterList>& plist);

  // delete the copy constructor
  //  Mesh(const Mesh& other) = delete;

  //
  // Build the cache, fine grained control
  // =============================================
  void cacheCellGeometry();  // cell centroid, volume
  void cacheCellFaces();
  void cacheCellEdges();
  void cacheCellNodes();
  void cacheCellCoordinates();

  void cacheFaceGeometry(); // centroid, area, normal
  void cacheFaceCells();
  void cacheFaceEdges();
  void cacheFaceNodes();
  void cacheFaceCoordinates();

  void cacheEdgeGeometry();  // edge centroid, length, vector
  void cacheEdgeCells();
  void cacheEdgeFaces();
  void cacheEdgeNodes();
  void cacheEdgeCoordinates();

  void cacheNodeCells();
  void cacheNodeFaces();
  void cacheNodeEdges();
  void cacheNodeCoordinates();

  // Parent entities may need to be cached too
  void cacheParentEntities();

  // Note that regions are cached on demand the first time they are requested,
  // but labeled sets must be pre-cached if the framework mesh is to be
  // destroyed.
  void precacheLabeledSets();

  // lumped caching
  void cacheDefault();
  void cacheAll();
  void recacheGeometry();
  void syncCache(); // call after all other cache* functions are done!

  // Columnar semi-structured meshes
  void buildColumns();
  void buildColumns(const std::vector<std::string>& regions);

  //
  // Baseline mesh functionality
  // =============================================
  // ----------------------
  // Accessors and Mutators
  // ----------------------
  // MPI_Comm for multi-node parallelism
  Comm_ptr_type getComm() const { return comm_; }
  void setComm(const Comm_ptr_type& comm) { comm_ = comm; }

  // Geometric model describes regions
  Teuchos::RCP<const AmanziGeometry::GeometricModel> getGeometricModel() const { return gm_; }
  void setGeometricModel(const Teuchos::RCP<const AmanziGeometry::GeometricModel>& gm) { gm_ = gm; }

  // The Host-only CPU based, slow (virtual) mesh framework
  Teuchos::RCP<const MeshFramework> getMeshFramework() const { return framework_mesh_; }
  Teuchos::RCP<MeshFramework> getMeshFramework() { return framework_mesh_; }
  void destroyFramework() { framework_mesh_ = Teuchos::null; }
  void setMeshFramework(const Teuchos::RCP<MeshFramework>& framework_mesh);

  // Algorithms for computing quantities rather than caching them
  Teuchos::RCP<const MeshAlgorithms> getAlgorithms() const { return algorithms_; }
  Teuchos::RCP<Teuchos::ParameterList> getParameterList() const { return plist_; }

  // Some meshes are subsets of or derived from a parent mesh.
  // Usually this is null, but some meshes may provide it.
  Teuchos::RCP<const Mesh> getParentMesh() const { return parent_; }
  void setParentMesh(const Teuchos::RCP<const Mesh>& parent);

  // Some meshes have a corresponding mesh that is better for visualization,
  // e.g. 3D surface meshes.
  const Mesh& getVisMesh() const;
  Teuchos::RCP<const Mesh> getVisMeshPtr() const { return vis_mesh_; }
  void setVisMesh(const Teuchos::RCP<const Mesh>& vis_mesh) { vis_mesh_ = vis_mesh; }

  // get the cache for use on device
  const MeshCache& getCache() const { return *cache_; }

  // basic info written to vo
  void printMeshStatistics() const;

  // ----------------
  // mesh properties
  // ----------------
  int getSpaceDimension() const { return data_.space_dim_; }
  int getManifoldDimension() const { return data_.manifold_dim_; }
  bool isOrdered() const { return data_.is_ordered_; }
  bool isLogical() const { return data_.is_logical_; }
  bool isSFM() const { return data_.is_sfm_; } // single face mesh -- special case

  bool hasNodes() const { return data_.has_nodes_; }
  bool hasEdges() const { return data_.has_edges_; }
  bool hasNodeFaces() const { return data_.has_node_faces_; }

  // -------------------
  // Map objects
  // -------------------
  //
  // maps define GIDs of each Entity_kind
  const Map_ptr_type& getMap(const Entity_kind kind, bool is_ghosted) const;

  Entity_GID getEntityGID(const Entity_kind kind, const Entity_ID lid) const;
  cEntity_GID_View getEntityGIDs(const Entity_kind kind, bool ghosted) const;
  Entity_ID getEntityLID(const Entity_kind kind, const Entity_GID gid, bool ghosted = true) const;


  // importers allow scatter/gather operations
  const Import_type& getImporter(const Entity_kind kind) const;

  // a list of ALL face LIDs that are on the boundary
  inline cEntity_ID_View getBoundaryFaces() const;

  // given a bounday face, return the corresponding face
  inline Entity_ID getBoundaryFaceFace(const Entity_ID bf) const;

  // a list of ALL node LIDs that are on the boundary
  inline cEntity_ID_View getBoundaryNodes() const;

  // given a bounday node, return the corresponding node
  inline Entity_ID getBoundaryNodeNode(const Entity_ID bn) const;

  // an importer from FACE-indexed objects to BOUNDARY_FACE-indexed objects
  //
  // Note this is not the same as getImporter(BOUNDARY_FACE), which
  // communicates BOUNDARY_FACE-indexed objects to other BOUNDARY_FACE-indexed
  // objects.
  const Import_type& getBoundaryFaceImporter() const;

  // an importer from NODE-indexed objects to BOUNDARY_NODE-indexed objects
  const Import_type& getBoundaryNodeImporter() const;

  // an importer from CELL-indexed objects to BOUNDARY_FACE-indexed objects
  const Import_type& getBoundaryFaceInternalCellImporter() const;

  // ----------------
  // sets of entities
  // ----------------
  // Sets are only stored on the Host mesh -- to use them on device, call:
  //   mesh.getSetEntities<MemSpace_kind::DEVICE(...)
  // before launching the kernel, then capture the resulting view in the lambda.
  bool isValidSetType(const AmanziGeometry::RegionType rtype, const Entity_kind kind) const;

  // NOTE: unclear exactly the semantic of what is and is not valid.  FIXME!
  bool isValidSetName(const std::string& name, const Entity_kind kind) const { return true; }

  int getSetSize(const std::string& region_name,
                 const Entity_kind kind,
                 const Parallel_kind ptype) const;

  template<MemSpace_kind MEM>
  auto // cEntity_ID_View on MEM
  getSetEntities(const std::string& region_name,
                 const Entity_kind kind,
                 const Parallel_kind ptype) const;

  template<MemSpace_kind MEM>
  auto // Kokkos::pair<cEntity_ID_View, cDouble_View> on MEM
  getSetEntitiesAndVolumeFractions(const std::string& region_name,
                                   const Entity_kind kind,
                                   const Parallel_kind ptype) const;

  // ----------------
  // Entity meta-data
  // ----------------
  Entity_ID getNumEntities(const Entity_kind kind, const Parallel_kind ptype) const;

  // corresponding entity in the parent mesh
  //
  // Note the kind refers to the kind in _this_ mesh -- for some lifted meshes,
  // this may not be the same as the entity kind in the parent mesh.  That
  // logic is left to the user of this class -- we simply store the IDs.
  Entity_ID getEntityParent(const Entity_kind kind, const Entity_ID entid) const;
  cEntity_ID_View getEntityParents(const Entity_kind kind) const;

  Cell_kind getCellKind(const Entity_ID c) const;
  Parallel_kind getParallelKind(const Entity_kind& kind, const Entity_ID id) const;

  //-----------------------
  // Geometry: coordinates
  //-----------------------
  AmanziGeometry::Point getNodeCoordinate(const Entity_ID n) const;

  void setNodeCoordinate(const Entity_ID n, const AmanziGeometry::Point& coord);
  void setNodeCoordinates(const cEntity_ID_View& nodes, const cPoint_View& new_coords);

  // coordinate views
  template <AccessPattern_kind = AccessPattern_kind::DEFAULT>
  cPoint_View getCellCoordinates(const Entity_ID c) const;

  template <AccessPattern_kind = AccessPattern_kind::DEFAULT>
  cPoint_View getFaceCoordinates(const Entity_ID f) const;

  template <AccessPattern_kind = AccessPattern_kind::DEFAULT>
  cPoint_View getEdgeCoordinates(const Entity_ID e) const;

  //-----------------------
  // Geometry: centroids
  //-----------------------
  template <AccessPattern_kind AP = AccessPattern_kind::DEFAULT>
  AmanziGeometry::Point getCellCentroid(const Entity_ID c) const;

  template <AccessPattern_kind = AccessPattern_kind::DEFAULT>
  AmanziGeometry::Point getFaceCentroid(const Entity_ID f) const;

  template <AccessPattern_kind = AccessPattern_kind::DEFAULT>
  AmanziGeometry::Point getEdgeCentroid(const Entity_ID e) const;

  // at run-time, calls one of get*Centroid
  template <AccessPattern_kind = AccessPattern_kind::DEFAULT>
  AmanziGeometry::Point
  getCentroid(const Entity_kind kind, const Entity_ID ent) const;

  // at compile-time, calls one of get*Centroid
  template <Entity_kind, AccessPattern_kind = AccessPattern_kind::DEFAULT>
  AmanziGeometry::Point getCentroid(const Entity_ID ent) const;

  //-----------------------
  // Geometry: extents
  //-----------------------
  template <AccessPattern_kind AP = AccessPattern_kind::DEFAULT>
  double getCellVolume(const Entity_ID c) const;

  template <AccessPattern_kind AP = AccessPattern_kind::DEFAULT>
  double getFaceArea(const Entity_ID f) const;

  template <AccessPattern_kind = AccessPattern_kind::DEFAULT>
  double getEdgeLength(const Entity_ID e) const;

  // at run-time, calls one of getCellVolume, getFaceArea, or getEdgeLength
  template <Entity_kind EK, AccessPattern_kind AP = AccessPattern_kind::DEFAULT>
  double getExtent(const Entity_ID e) const;

  // at compile-time, calls one of getCellVolume, getFaceArea, or getEdgeLength
  template <AccessPattern_kind AP = AccessPattern_kind::DEFAULT>
  double getExtent(const Entity_kind kind, const Entity_ID e) const;

  //-----------------------
  // Geometry: other
  //-----------------------
  // Normal vector of a face
  template <AccessPattern_kind = AccessPattern_kind::DEFAULT>
  AmanziGeometry::Point getFaceNormal(const Entity_ID f) const;

  // Normal vector and natural direction of a face, outward with respect to a
  // cell.
  //
  // The vector is normalized and then weighted by the area of the face.
  //
  // The orientation is 1 if the outward normal is the same direction as the
  // natural normal, -1 if in opposite directions, and 0 if there is no natural
  // normal.
  template <AccessPattern_kind = AccessPattern_kind::DEFAULT>
  AmanziGeometry::Point
  getFaceNormal(const Entity_ID f, const Entity_ID c, int* orientation = nullptr) const;

  // Vector describing the edge, where the length is the edge length.
  //
  // Orientation is the natural orientation, e.g. that it points from node 0 to
  // node 1 with respect to edge_node adjacency information.
  template <AccessPattern_kind = AccessPattern_kind::DEFAULT>
  AmanziGeometry::Point getEdgeVector(const Entity_ID e) const;

  bool isPointInCell(const AmanziGeometry::Point& p, const Entity_ID cellid) const;

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
  size_type getCellNumFaces(const Entity_ID c) const;


  // note, no AccessPattern_kind -- only works on cached
  const Entity_ID& getCellFace(const Entity_ID c, const size_type i) const;

  cEntity_ID_View getCellFaces(const Entity_ID c) const;
  void getCellFaces(const Entity_ID c, cEntity_ID_View& cfaces) const;

  Kokkos::pair<cEntity_ID_View, cDirection_View> getCellFacesAndDirections(const Entity_ID c) const;
  void getCellFacesAndDirs(const Entity_ID c,
                           cEntity_ID_View& faces,
                           cDirection_View* const dirs) const;

  Kokkos::pair<cEntity_ID_View, cPoint_View> getCellFacesAndBisectors(const Entity_ID c) const;
  void getCellFacesAndBisectors(const Entity_ID c, cEntity_ID_View& cfaces,
          cPoint_View* const bisectors) const;


  //
  // Downward adjacency
  //
  // Get edges of a cell
  size_type getCellNumEdges(const Entity_ID c) const;

  const Entity_ID& getCellEdge(const Entity_ID c, const size_type i) const;

  template <AccessPattern_kind AP = AccessPattern_kind::DEFAULT>
  cEntity_ID_View getCellEdges(const Entity_ID c) const;

  // Get nodes of a cell.
  size_type getCellNumNodes(const Entity_ID c) const;

  const Entity_ID& getCellNode(const Entity_ID c, const size_type i) const;

  template <AccessPattern_kind = AccessPattern_kind::DEFAULT>
  cEntity_ID_View getCellNodes(const Entity_ID c) const;

  template <AccessPattern_kind = AccessPattern_kind::DEFAULT>
  void getCellNodes(const Entity_ID c, cEntity_ID_View& nodes) const;

  std::size_t getCellMaxNodes() const;

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
  size_type getFaceNumEdges(const Entity_ID f) const;

  const Entity_ID& getFaceEdge(const Entity_ID f, const size_type i) const;

  cEntity_ID_View getFaceEdges(const Entity_ID f) const;

  Kokkos::pair<cEntity_ID_View, cDirection_View>
  getFaceEdgesAndDirections(const Entity_ID f) const;

  void getFaceEdgesAndDirs(const Entity_ID f,
          cEntity_ID_View& edges,
          cDirection_View* const dirs) const;

  // Get nodes of face
  //
  // In 3D, the nodes of the face are returned in ccw order consistent
  // with the face normal.
  size_type getFaceNumNodes(const Entity_ID f) const;

  const Entity_ID& getFaceNode(const Entity_ID f, const size_type i) const;

  cEntity_ID_View getFaceNodes(const Entity_ID f) const;

  // Get nodes of edge
  size_type getEdgeNumNodes(const Entity_ID e) const;

  const Entity_ID& getEdgeNode(const Entity_ID e, const size_type i) const;

  cEntity_ID_View getEdgeNodes(const Entity_ID e) const;
  void getEdgeNodes(const Entity_ID e, cEntity_ID_View& nodes) const;

  // max sizes for discretizations
  std::size_t getCellMaxFaces() const;
  std::size_t getCellMaxEdges() const;

  //-------------------
  // Upward adjacencies
  //-------------------
  // Face --> Cell
  // The cells are returned in no particular order. Also, the order of cells
  // is not guaranteed to be the same for corresponding faces on different
  // processors
  size_type getFaceNumCells(const Entity_ID f) const;

  const Entity_ID& getFaceCell(const Entity_ID f, const size_type i) const;

  cEntity_ID_View getFaceCells(const Entity_ID f) const;
  void getFaceCells(const Entity_ID f, cEntity_ID_View& cells) const;

  // Edge --> Cell
  //
  // The order of cells is not guaranteed to be the same for corresponding
  // edges on different processors
  size_type getEdgeNumCells(const Entity_ID f) const;

  const Entity_ID& getEdgeCell(const Entity_ID f, const size_type i) const;

  template <AccessPattern_kind AP = AccessPattern_kind::DEFAULT>
  cEntity_ID_View getEdgeCells(const Entity_ID e) const;

  // Edge --> Face
  // NOTE: The order of faces is not guaranteed to be the same for
  // corresponding edges on different processors
  size_type getEdgeNumFaces(const Entity_ID f) const;

  const Entity_ID& getEdgeFace(const Entity_ID f, const size_type i) const;

  cEntity_ID_View getEdgeFaces(const Entity_ID e) const;

  // Node --> Cell
  // NOTE: The order of cells is not guaranteed to be the same for
  // corresponding nodes on different processors
  template <AccessPattern_kind AP = AccessPattern_kind::DEFAULT>
  cEntity_ID_View
  getNodeCells(const Entity_ID n) const;

  // Node --> Face
  // Faces of type parallel 'ptype' connected to a node
  // NOTE: The order of faces is not guarnateed to be the same for
  // corresponding nodes on different processors
  cEntity_ID_View getNodeFaces(const Entity_ID n) const;

  // Node --> Edge
  //
  // The order of edges is not guaranteed to be the same for corresponding
  // node on different processors
  cEntity_ID_View getNodeEdges(const Entity_ID n) const;

 protected:
  // common error messaging
  void throwAccessError_(const std::string& func_name) const;

 public:
  Teuchos::RCP<MeshColumns> columns;

 protected:
  Comm_ptr_type comm_;
  Teuchos::RCP<const AmanziGeometry::GeometricModel> gm_;
  Teuchos::RCP<Teuchos::ParameterList> plist_;

  Teuchos::RCP<MeshFramework> framework_mesh_;
  Teuchos::RCP<MeshAlgorithms> algorithms_;
  Teuchos::RCP<MeshMaps> maps_;
  Teuchos::RCP<MeshSets> sets_;
  Teuchos::RCP<MeshSetVolumeFractions> set_vol_fracs_;

  Teuchos::RCP<const Mesh> parent_;
  Teuchos::RCP<const Mesh> vis_mesh_;

  MeshCacheData data_;
  Teuchos::RCP<MeshCache> cache_;

  // friend MeshAudit_Sets so it can check any existing MeshSets
  // friend MeshAudit_Sets<Mesh>;
};


} // namespace AmanziMesh
} // namespace Amanzi
