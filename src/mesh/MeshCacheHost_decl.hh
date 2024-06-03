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
#include "GeometricModel.hh"

#include "MeshDefs.hh"
#include "MeshSets.hh"
#include "MeshColumns.hh"
#include "MeshCacheBase.hh"

namespace Amanzi {
namespace AmanziMesh {

class MeshFramework;
struct MeshAlgorithms;
struct MeshCacheDevice; 
  
// forward declaration of MeshAudit to friend
template <class Mesh_type> class MeshAudit_Sets;

struct MeshCacheHost : public MeshCacheBase {
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

  MeshCacheHost() : MeshCacheBase() {}

  //
  // To be used by non-framework meshes
  //
  MeshCacheHost(const Teuchos::RCP<Teuchos::ParameterList>& plist) : MeshCacheBase(plist) {}

  //
  // Standard constructor, used by the factory
  //
  MeshCacheHost(const Teuchos::RCP<MeshFramework>& framework_mesh,
            const Teuchos::RCP<MeshAlgorithms>& algorithms,
            const Teuchos::RCP<Teuchos::ParameterList>& plist)
    : MeshCacheBase(framework_mesh, algorithms, plist)
  {
    setMeshFramework(framework_mesh);
  }

  // copy constructor
  MeshCacheHost(const MeshCacheHost& other) = default;

  //
  // Memory-transfer constructor -- used for getting HOST-view meshes from
  // DEVICE-view meshes, and visa versa.
  //
  MeshCacheHost(const MeshCacheDevice& other);

  void setMeshFramework(const Teuchos::RCP<MeshFramework>& framework_mesh);
  
  Entity_GID getEntityGID(const Entity_kind kind, const Entity_ID lid) const
  {
    return getMap(kind, true)->getGlobalElement(lid);
  }
  cEntity_GID_View getEntityGIDs(const Entity_kind kind, bool ghosted) const
  {
    return getMap(kind, ghosted)->getMyGlobalIndices();
  }
  Entity_ID getEntityLID(const Entity_kind kind, const Entity_GID gid, bool ghosted=true) const
  {
    return getMap(kind, ghosted)->getLocalElement(gid);
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
  Teuchos::RCP<const MeshCacheHost> getParentMesh() const { return parent_; }
  void setParentMesh(const Teuchos::RCP<const MeshCacheHost>& parent);

  // Some meshes have a corresponding mesh that is better for visualization.
  const MeshCacheHost& getVisMesh() const
  {
    if (vis_mesh_.get()) return *vis_mesh_;
    return *this;
  }
  Teuchos::RCP<const MeshCacheHost> getVisMeshPtr() const { return vis_mesh_; }
  void setVisMesh(const Teuchos::RCP<const MeshCacheHost>& vis_mesh) { vis_mesh_ = vis_mesh; }


  // -------------------
  // Access map objects
  // -------------------
  //
  // a list of ALL face LIDs that are on the boundary
  cEntity_ID_View getBoundaryFaces() const  { return maps_.getBoundaryFaces<MEM>(); }

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
  
  
  Entity_ID getEntityParent(const Entity_kind kind, const Entity_ID entid) const;

  cEntity_ID_View getEntityParents(const Entity_kind kind) const;

  
  
  Cell_kind getCellType(const Entity_ID c) const;

  
  Parallel_kind getParallelType(const Entity_kind& kind, const Entity_ID id) const;
  //---------------------
  // Geometry
  //---------------------
  // node locations
  template <AccessPattern_kind = AccessPattern_kind::DEFAULT>
   AmanziGeometry::Point getNodeCoordinate(const Entity_ID n) const;

  
  void setNodeCoordinate(const Entity_ID n, const AmanziGeometry::Point& coord);
  void setNodeCoordinates(const cEntity_ID_View& nodes, const cPoint_View& new_coords);

  // coordinate views
  template <AccessPattern_kind = AccessPattern_kind::DEFAULT>
   cPoint_View getEdgeCoordinates(const Entity_ID e) const;

  template <AccessPattern_kind = AccessPattern_kind::DEFAULT>
   cPoint_View getFaceCoordinates(const Entity_ID f) const;

  template <AccessPattern_kind = AccessPattern_kind::DEFAULT>
   cPoint_View getCellCoordinates(const Entity_ID c) const;

  // cell centroids
  template <AccessPattern_kind AP = AccessPattern_kind::DEFAULT>
   AmanziGeometry::Point getCellCentroid(const Entity_ID c) const;

  // face centroids
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

  // extent
  template <AccessPattern_kind AP = AccessPattern_kind::DEFAULT>
   // double
    double
    getCellVolume(const Entity_ID c) const;

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

  // // Normal vector of a face  
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
   size_type getCellNumFaces(const Entity_ID c) const;

  // note, no AccessPattern_kind -- as this creates a view there is no need
  template <AccessPattern_kind AP = AccessPattern_kind::DEFAULT>
   cEntity_ID_View getCellFaces(const Entity_ID c) const;

  // note, no AccessPattern_kind -- only works on cached
  
  const Entity_ID& getCellFace(const Entity_ID c, const size_type i) const;

  
  Kokkos::pair<cEntity_ID_View, cDirection_View> getCellFacesAndDirections(const Entity_ID c) const;

  
  
  Kokkos::pair<cEntity_ID_View, cPoint_View> getCellFacesAndBisectors(const Entity_ID c) const;

  template <AccessPattern_kind AP = AccessPattern_kind::DEFAULT>
   void getCellFaces(const Entity_ID c, cEntity_ID_View& faces) const;

  
  void
  getCellFacesAndDirs(const Entity_ID c, cEntity_ID_View& faces, cDirection_View* const dirs) const;

  
  void getCellFacesAndBisectors(const Entity_ID c,
                                cEntity_ID_View& faces,
                                cPoint_View* const bisectors) const;

  // //
  // // Downward adjacency -- edges of a cell
  // //
  template <AccessPattern_kind = AccessPattern_kind::DEFAULT>
   size_type getCellNumEdges(const Entity_ID c) const;

  template <AccessPattern_kind AP = AccessPattern_kind::DEFAULT>
   cEntity_ID_View getCellEdges(const Entity_ID c) const;

  
  const Entity_ID& getCellEdge(const Entity_ID c, const size_type i) const;

  //[[deprecated("Prefer to use non-void variant that returns edges directly")]]
  template <AccessPattern_kind AP = AccessPattern_kind::DEFAULT>
   void getCellEdges(const Entity_ID c, cEntity_ID_View& edges) const;

  // Get nodes of a cell.
  template <AccessPattern_kind = AccessPattern_kind::DEFAULT>
   size_type getCellNumNodes(const Entity_ID c) const;

  
  cEntity_ID_View getCellNodes(const Entity_ID c) const;

  
  Entity_ID getCellNode(const Entity_ID c, const size_type i) const;

  //[[deprecated("Prefer to use non-void variant that returns nodes directly")]]
  
  void getCellNodes(const Entity_ID c, cEntity_ID_View& nodes) const;

  
  Entity_ID getCellCellBelow(const Entity_ID cellid) const;

  
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
  // 
  // size_type getFaceNumEdges(const Entity_ID f) const;

  
  cEntity_ID_View getFaceEdges(const Entity_ID f) const;

  
  const Entity_ID& getFaceEdge(const Entity_ID f, const size_type i) const;

  
  
  Kokkos::pair<cEntity_ID_View, cDirection_View> getFaceEdgesAndDirections(const Entity_ID f) const;

  //[[deprecated("Prefer to use non-void variant that returns edges directly")]]
  
  void getFaceEdges(const Entity_ID f, cEntity_ID_View& fedges) const;

  //[[deprecated("Prefer to use non-void variant that returns edges directly")]]
  
  void getFaceEdgesAndDirs(const Entity_ID f,
                           cEntity_ID_View& edges,
                           cDirection_View* const dirs = nullptr) const;

  
  std::vector<int> getFaceCellEdgeMap(const Entity_ID faceid, const Entity_ID cellid) const;

  // // Get nodes of face
  // //
  // // In 3D, the nodes of the face are returned in ccw order consistent
  // // with the face normal.
  // template<AccessPattern_kind = AccessPattern_kind::DEFAULT>
  // 
  // size_type getFaceNumNodes(const Entity_ID f) const;

  
  cEntity_ID_View getFaceNodes(const Entity_ID f) const;

  
  const Entity_ID& getFaceNode(const Entity_ID f, const size_type i) const;

  
  void getFaceNodes(const Entity_ID f, cEntity_ID_View& nodes) const;

  // //
  // // NOT CURRENTLY IMPLEMENTED, here to satisfy the interface
  // //
  
  cPoint_View getFaceHOCoordinates(const Entity_ID f) const { return cPoint_View(); }

  // // Get nodes of edge
  // template<AccessPattern_kind = AccessPattern_kind::DEFAULT>
  // 
  // size_type getEdgeNumNodes(const Entity_ID e) const;

  template <AccessPattern_kind AP = AccessPattern_kind::DEFAULT>
   cEntity_ID_View getEdgeNodes(const Entity_ID e) const;

  template <AccessPattern_kind AP = AccessPattern_kind::DEFAULT>
   Entity_ID getEdgeNode(const Entity_ID e, const size_type i) const;

  // //
  // // NOT CURRENTLY IMPLEMENTED, here to satisfy the interface
  // //
  // // NOTE: user code is incorrect too, and requires this to return
  // // non-const.  When or if this gets implemented, it needs to return const so
  // // that user code does not change the entries in the View!
  
  Point_View getEdgeHOCoordinates(const Entity_ID e) const { return Point_View(); }

  //[[deprecated("Prefer to use non-void variant that returns nodes directly")]]
  template <AccessPattern_kind AP = AccessPattern_kind::DEFAULT>
   void getEdgeNodes(const Entity_ID e, cEntity_ID_View& nodes) const;

  //-------------------
  // Upward adjacencies
  //-------------------
  // The cells are returned in no particular order. Also, the order of cells
  // is not guaranteed to be the same for corresponding faces on different
  // processors

  // template<AccessPattern_kind = AccessPattern_kind::DEFAULT>
  // 
  // size_type getFaceNumCells(const Entity_ID f, const Parallel_kind ptype) const;

  template <AccessPattern_kind AP = AccessPattern_kind::DEFAULT>
   size_type getFaceNumCells(const Entity_ID f,
                                                   const Parallel_kind ptype) const;

  template <AccessPattern_kind AP = AccessPattern_kind::DEFAULT>
  cEntity_ID_View getFaceCells(const Entity_ID f) const;

  const Entity_ID& getFaceCell(const Entity_ID f, const size_type i) const;

  template <AccessPattern_kind AP = AccessPattern_kind::DEFAULT>
   void getFaceCells(const Entity_ID f, cEntity_ID_View& cells) const;

  
  std::size_t getCellMaxFaces() const;

  
  std::size_t getCellMaxEdges() const;

  // Cells of a given Parallel_kind connected to an edge
  //
  // The order of cells is not guaranteed to be the same for corresponding
  // edges on different processors
  template <AccessPattern_kind AP = AccessPattern_kind::DEFAULT>
   cEntity_ID_View getEdgeCells(const Entity_ID e) const;


  //[[deprecated("Prefer to use non-void variant that returns edges directly")]]
  template <AccessPattern_kind AP = AccessPattern_kind::DEFAULT>
   void getEdgeCells(const Entity_ID e, cEntity_ID_View& cells) const;

  // Faces of type 'ptype' connected to an edge
  // NOTE: The order of faces is not guaranteed to be the same for
  // corresponding edges on different processors
  
  template <AccessPattern_kind AP = AccessPattern_kind::DEFAULT>
   cEntity_ID_View getEdgeFaces(const Entity_ID e) const;

  //[[deprecated("Prefer to use non-void variant that returns edges directly")]]
  template <AccessPattern_kind AP = AccessPattern_kind::DEFAULT>
   void getEdgeFaces(const Entity_ID edgeid, cEntity_ID_View& faces) const;

  // Cells of type 'ptype' connected to a node
  // NOTE: The order of cells is not guaranteed to be the same for
  // corresponding nodes on different processors
  template <AccessPattern_kind AP = AccessPattern_kind::DEFAULT>
   cEntity_ID_View
  getNodeCells(const Entity_ID n, const Parallel_kind ptype = Parallel_kind::OWNED) const;

  //[[deprecated("Prefer to use non-void variant that returns edges directly")]]
  template <AccessPattern_kind AP = AccessPattern_kind::DEFAULT>
   void
  getNodeCells(const Entity_ID n, const Parallel_kind ptype, cEntity_ID_View& cells) const;

  // Faces of type parallel 'ptype' connected to a node
  // NOTE: The order of faces is not guarnateed to be the same for
  // corresponding nodes on different processors
  
  template <AccessPattern_kind AP = AccessPattern_kind::DEFAULT>
   cEntity_ID_View getNodeFaces(const Entity_ID n) const;

  //[[deprecated("Prefer to use non-void variant that returns edges directly")]]
  
  template <AccessPattern_kind AP = AccessPattern_kind::DEFAULT>
   void getNodeFaces(const Entity_ID n, cEntity_ID_View& faces) const;

  void PrintMeshStatistics() const;

  bool isPointInCell(const AmanziGeometry::Point& p, const Entity_ID cellid) const;

  // // Edges of type 'ptype' connected to a node
  // //
  // // The order of edges is not guaranteed to be the same for corresponding
  // // node on different processors
  // 
  // View_type<const Entity_ID,MEM> getNodeEdges(const Entity_ID n,
  //         const Parallel_kind ptype) const {
  //   Errors::Message msg("MeshCacheHost::getNodeEdges not implemented");
  //   Exceptions::amanzi_throw(msg);
  // }


  // //[[deprecated("Prefer to use non-void variant that returns edges directly")]]
  // 
  // void getNodeEdges(const Entity_ID n,
  //                   const Parallel_kind ptype,
  //                   cEntity_ID_View& edgeids) const {
  //   Errors::Message msg("MeshCacheHost::getNodeEdges not implemented");
  //   Exceptions::amanzi_throw(msg);
  // }
  inline void buildColumns() {
    columns.initialize(*this);
  }
  inline void buildColumns(const std::vector<std::string>& regions)
  {
    columns.initialize(*this, regions);
  }

  void recacheGeometry();

 protected:
  // common error messaging
  void throwAccessError_(const std::string& func_name) const;

  // these are kept here, not in base, because they are device-specific
  Teuchos::RCP<const MeshCacheHost> parent_;
  Teuchos::RCP<const MeshCacheHost> vis_mesh_;

  // friend MeshAudit_Sets so it can check any existing MeshSets
  friend MeshAudit_Sets<MeshCacheHost>;
};


// -----------------------------------------------------------------------------
// Caching algorithms
// -----------------------------------------------------------------------------
inline void
cacheDefault(MeshCacheHost& mesh)
{
  if (mesh.hasNodes()) { mesh.cacheNodeCoordinates(); }
  mesh.cacheCellFaces();
  mesh.cacheFaceCells();
  if (mesh.hasNodes()) { mesh.cacheFaceNodes(); }
  mesh.cacheCellGeometry();
  mesh.cacheFaceGeometry();
  if (mesh.hasEdges()) {
    mesh.cacheFaceEdges();
    mesh.cacheEdgeFaces();
    mesh.cacheEdgeGeometry();
  }
}


inline void
cacheAll(MeshCacheHost& mesh)
{
  mesh.cacheCellFaces();
  mesh.cacheFaceCells();
  mesh.cacheCellCoordinates();
  mesh.cacheFaceCoordinates();
  if (mesh.hasNodes()) {
    mesh.cacheNodeCoordinates(); 
    mesh.cacheCellNodes();
    mesh.cacheNodeCells();
    mesh.cacheFaceNodes();
    mesh.cacheNodeFaces();
  }
  if (mesh.hasEdges()) {
    mesh.cacheCellEdges();
    mesh.cacheEdgeCells();
    mesh.cacheNodeEdges();
    mesh.cacheEdgeNodes();
    mesh.cacheEdgeCoordinates();
  }
}

using MeshHost = Amanzi::AmanziMesh::MeshCacheHost;
  
} // namespace AmanziMesh
} // namespace Amanzi
