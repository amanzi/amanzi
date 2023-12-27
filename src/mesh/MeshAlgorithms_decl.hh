/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/
// Functions for commonly used Mesh operations and algorithms
/*

  Note that these are distinct from Mesh_Helpers.  Helpers are functions most
  commonly used within the mesh library, while these are functions that are
  commonly used by clients of the mesh library.

*/

#pragma once

#include "MeshCache_decl.hh"
#include "MeshHelpers.hh"

namespace Amanzi {
namespace AmanziMesh {

// -----------------------------------------------------------------------------
// Derived topological relationships
// -----------------------------------------------------------------------------

//
// Given a boundary face ID, get the corresponding face ID
//
template <class Mesh_type>
KOKKOS_INLINE_FUNCTION Entity_ID
getBoundaryFaceFace(const Mesh_type& mesh, Entity_ID bf);


//
// Given a face ID, get the corresponding boundary face ID (assuming it is a bf)
//
template <MemSpace_kind MEM>
Entity_ID
getFaceOnBoundaryBoundaryFace(const MeshCache<MEM>& mesh, Entity_ID f);


//
// Given a boundary face ID, get the cell internal to that face.
//
template <class Mesh_type>
KOKKOS_INLINE_FUNCTION Entity_ID
getBoundaryFaceInternalCell(const Mesh_type& mesh, Entity_ID bf);


//
// Given a face ID, and assuming it is a boundary face, get the cell internal.
//
template <class Mesh_type>
KOKKOS_INLINE_FUNCTION Entity_ID
getFaceOnBoundaryInternalCell(const Mesh_type& mesh, Entity_ID f);

//
// Given a boundary face f, return the exterior normal. If f is an interior face,
// dir = 0 and normal orientation is not be reliable in parallel algorithms
//
template <class Mesh_type>
KOKKOS_INLINE_FUNCTION AmanziGeometry::Point
getFaceNormalExterior(const Mesh_type& mesh, int f, int* dir);


// -----------------------------------------------------------------------------
// Derived geometric relationships
// -----------------------------------------------------------------------------

//
// This class provides default, virtual algorithms for computing geometric
// quantities given nodal coordinates and topological information.
//
struct MeshAlgorithms {
  virtual ~MeshAlgorithms() = default;

  // lumped things for more efficient calculation
  // KOKKOS_INLINE_FUNCTION
  // virtual std::pair<double, AmanziGeometry::Point>
  // computeCellGeometry(const Mesh& mesh, const Entity_ID c) const;

  virtual inline std::pair<double, AmanziGeometry::Point>
  computeCellGeometry(const MeshHost& mesh, const Entity_ID c) const;

  // KOKKOS_INLINE_FUNCTION
  // virtual std::tuple<double, AmanziGeometry::Point, typename Mesh::cPoint_View>
  // computeFaceGeometry(const Mesh& mesh, const Entity_ID f) const;

  virtual inline std::tuple<double, AmanziGeometry::Point, typename Mesh::cPoint_View>
  computeFaceGeometry(const MeshHost& mesh, const Entity_ID f) const;

  // KOKKOS_INLINE_FUNCTION
  // virtual std::pair<AmanziGeometry::Point, AmanziGeometry::Point>
  // computeEdgeGeometry(const Mesh& mesh, const Entity_ID e) const;

  virtual inline std::pair<AmanziGeometry::Point, AmanziGeometry::Point>
  computeEdgeGeometry(const MeshHost& mesh, const Entity_ID e) const;

  // KOKKOS_INLINE_FUNCTION
  // virtual double computeCellVolume(const Mesh& mesh, const Entity_ID c) const;

  virtual inline double computeCellVolume(const MeshHost& mesh, const Entity_ID c) const;

  // KOKKOS_INLINE_FUNCTION
  // virtual AmanziGeometry::Point computeCellCentroid(const Mesh& mesh, const Entity_ID c) const;

  virtual inline AmanziGeometry::Point
  computeCellCentroid(const MeshHost& mesh, const Entity_ID c) const;

  // KOKKOS_INLINE_FUNCTION
  // virtual double computeFaceArea(const Mesh& mesh, const Entity_ID f) const;

  virtual inline double computeFaceArea(const MeshHost& mesh, const Entity_ID f) const;

  // KOKKOS_INLINE_FUNCTION
  // virtual AmanziGeometry::Point computeFaceCentroid(const Mesh& mesh, const Entity_ID f) const;

  virtual inline AmanziGeometry::Point
  computeFaceCentroid(const MeshHost& mesh, const Entity_ID f) const;

  // KOKKOS_INLINE_FUNCTION
  // virtual AmanziGeometry::Point computeFaceNormal(const Mesh& mesh, const Entity_ID f,
  //         const Entity_ID c, int * const orientation) const;

  virtual inline AmanziGeometry::Point computeFaceNormal(const MeshHost& mesh,
                                                         const Entity_ID f,
                                                         const Entity_ID c,
                                                         int* const orientation) const;

  // KOKKOS_INLINE_FUNCTION
  // virtual double computeEdgeLength(const Mesh& mesh, const Entity_ID e) const;

  virtual inline double computeEdgeLength(const MeshHost& mesh, const Entity_ID e) const;

  // KOKKOS_INLINE_FUNCTION
  // virtual AmanziGeometry::Point
  // computeEdgeVector(const Mesh& mesh, const Entity_ID e, const Entity_ID n, int * const orientation) const;

  virtual inline AmanziGeometry::Point computeEdgeVector(const MeshHost& mesh,
                                                         const Entity_ID e,
                                                         const Entity_ID n,
                                                         int* const orientation) const;

  // KOKKOS_INLINE_FUNCTION
  // virtual AmanziGeometry::Point
  // computeEdgeCentroid(const Mesh& mesh, const Entity_ID e) const;

  virtual inline AmanziGeometry::Point
  computeEdgeCentroid(const MeshHost& mesh, const Entity_ID e) const;

  // KOKKOS_INLINE_FUNCTION
  // virtual void computeCellFacesAndBisectors(const Mesh& mesh, const Entity_ID cellid,
  //         typename Mesh::cEntity_ID_View& faceids, typename Mesh::cPoint_View * const bisectors) const;

  virtual inline void
  computeCellFacesAndBisectors(const MeshHost& mesh,
                               const Entity_ID cellid,
                               typename MeshHost::cEntity_ID_View& faceids,
                               typename MeshHost::cPoint_View* const bisectors) const;
};


// -----------------------------------------------------------------------------
// Collective operations
// -----------------------------------------------------------------------------

//
// Given a vector on faces, import to vector on boundary faces
//
template <MemSpace_kind MEM>
void
copyFacesToBoundaryFaces(const MeshCache<MEM>& mesh,
                         const Epetra_MultiVector& faces,
                         Epetra_MultiVector& boundary_faces);

//
// Given a vector on faces, import to vector on boundary faces
//
template <MemSpace_kind MEM>
void
copyBoundaryFacesToFaces(const MeshCache<MEM>& mesh,
                         const Epetra_MultiVector& boundary_faces,
                         Epetra_MultiVector& faces);

//
// Given a vector on cells, set the boundary_face entries by their internal cell
//
template <MemSpace_kind MEM>
void
copyCellsToBoundaryFaces(const MeshCache<MEM>& mesh,
                         const Epetra_MultiVector& cells,
                         Epetra_MultiVector& boundary_faces);


// -----------------------------------------------------------------------------
// Map manipulation and creation
// -----------------------------------------------------------------------------

//
// Creates a pair of maps, <ALL, OWNED>, for a given entity_kind.
// Uses GIDs provided by the Mesh object.
//
template <class Mesh_type>
std::pair<Map_ptr_type, Map_ptr_type>
createMapsFromMeshGIDs(const Mesh_type& mesh, const Entity_kind kind);

//
// Given discontinuous maps, make continuous maps
//
inline std::pair<Map_ptr_type, Map_ptr_type>
createContiguousMaps(const Map_ptr_type& ghosted, const Map_ptr_type& owned);

//
// Creates a pair of maps, <ALL, OWNED>, for a given entity_kind.
// Uses a natural ordering of GIDs, proc 0 == 0...n, proc 1 = n..., etc.
//
template <class Mesh_type>
std::pair<Map_ptr_type, Map_ptr_type>
createMapsFromContiguousGIDs(const Mesh_type& mesh, const Entity_kind kind);


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
// Mesh deformation
// -----------------------------------------------------------------------------
template <class Mesh_type>
void
deform(Mesh_type& mesh,
       const typename Mesh_type::cEntity_ID_View& nodeids,
       const typename Mesh_type::cPoint_View& newpos);


} // namespace AmanziMesh
} // namespace Amanzi
