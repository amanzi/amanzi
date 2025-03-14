/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Rao Garimella
           Ethan Coon (coonet@ornl.gov)
*/

/*
  Mesh operations that are commonly needed by users of the library, both
  collective and "local" calculations.  These only use the public API of the
  Mesh, and so are implemented as non-member functions.

  While they can be called by developers of the Mesh library, typically they
  are used by clients of the library -- they are intended for external use.

  Note, these are templated on Mesh_type because they may be used by either
  MeshFramework objects or MeshCache objects.

*/

#pragma once

#include "Point.hh"
#include "MeshDefs.hh"

class Epetra_MultiVector;

namespace Amanzi {
namespace AmanziMesh {

// -----------------------------------------------------------------------------
// Derived topological relationships
//
// Non-collective functions
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
getFaceNormalExterior(const Mesh_type& mesh, Entity_ID f, int* dir);


//
// Get the directional int for a face that is on the boundary.
//
template <class Mesh_type>
KOKKOS_INLINE_FUNCTION int
getBoundaryDirection(const Mesh_type& mesh, Entity_ID f);


//
// Given a cell and a face bordering that cell, return the cell on the other
// side of that face.
//
template <class Mesh_type>
KOKKOS_INLINE_FUNCTION Entity_ID
getFaceAdjacentCell(const Mesh_type& mesh, int c, int f);


//
// Given a cell, return all neighboring cells.
//
template <class Mesh_type>
typename Mesh_type::Entity_ID_View
getCellFaceAdjacentCells(const Mesh_type& mesh, Entity_ID c, Parallel_kind ptype);


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
// Mesh deformation
// -----------------------------------------------------------------------------
template <class Mesh_type>
void
deform(Mesh_type& mesh,
       const typename Mesh_type::cEntity_ID_View& nodeids,
       const typename Mesh_type::cPoint_View& newpos);


} // namespace AmanziMesh
} // namespace Amanzi
