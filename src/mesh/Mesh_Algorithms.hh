/*
  Copyright 2010-201x held jointly by LANL, ORNL, LBNL, and PNNL.
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

#include "Epetra_MultiVector.h"
#include "Mesh.hh"

namespace Amanzi {
namespace AmanziMesh {

// -----------------------------------------------------------------------------
// Given a boundary face ID, get the corresponding face ID
// -----------------------------------------------------------------------------
Entity_ID
getBoundaryFaceFace(const Mesh& mesh, Entity_ID bf);

// -----------------------------------------------------------------------------
// Given a face ID, get the corresponding boundary face ID (assuming it is a bf)
// -----------------------------------------------------------------------------
Entity_ID
getFaceOnBoundaryBoundaryFace(const Mesh& mesh, Entity_ID f);

// -----------------------------------------------------------------------------
// Given a boundary face ID, get the cell internal to that face.
// -----------------------------------------------------------------------------
Entity_ID
getBoundaryFaceInternalCell(const Mesh& mesh, Entity_ID bf);

// -----------------------------------------------------------------------------
// Given a face ID, and assuming it is a boundary face, get the cell internal.
// -----------------------------------------------------------------------------
Entity_ID
getFaceOnBoundaryInternalCell(const Mesh& mesh, Entity_ID f);

// -----------------------------------------------------------------------------
// Given a vector on faces, import to vector on boundary faces
// -----------------------------------------------------------------------------
void
copyFacesToBoundaryFaces(const Mesh& mesh,
                         const Epetra_MultiVector& faces,
                         Epetra_MultiVector& boundary_faces);

// -----------------------------------------------------------------------------
// Given a vector on faces, import to vector on boundary faces
// -----------------------------------------------------------------------------
void
copyBoundaryFacesToFaces(const Mesh& mesh,
                         const Epetra_MultiVector& boundary_faces,
                         Epetra_MultiVector& faces);

// -----------------------------------------------------------------------------
// Given a vector on cells, set the boundary_face entries by their internal cell
// -----------------------------------------------------------------------------
void
copyCellsToBoundaryFaces(const Mesh& mesh,
                         const Epetra_MultiVector& cells,
                         Epetra_MultiVector& boundary_faces);

// -----------------------------------------------------------------------------
// Given a boundary face f, return the exterior normal. If f is an interior face,
// dir = 0 and normal orientation is not be reliable in parallel algorithms
// -----------------------------------------------------------------------------
AmanziGeometry::Point
getFaceNormalExterior(const Mesh& mesh, int f, int* dir);




} // namespace AmanziMesh
} // namespace Amanzi
