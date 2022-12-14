/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

//
// Mesh
//
// Helper functions for Mesh operations and algorithms
//
#pragma once

#include "Epetra_MultiVector.h"
#include "Mesh.hh"

namespace Amanzi {
namespace AmanziMesh {

// -----------------------------------------------------------------------------
// Given a boundary face ID, get the corresponding face ID
// -----------------------------------------------------------------------------
AmanziMesh::Entity_ID
getBoundaryFaceFace(const AmanziMesh::Mesh& mesh, AmanziMesh::Entity_ID bf);

// -----------------------------------------------------------------------------
// Given a face ID, get the corresponding boundary face ID (assuming it is a bf)
// -----------------------------------------------------------------------------
AmanziMesh::Entity_ID
getFaceOnBoundaryBoundaryFace(const AmanziMesh::Mesh& mesh, AmanziMesh::Entity_ID f);

// -----------------------------------------------------------------------------
// Given a boundary face ID, get the cell internal to that face.
// -----------------------------------------------------------------------------
AmanziMesh::Entity_ID
getBoundaryFaceInternalCell(const AmanziMesh::Mesh& mesh, AmanziMesh::Entity_ID bf);

// -----------------------------------------------------------------------------
// Given a face ID, and assuming it is a boundary face, get the cell internal.
// -----------------------------------------------------------------------------
AmanziMesh::Entity_ID
getFaceOnBoundaryInternalCell(const AmanziMesh::Mesh& mesh, AmanziMesh::Entity_ID f);

// -----------------------------------------------------------------------------
// Given a vector on faces, import to vector on boundary faces
// -----------------------------------------------------------------------------
void
copyFacesToBoundaryFaces(const AmanziMesh::Mesh& mesh,
                         const Epetra_MultiVector& faces,
                         Epetra_MultiVector& boundary_faces);

// -----------------------------------------------------------------------------
// Given a vector on faces, import to vector on boundary faces
// -----------------------------------------------------------------------------
void
copyBoundaryFacesToFaces(const AmanziMesh::Mesh& mesh,
                         const Epetra_MultiVector& boundary_faces,
                         Epetra_MultiVector& faces);

// -----------------------------------------------------------------------------
// Given a vector on cells, set the boundary_face entries by their internal cell
// -----------------------------------------------------------------------------
void
copyCellsToBoundaryFaces(const AmanziMesh::Mesh& mesh,
                         const Epetra_MultiVector& cells,
                         Epetra_MultiVector& boundary_faces);


// -----------------------------------------------------------------------------
// Given a boundary face f, return the exterior normal. If f is an interior face,
// dir = 0 and normal orientation is not be reliable in parallel algorithms
// -----------------------------------------------------------------------------
AmanziGeometry::Point
getFaceNormalExterior(const AmanziMesh::Mesh& mesh, int f, int* dir);


// -----------------------------------------------------------------------------
// Given a cell c and face f, returns the neighbooring cell
// -----------------------------------------------------------------------------
int
cell_get_face_adj_cell(const AmanziMesh::Mesh& mesh, int c, int f);

} // namespace AmanziMesh
} // namespace Amanzi
