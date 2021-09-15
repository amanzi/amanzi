//
// Mesh
//
// Copyright 2010-201x held jointly by LANL, ORNL, LBNL, and PNNL.
// Amanzi is released under the three-clause BSD License.
// The terms of use and "as is" disclaimer for this license are
// provided in the top-level COPYRIGHT file.
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



} // namespace AmanziMesh
} // namespace Amanzi
