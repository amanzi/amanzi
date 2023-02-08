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
// Copyright 2010-201x held jointly by LANL, ORNL, LBNL, and PNNL.
// Amanzi is released under the three-clause BSD License.
// The terms of use and "as is" disclaimer for this license are
// provided in the top-level COPYRIGHT file.
//
// Algorithms for commonly used Mesh operations -- for the Mesh user.
//

#include "Mesh_Algorithms.hh"

namespace Amanzi {
namespace AmanziMesh {

// -----------------------------------------------------------------------------
// Given a boundary face ID, get the corresponding face ID
// -----------------------------------------------------------------------------
Entity_ID
getBoundaryFaceFace(const Mesh& mesh, Entity_ID bf)
{
  const auto& fmap = mesh.getMap(AmanziMesh::Entity_kind::FACE, true);
  const auto& bfmap = mesh.getMap(AmanziMesh::Entity_kind::BOUNDARY_FACE, true);
  return fmap.LID(bfmap.GID(bf));
}

// -----------------------------------------------------------------------------
// Given a face ID, get the corresponding boundary face ID (assuming it is a bf)
// -----------------------------------------------------------------------------
Entity_ID
getFaceOnBoundaryBoundaryFace(const Mesh& mesh, Entity_ID f)
{
  const auto& fmap = mesh.getMap(AmanziMesh::Entity_kind::FACE, true);
  const auto& bfmap = mesh.getMap(AmanziMesh::Entity_kind::BOUNDARY_FACE, true);
  return bfmap.LID(fmap.GID(f));
}

// -----------------------------------------------------------------------------
// Given a boundary face ID, get the cell internal to that face.
// -----------------------------------------------------------------------------
Entity_ID
getBoundaryFaceInternalCell(const Mesh& mesh, Entity_ID bf)
{
  return getFaceOnBoundaryInternalCell(mesh, getBoundaryFaceFace(mesh, bf));
}


// -----------------------------------------------------------------------------
// Given a face ID, and assuming it is a boundary face, get the cell internal.
// -----------------------------------------------------------------------------
Entity_ID
getFaceOnBoundaryInternalCell(const Mesh& mesh, Entity_ID f)
{
  Entity_ID_View cells;
  mesh.getFaceCells(f, Parallel_kind::ALL, cells);
  if (cells.size() != 1) {
    AmanziGeometry::Point fc = mesh.getFaceCentroid(f);
    std::stringstream msgs;
    msgs << "getFaceOnBoundaryInternalCell called with non-internal face GID "
        << mesh.getMap(AmanziMesh::Entity_kind::FACE, true).GID(f) << " at " << fc;
    Errors::Message msg(msgs.str());
    Exceptions::amanzi_throw(msg);
  }
  return cells[0];
}

// -----------------------------------------------------------------------------
// Given a vector on faces, import to vector on boundary faces
// -----------------------------------------------------------------------------
void
copyFacesToBoundaryFaces(const Mesh& mesh,
                         const Epetra_MultiVector& faces,
                         Epetra_MultiVector& boundary_faces)
{
  int ierr = boundary_faces.Import(faces, mesh.getBoundaryFaceImporter(), Insert);
  AMANZI_ASSERT(!ierr);
}

// -----------------------------------------------------------------------------
// Given a vector on faces, import to vector on boundary faces
// -----------------------------------------------------------------------------
void
copyBoundaryFacesToFaces(const Mesh& mesh,
                         const Epetra_MultiVector& boundary_faces,
                         Epetra_MultiVector& faces)
{
  int ierr = faces.Export(boundary_faces, mesh.getBoundaryFaceImporter(), Insert);
  AMANZI_ASSERT(!ierr);
}

// -----------------------------------------------------------------------------
// Given a vector on cells, set the boundary_face entries by their internal cell
// -----------------------------------------------------------------------------
void
copyCellsToBoundaryFaces(const Mesh& mesh,
                         const Epetra_MultiVector& cells,
                         Epetra_MultiVector& boundary_faces)
{
  AMANZI_ASSERT(cells.NumVectors() == boundary_faces.NumVectors());
  for (Entity_ID bf=0; bf!=boundary_faces.MyLength(); ++bf) {
    Entity_ID c = getBoundaryFaceInternalCell(mesh, bf);
    for (int i=0; i!=boundary_faces.NumVectors(); ++i) {
      boundary_faces[i][bf] = cells[i][c];
    }
  }
}

// -----------------------------------------------------------------------------
// Exterior boundary normal: dir = 0 for internal face
// -----------------------------------------------------------------------------
AmanziGeometry::Point
getFaceNormalExterior(const Mesh& mesh, int f, int* dir)
{
auto cells = mesh.getFaceCells(f, Parallel_kind::ALL);

  auto normal = mesh.getFaceNormal(f, cells[0], dir);
  if (cells.size() > 1) *dir = 0;

  return normal;
}



} // namespace AmanziMesh
} // namespace Amanzi
