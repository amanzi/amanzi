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

#include "Mesh_Algorithms.hh"

namespace Amanzi {
namespace AmanziMesh {


// -----------------------------------------------------------------------------
// Given a boundary face ID, get the corresponding face ID
// -----------------------------------------------------------------------------
AmanziMesh::Entity_ID
getBoundaryFaceFace(const AmanziMesh::Mesh& mesh, AmanziMesh::Entity_ID bf)
{
  const auto& fmap = mesh.face_map(true);
  const auto& bfmap = mesh.exterior_face_map(true);
  return fmap.LID(bfmap.GID(bf));
}

// -----------------------------------------------------------------------------
// Given a face ID, get the corresponding boundary face ID (assuming it is a bf)
// -----------------------------------------------------------------------------
AmanziMesh::Entity_ID
getFaceOnBoundaryBoundaryFace(const AmanziMesh::Mesh& mesh, AmanziMesh::Entity_ID f)
{
  const auto& fmap = mesh.face_map(true);
  const auto& bfmap = mesh.exterior_face_map(true);
  return bfmap.LID(fmap.GID(f));
}

// -----------------------------------------------------------------------------
// Given a boundary face ID, get the cell internal to that face.
// -----------------------------------------------------------------------------
AmanziMesh::Entity_ID
getBoundaryFaceInternalCell(const AmanziMesh::Mesh& mesh, AmanziMesh::Entity_ID bf)
{
  return getFaceOnBoundaryInternalCell(mesh, getBoundaryFaceFace(mesh, bf));
}


// -----------------------------------------------------------------------------
// Given a face ID, and assuming it is a boundary face, get the cell internal.
// -----------------------------------------------------------------------------
AmanziMesh::Entity_ID
getFaceOnBoundaryInternalCell(const AmanziMesh::Mesh& mesh, AmanziMesh::Entity_ID f)
{
  AmanziMesh::Entity_ID_List cells;
  mesh.face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
  if (cells.size() != 1) {
    AmanziGeometry::Point fc = mesh.face_centroid(f);
    std::stringstream msgs;
    msgs << "getFaceOnBoundaryInternalCell called with non-internal face GID "
        << mesh.face_map(true).GID(f) << " at " << fc;
    Errors::Message msg(msgs.str());
    Exceptions::amanzi_throw(msg);
  }
  return cells[0];
}

// -----------------------------------------------------------------------------
// Given a vector on faces, import to vector on boundary faces
// -----------------------------------------------------------------------------
void
copyFacesToBoundaryFaces(const AmanziMesh::Mesh& mesh,
                         const Epetra_MultiVector& faces,
                         Epetra_MultiVector& boundary_faces)
{
  int ierr = boundary_faces.Import(faces, mesh.exterior_face_importer(), Insert);
  AMANZI_ASSERT(!ierr);
}

// -----------------------------------------------------------------------------
// Given a vector on faces, import to vector on boundary faces
// -----------------------------------------------------------------------------
void
copyBoundaryFacesToFaces(const AmanziMesh::Mesh& mesh,
                         const Epetra_MultiVector& boundary_faces,
                         Epetra_MultiVector& faces)
{
  int ierr = faces.Export(boundary_faces, mesh.exterior_face_importer(), Insert);
  AMANZI_ASSERT(!ierr);
}

// -----------------------------------------------------------------------------
// Given a vector on cells, set the boundary_face entries by their internal cell
// -----------------------------------------------------------------------------
void
copyCellsToBoundaryFaces(const AmanziMesh::Mesh& mesh,
                         const Epetra_MultiVector& cells,
                         Epetra_MultiVector& boundary_faces)
{
  AMANZI_ASSERT(cells.NumVectors() == boundary_faces.NumVectors());
  for (AmanziMesh::Entity_ID bf=0; bf!=boundary_faces.MyLength(); ++bf) {
    AmanziMesh::Entity_ID c = getBoundaryFaceInternalCell(mesh, bf);
    for (int i=0; i!=boundary_faces.NumVectors(); ++i) {
      boundary_faces[i][bf] = cells[i][c];
    }
  }
}


// -----------------------------------------------------------------------------
// Exterior boundary normal: dir = 0 for internal face
// -----------------------------------------------------------------------------
AmanziGeometry::Point
getFaceNormalExterior(const AmanziMesh::Mesh& mesh, int f, int* dir)
{
  Entity_ID_List cells;
  mesh.face_get_cells(f, Parallel_type::ALL, &cells);

  auto normal = mesh.face_normal(f, false, cells[0], dir);
  if (cells.size() > 1) *dir = 0;

  return normal;
}


// -----------------------------------------------------------------------------
// Given a cell c and face f, returns the neighbooring cell
// -----------------------------------------------------------------------------
int
cell_get_face_adj_cell(const AmanziMesh::Mesh& mesh, int c, int f)
{
  Entity_ID_List cells;
  mesh.face_get_cells(f, Parallel_type::ALL, &cells);

  if (cells.size() == 2)
    return cells[0] + cells[1] - c;

  return -1;
}


} // namespace AmanziMesh
} // namespace Amanzi
