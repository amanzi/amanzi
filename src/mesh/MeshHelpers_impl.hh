/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Rao Garimella
           Ethan Coon (coonet@ornl.gov)
*/

#pragma once

#include "Epetra_MultiVector.h"

#include "Point.hh"
#include "Geometry.hh"
#include "MeshDefs.hh"
#include "ViewUtils.hh"
#include "MeshUtils.hh"
#include "MeshCache_decl.hh"

#include "MeshHelpers_decl.hh"

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
getBoundaryFaceFace(const Mesh_type& mesh, Entity_ID bf)
{
  return mesh.getBoundaryFaces()(bf);
}


//
// Given a face ID, get the corresponding boundary face ID (assuming it is a bf)
//
template <MemSpace_kind MEM>
Entity_ID
getFaceOnBoundaryBoundaryFace(const MeshCache<MEM>& mesh, Entity_ID f)
{
  // should this be deprecated?  It seems likely to not perform well
  const auto& fmap = mesh.getMap(AmanziMesh::Entity_kind::FACE, true);
  const auto& bfmap = mesh.getMap(AmanziMesh::Entity_kind::BOUNDARY_FACE, true);
  return bfmap.LID(fmap.GID(f));
}


//
// Given a boundary face ID, get the cell internal to that face.
//
template <class Mesh_type>
KOKKOS_INLINE_FUNCTION Entity_ID
getBoundaryFaceInternalCell(const Mesh_type& mesh, Entity_ID bf)
{
  return getFaceOnBoundaryInternalCell(mesh, getBoundaryFaceFace(mesh, bf));
}


//
// Given a face ID, and assuming it is a boundary face, get the cell internal.
//
template <class Mesh_type>
KOKKOS_INLINE_FUNCTION Entity_ID
getFaceOnBoundaryInternalCell(const Mesh_type& mesh, Entity_ID f)
{
  if constexpr (std::is_same<Mesh_type, Mesh>::value) {
    assert(mesh.getFaceNumCells(f, Parallel_kind::ALL) == 1 &&
           "getFaceOnBoundaryInternalCell() requires a face on the boundary.");
  } else {
    if (mesh.getFaceNumCells(f, Parallel_kind::ALL) != 1) {
      AmanziGeometry::Point fc = mesh.getFaceCentroid(f);
      std::stringstream msgs;
      msgs << "getFaceOnBoundaryInternalCell called with non-internal face GID "
           << mesh.getMap(AmanziMesh::Entity_kind::FACE, true)->getGlobalElement(f) << " at " << fc;
      Errors::Message msg(msgs.str());
      Exceptions::amanzi_throw(msg);
    }
  }
  return mesh.getFaceCell(f, 0);
}


//
// Exterior boundary normal: dir = 0 for internal face
//
template <class Mesh_type>
KOKKOS_INLINE_FUNCTION AmanziGeometry::Point
getFaceNormalExterior(const Mesh_type& mesh, int f, int* dir)
{
  auto cells = mesh.getFaceCells(f);
  auto normal = mesh.getFaceNormal(f, cells[0], dir);
  if (cells.size() > 1) *dir = 0;
  return normal;
}


//
// Get the directional int for a face that is on the boundary.
//
template <class Mesh_type>
KOKKOS_INLINE_FUNCTION int
getBoundaryDirection(const Mesh_type& mesh, Entity_ID f)
{
  int orientation;
  getFaceNormalExterior(mesh, f, &orientation);
  return orientation;
}


//
// Given a cell and a face bordering that cell, return the cell on the other
// side of that face.
//
template <class Mesh_type>
int
getFaceAdjacentCell(const Mesh_type& mesh, int c, int f)
{
  auto cells = mesh.getFaceCells(f);
  if (cells.size() == 2) return cells[0] + cells[1] - c;
  return -1;
}


//
// Given a cell, return all neighboring cells.
//
template <class Mesh_type>
typename Mesh_type::Entity_ID_View
getCellFaceAdjacentCells(const Mesh_type& mesh, Entity_ID c, Parallel_kind ptype)
{
  static_assert(!std::is_same<Mesh_type, MeshCache<MemSpace_kind::DEVICE>>::value);
  auto cfaces = mesh.getCellFaces(c);
  typename Mesh_type::Entity_ID_View adj_cells;
  Entity_ID_List vadj_cells;
  for (const auto& f : cfaces) {
    auto fcells = mesh.getFaceCells(f);
    for (const auto& fc : fcells) {
      if (c != fc && std::find(vadj_cells.begin(), vadj_cells.end(), fc) == vadj_cells.end())
        vadj_cells.push_back(fc);
    }
  }
  vectorToView(adj_cells, vadj_cells);
  return adj_cells;
}


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
                         Epetra_MultiVector& boundary_faces)
{
  int ierr = boundary_faces.Import(faces, mesh.getBoundaryFaceImporter(), Insert);
  AMANZI_ASSERT(!ierr);
}


//
// Given a vector on faces, import to vector on boundary faces
//
template <MemSpace_kind MEM>
void
copyBoundaryFacesToFaces(const MeshCache<MEM>& mesh,
                         const Epetra_MultiVector& boundary_faces,
                         Epetra_MultiVector& faces)
{
  int ierr = faces.Export(boundary_faces, mesh.getBoundaryFaceImporter(), Insert);
  AMANZI_ASSERT(!ierr);
}


//
// Given a vector on cells, set the boundary_face entries by their internal cell
//
template <MemSpace_kind MEM>
void
copyCellsToBoundaryFaces(const MeshCache<MEM>& mesh,
                         const Epetra_MultiVector& cells,
                         Epetra_MultiVector& boundary_faces)
{
  AMANZI_ASSERT(cells.NumVectors() == boundary_faces.NumVectors());
  for (Entity_ID bf = 0; bf != boundary_faces.MyLength(); ++bf) {
    Entity_ID c = getBoundaryFaceInternalCell(mesh, bf);
    for (int i = 0; i != boundary_faces.NumVectors(); ++i) { boundary_faces[i][bf] = cells[i][c]; }
  }
}


// -----------------------------------------------------------------------------
// Mesh deformation
// -----------------------------------------------------------------------------
template <class Mesh_type>
void
deform(Mesh_type& mesh,
       const typename Mesh_type::cEntity_ID_View& nodeids,
       const typename Mesh_type::cPoint_View& newpos)
{
  mesh.setNodeCoordinates(nodeids, newpos);
}

} // namespace AmanziMesh
} // namespace Amanzi
