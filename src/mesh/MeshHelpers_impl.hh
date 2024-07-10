/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Rao Garimella
           Ethan Coon (coonet@ornl.gov)
*/

#pragma once

#include "Kokkos_StdAlgorithms.hpp"
#include "Tpetra_MultiVector.hpp"

#include "Point.hh"
#include "Geometry.hh"
#include "MeshDefs.hh"
#include "ViewUtils.hh"

#include "MeshHelpers_decl.hh"
#include "Mesh_decl.hh"
#include "MeshCache_decl.hh"

namespace Amanzi {
namespace AmanziMesh {

// -----------------------------------------------------------------------------
// Derived topological relationships
//
// Non-collective functions
// -----------------------------------------------------------------------------

//
// Given a face ID, get the corresponding boundary face ID (assuming it is a bf)
//
inline Entity_ID
getFaceOnBoundaryBoundaryFace(const Mesh& mesh, Entity_ID f)
{
  // should this be deprecated?  It seems likely to not perform well
  const auto& fmap = mesh.getMap(AmanziMesh::Entity_kind::FACE, true);
  const auto& bfmap = mesh.getMap(AmanziMesh::Entity_kind::BOUNDARY_FACE, true);
  return bfmap->getLocalElement(fmap->getGlobalElement(f));
}

KOKKOS_INLINE_FUNCTION
Entity_ID
getFaceOnBoundaryBoundaryFace(const MeshCache& mesh, Entity_ID f)
{
  // This implementation uses Kokkos::find, but that is host-only.  Instead,
  // there is a Kokkos::find based on Teams that we could probably use, but for
  // now we use a very simple hand-rolled for-loop based find.
  // auto iter = Kokkos::Experimental::find(DefaultExecutionSpace(), mesh.data.boundary_faces.view_device(), f);
  // assert(iter != Kokkos::Experimental::end(mesh.data.boundary_faces.view_device()));
  // return static_cast<Entity_ID>(iter - Kokkos::Experimental::begin(mesh.data.boundary_faces.view_device()));
  return find(mesh.data.boundary_faces.view_device(), f);
}


//
// Given a boundary face ID, get the cell internal to that face.
//
inline Entity_ID
getBoundaryFaceInternalCell(const Mesh& mesh, Entity_ID bf)
{
  return getFaceOnBoundaryInternalCell(mesh, mesh.getBoundaryFaceFace(bf));
}

KOKKOS_INLINE_FUNCTION
Entity_ID
getBoundaryFaceInternalCell(const MeshCache& mesh, Entity_ID bf)
{
  return getFaceOnBoundaryInternalCell(mesh, mesh.getBoundaryFaceFace(bf));
}

//
// Given a face ID, and assuming it is a boundary face, get the cell internal.
//
inline Entity_ID
getFaceOnBoundaryInternalCell(const Mesh& mesh, Entity_ID f) {
  assert(mesh.getFaceNumCells(f) == 1);
  return mesh.getFaceCell(f, 0);
}

KOKKOS_INLINE_FUNCTION
Entity_ID
getFaceOnBoundaryInternalCell(const MeshCache& mesh, const Entity_ID f)
{
  assert(mesh.getFaceNumCells(f) == 1);
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
  static_assert(Mesh_type::MEM == MemSpace_kind::HOST);

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
template <class Mesh_type>
void
copyFacesToBoundaryFaces(const Mesh_type& mesh,
                         const MultiVector_type& faces,
                         MultiVector_type& boundary_faces)
{
  boundary_faces.doImport(faces, mesh.getBoundaryFaceFaceImporter(), Tpetra::INSERT);
}


//
// Given a vector on faces, import to vector on boundary faces
//
template <class Mesh_type>
void
copyBoundaryFacesToFaces(const Mesh_type& mesh,
                         const MultiVector_type& boundary_faces,
                         MultiVector_type& faces)
{
  faces.doExport(boundary_faces, mesh.getBoundaryFaceImporter(), Tpetra::INSERT);
}


//
// Given a vector on cells, set the boundary_face entries by their internal cell
//
template <class Mesh_type>
void
copyCellsToBoundaryFaces(const Mesh_type& mesh,
                         const MultiVector_type& cells,
                         MultiVector_type& boundary_faces)
{
  boundary_faces.doImport(cells, mesh.getBoundaryFaceInternalCellImporter(), Tpetra::INSERT);
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
