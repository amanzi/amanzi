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
#include "MeshCache_impl.hh"

namespace Amanzi {
namespace AmanziMesh {

  std::pair<double, AmanziGeometry::Point>
  MeshFrameworkAlgorithms::computeCellGeometry(const Mesh& mesh, const Entity_ID c) const
  {
    return MeshAlgorithms::computeCellGeometry(mesh, c);
  }

  std::tuple<double, AmanziGeometry::Point, cPoint_View>
  MeshFrameworkAlgorithms::computeFaceGeometry(const Mesh& mesh, const Entity_ID f) const
  {
    return MeshAlgorithms::computeFaceGeometry(mesh, f);
  }

  std::pair<AmanziGeometry::Point, AmanziGeometry::Point>
  MeshFrameworkAlgorithms::computeEdgeGeometry(const Mesh& mesh, const Entity_ID c) const
  {
    return MeshAlgorithms::computeEdgeGeometry(mesh, c);
  }


  double MeshFrameworkAlgorithms::getCellVolume(const Mesh& mesh, const Entity_ID c) const 
  { 
    return MeshAlgorithms::computeCellGeometry(mesh, c).first; 
  }
  
  AmanziGeometry::Point MeshFrameworkAlgorithms::getCellCentroid(const Mesh& mesh, const Entity_ID c) const
  {
    return MeshAlgorithms::computeCellGeometry(mesh, c).second;
  }

  double MeshFrameworkAlgorithms::getFaceArea(const Mesh& mesh, const Entity_ID f) const
  {
    return std::get<0>(MeshAlgorithms::computeFaceGeometry(mesh, f));
  }

  AmanziGeometry::Point MeshFrameworkAlgorithms::getFaceCentroid(const Mesh& mesh, const Entity_ID f) const
  {
    return std::get<1>(MeshAlgorithms::computeFaceGeometry(mesh, f));
  }

  AmanziGeometry::Point MeshFrameworkAlgorithms::getFaceNormal(
      const Mesh& mesh, const Entity_ID f, const Entity_ID c, int * const orientation) const
  {
    auto geom = MeshAlgorithms::computeFaceGeometry(mesh, f);

    cEntity_ID_View fcells;
    mesh.getFaceCells(f, Parallel_kind::ALL, fcells);
    if (orientation) *orientation = 0;
    Entity_ID cc = (c < 0) ? fcells[0] : c;

    int i = std::find(fcells.begin(), fcells.end(), cc) - fcells.begin();
    AmanziGeometry::Point normal = std::get<2>(geom)[i];

    if (mesh.getSpaceDimension() == mesh.getManifoldDimension()) {
      if (c < 0) {
        normal *= MeshAlgorithms::getFaceDirectionInCell(mesh, f, cc);
      } else if (orientation) {
        *orientation = MeshAlgorithms::getFaceDirectionInCell(mesh, f, cc);
      }
    } else if (c < 0) {
      Errors::Message msg("MeshFramework: asking for the natural normal of a submanifold mesh is not valid.");
      Exceptions::amanzi_throw(msg);
    }
    return normal;
  }

  double
  MeshFrameworkAlgorithms::getEdgeLength(const Mesh& mesh, const Entity_ID e) const
  {
    return AmanziGeometry::norm(MeshAlgorithms::computeEdgeGeometry(mesh, e).first);
  }

  AmanziGeometry::Point 
  MeshFrameworkAlgorithms::getEdgeVector(const Mesh& mesh, const Entity_ID e, const Entity_ID n, int * const orientation) const
  {
    auto geom = MeshAlgorithms::computeEdgeGeometry(mesh, e);
    if (n >= 0) {
      cEntity_ID_View nodes;
      mesh.getEdgeNodes(e, nodes);
      if (n == nodes[0]) {
        if (orientation) *orientation = 1;
        return geom.first;
      } else if (n == nodes[1]) {
        if (orientation) *orientation = -1;
        return -geom.first;
      } else {
        AMANZI_ASSERT(0);
      }
    }
    return geom.first;
  }

  AmanziGeometry::Point
  MeshFrameworkAlgorithms::getEdgeCentroid(const Mesh& mesh, const Entity_ID e) const
  {
    return MeshAlgorithms::computeEdgeGeometry(mesh, e).second;
  }

  void MeshFrameworkAlgorithms::getCellFacesAndBisectors(const Mesh& mesh, const Entity_ID cellid,
          cEntity_ID_View& faceids, cPoint_View * const bisectors) const
  {
    mesh.getCellFaces(cellid, faceids);
    if (bisectors)
      *bisectors = MeshAlgorithms::computeBisectors(mesh, cellid, faceids);
  }


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
  cEntity_ID_View cells;
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
