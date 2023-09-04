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
//#include "Mesh.hh"
#include "MeshCache_decl.hh"
#include "Mesh_Helpers.hh"

namespace Amanzi {
namespace AmanziMesh {

using Mesh = MeshCache<MemSpace_kind::HOST>;

//
// This class provides default, virtual algorithms for computing geometric
// quantities given nodal coordinates and topological information.
//
// Split into two classes to aid in deletion of the MeshFramework class, while
// keeping the MeshFrameworkAlgorithms class around for use by the Cache.
//
struct MeshFrameworkAlgorithms {
  virtual ~MeshFrameworkAlgorithms(){};

  // lumped things for more efficient calculation
  virtual std::pair<double, AmanziGeometry::Point>
  computeCellGeometry(const Mesh& mesh, const Entity_ID c) const;

  virtual std::tuple<double, AmanziGeometry::Point, cPoint_View>
  computeFaceGeometry(const Mesh& mesh, const Entity_ID f) const;

  virtual std::pair<AmanziGeometry::Point, AmanziGeometry::Point>
  computeEdgeGeometry(const Mesh& mesh, const Entity_ID e) const;

  virtual double getCellVolume(const Mesh& mesh, const Entity_ID c) const;
  virtual AmanziGeometry::Point getCellCentroid(const Mesh& mesh, const Entity_ID c) const;

  virtual double getFaceArea(const Mesh& mesh, const Entity_ID f) const;
  virtual AmanziGeometry::Point getFaceCentroid(const Mesh& mesh, const Entity_ID f) const;
  virtual AmanziGeometry::Point getFaceNormal(const Mesh& mesh,
                                              const Entity_ID f,
                                              const Entity_ID c,
                                              int* const orientation) const;
  virtual double getEdgeLength(const Mesh& mesh, const Entity_ID e) const;
  virtual AmanziGeometry::Point getEdgeVector(const Mesh& mesh,
                                              const Entity_ID e,
                                              const Entity_ID n,
                                              int* const orientation) const;

  virtual AmanziGeometry::Point getEdgeCentroid(const Mesh& mesh, const Entity_ID e) const;

  virtual void getCellFacesAndBisectors(const Mesh& mesh,
                                        const Entity_ID cellid,
                                        cEntity_ID_View& faceids,
                                        cPoint_View* const bisectors) const;
};


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
