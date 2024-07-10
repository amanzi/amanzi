/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/
/*

MeshAlgorithms are algorithms that are specific to the MeshFramework
implementation.  They extend the MeshFramework API and provide the COMPUTE
algorithm for the MeshCache implementaiton.  They are kept separate from
the MeshFramework API because then the MeshFramework API can be deleted and
these can be used with the MeshCache instead.

Note that the "virtual" algorithms are kept in a separate class to enable
overriding them for special meshes (MeshLogical).

*/

#include "MeshAlgorithms.hh"
#include "MeshInternals.hh"
#include "Mesh.hh"

namespace Amanzi {
namespace AmanziMesh {


AmanziGeometry::Point
MeshAlgorithms::computeFaceNormal(const Mesh& mesh,
                                  const Entity_ID f,
                                  const Entity_ID c,
                                  int* const orientation) const
{
  auto geom = Impl::computeFaceGeometry(mesh, f);

  Mesh::cEntity_ID_View fcells = mesh.getFaceCells(f);
  if (orientation) *orientation = 0;

  Entity_ID cc;
  std::size_t i;
  if (c < 0) {
    cc = fcells[0];
    i = 0;
  } else {
    cc = c;
    auto ncells = fcells.size();
    for (i = 0; i != ncells; ++i)
      if (fcells[i] == cc) break;
  }
  AmanziGeometry::Point normal = std::get<2>(geom)[i];

  if (mesh.getSpaceDimension() == mesh.getManifoldDimension()) {
    if (c < 0) {
      assert(orientation == nullptr);
      normal *= Impl::getFaceDirectionInCell(mesh, f, cc);
    } else if (orientation) {
      *orientation = Impl::getFaceDirectionInCell(mesh, f, cc);
    }
  } else {
    // manifold case
    if (c < 0) {
      assert(orientation == nullptr);

      if (fcells.size() != 2) {
        normal *= Impl::getFaceDirectionInCell(mesh, f, cc);
      } else {
        // average normals oriented from lower to higher GIDs
        int pos_i = mesh.getEntityGID(Entity_kind::CELL, fcells[0]) >
                        mesh.getEntityGID(Entity_kind::CELL, fcells[1]) ?
                      0 :
                      1;
        normal = (std::get<2>(geom)[1 - pos_i] - std::get<2>(geom)[pos_i]) / 2;
      }
    } else if (orientation) {
      *orientation = Impl::getFaceDirectionInCell(mesh, f, cc);
    }
  }

  if (orientation) assert(*orientation != 0);
  return normal;
}


AmanziGeometry::Point
MeshAlgorithms::computeEdgeVector(const Mesh& mesh,
                                  const Entity_ID e,
                                  const Entity_ID n,
                                  int* const orientation) const
{
  auto geom = Impl::computeEdgeGeometry(mesh, e);
  if (n >= 0) {
    auto nodes = mesh.getEdgeNodes(e);
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


void
MeshAlgorithms::computeCellFacesAndBisectors(const Mesh& mesh,
                                             const Entity_ID cellid,
                                             typename Mesh::cEntity_ID_View& faceids,
                                             typename Mesh::cPoint_View* const bisectors) const
{
  mesh.getCellFaces(cellid, faceids);
  if (bisectors) *bisectors = Impl::computeBisectors(mesh, cellid, faceids);
}


} // namespace AmanziMesh
} // namespace Amanzi
