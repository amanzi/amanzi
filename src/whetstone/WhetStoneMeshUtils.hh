/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Additional mesh routines.
*/

#ifndef AMANZI_WHETSTONE_MESH_UTILS_HH_
#define AMANZI_WHETSTONE_MESH_UTILS_HH_

#include "Teuchos_RCP.hpp"

// Amanzi
#include "MeshDefs.hh"
#include "Mesh.hh"
#include "Point.hh"

#include "WhetStoneDefs.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Return centroid weights for 2D and 3D polygons. We take list of
* nodes as input since it is not cached by the mesh.
* NOTE: polygon must be star-shaped w.r.t. its geometric center.
****************************************************************** */
inline void
PolygonCentroidWeights(const AmanziMesh::Mesh& mesh,
                       const AmanziMesh::cEntity_ID_View& nodes,
                       double area,
                       AmanziMesh::Double_List& weights)
{
  int d = mesh.getSpaceDimension();
  int nnodes = nodes.size();

  weights.assign(nnodes, 1.0 / (3 * nnodes));
  AmanziGeometry::Point p1(d), p2(d), p3(d), xg(d);

  // geometric center
  for (int i = 0; i < nnodes; ++i) {
    p1 = mesh.getNodeCoordinate(nodes[i]);
    xg += p1;
  }
  xg /= nnodes;

  // corner volume contributions
  for (int i1 = 0; i1 < nnodes; ++i1) {
    int i2 = (i1 + 1) % nnodes;
    p1 = mesh.getNodeCoordinate(nodes[i1]);
    p2 = mesh.getNodeCoordinate(nodes[i2]);

    p3 = (p1 - xg) ^ (p2 - xg);
    double tmp = norm(p3) / (6 * area);

    weights[i1] += tmp;
    weights[i2] += tmp;
  }
}


/* ******************************************************************
// Faces of ptype of cell c that are connected to node v.
****************************************************************** */
inline void
node_get_cell_faces(const AmanziMesh::Mesh& mesh,
                    const AmanziMesh::Entity_ID v,
                    const AmanziMesh::Entity_ID c,
                    const AmanziMesh::Parallel_kind ptype,
                    AmanziMesh::Entity_ID_View* faces)
{
  AmanziMesh::Entity_ID_List vfaces;
  int nfaces_owned =
    mesh.getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED);

  const auto& faces_tmp = mesh.getCellFaces(c);
  int nfaces = faces_tmp.size();

  for (int i = 0; i < nfaces; ++i) {
    int f = faces_tmp[i];
    if (ptype == AmanziMesh::Parallel_kind::OWNED && f >= nfaces_owned) continue;

    auto nodes = mesh.getFaceNodes(f);
    int nnodes = nodes.size();
    for (int k = 0; k < nnodes; ++k) {
      if (nodes[k] == v) {
        vfaces.push_back(f);
        break;
      }
    }
  }
  vectorToView(*faces, vfaces);
}


/* ******************************************************************
* Extension of Mesh API: entities in cell
****************************************************************** */
inline void
cell_get_entities(const AmanziMesh::Mesh& mesh,
                  int c,
                  const AmanziMesh::Entity_kind kind,
                  AmanziMesh::cEntity_ID_View* entities)
{
  if (kind == AmanziMesh::Entity_kind::FACE) {
    *entities = mesh.getCellFaces(c);
  } else if (kind == AmanziMesh::Entity_kind::EDGE) {
    *entities = mesh.getCellEdges(c);
  } else if (kind == AmanziMesh::Entity_kind::NODE) {
    *entities = mesh.getCellNodes(c);
  } else if (kind == AmanziMesh::Entity_kind::CELL) {
    AmanziMesh::Entity_ID_View lentities;
    Kokkos::resize(lentities, 1);
    lentities[0] = c;
    *entities = lentities;
  } else {
    AmanziMesh::Entity_ID_View lentities;
    *entities = lentities;
  }
}


/* ******************************************************************
* Extension of Mesh API: entities in cell
****************************************************************** */
inline AmanziGeometry::Point
getNodeUnitNormal(const AmanziMesh::Mesh& mesh, int v)
{
  Entity_ID_List faces, cells;
  mesh.node_get_faces(v, AmanziMesh::Parallel_type::ALL, &faces);
  int nfaces = faces.size();

  int d = mesh.space_dimension();
  AmanziGeometry::Point normal(d);

  for (int n = 0; n < nfaces; ++n) {
    int f = faces[n];
    mesh.face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
    if (cells.size() == 1) normal += mesh.face_normal(faces[n]);
  }
  return normal / norm(normal);
}

} // namespace WhetStone
} // namespace Amanzi

#endif
