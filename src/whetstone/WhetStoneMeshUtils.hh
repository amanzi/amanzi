/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Additional mesh routines.
*/

#ifndef AMANZI_WHETSTONE_MESH_UTILS_HH_
#define AMANZI_WHETSTONE_MESH_UTILS_HH_

#include "Teuchos_RCP.hpp"

// Amanzi
#include "Mesh.hh"
#include "MeshDefs.hh"
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
                       const Kokkos::View<AmanziMesh::Entity_ID*>& nodes,
                       double area, std::vector<double>& weights)
{
  int d = mesh.space_dimension();
  int nnodes = nodes.extent(0);

  weights.assign(nnodes, 1.0 / (3 * nnodes));
  AmanziGeometry::Point p1(d), p2(d), p3(d), xg(d);

  // geometric center
  for (int i = 0; i < nnodes; ++i) {
    mesh.node_get_coordinates(nodes(i), &p1);
    xg += p1;
  }
  xg /= nnodes;

  // corner volume contributions
  for (int i1 = 0; i1 < nnodes; ++i1) {
    int i2 = (i1 + 1) % nnodes;
    mesh.node_get_coordinates(nodes(i1), &p1);
    mesh.node_get_coordinates(nodes(i2), &p2);

    p3 = (p1 - xg) ^ (p2 - xg);
    double tmp = norm(p3) / (6 * area);

    weights[i1] += tmp;
    weights[i2] += tmp;
  }
}


/* ******************************************************************
 * Extension of Mesh API.
 ****************************************************************** */
inline int
cell_get_face_adj_cell(const AmanziMesh::Mesh& mesh, int c, int f)
{
  Kokkos::View<Amanzi::AmanziMesh::Entity_ID*> cells;
  mesh.face_get_cells(f, Parallel_type::ALL, cells);

  if (cells.extent(0) == 2) return cells(0) + cells(1) - c;

  return -1;
}


/* ******************************************************************
 * Exterior boundary normal: dir = 0 for internal face
 ****************************************************************** */
inline AmanziGeometry::Point
face_normal_exterior(const AmanziMesh::Mesh& mesh, int f, int* dir)
{
  Kokkos::View<Amanzi::AmanziMesh::Entity_ID*> cells;
  mesh.face_get_cells(f, Amanzi::AmanziMesh::Parallel_type::ALL, cells);

  auto normal = mesh.face_normal(f, false, cells(0), dir);
  if (cells.extent(0) > 1) *dir = 0;

  return normal;
}


/* ******************************************************************
 * Geometric center of a mesh cell
 ****************************************************************** */
inline AmanziGeometry::Point
cell_geometric_center(const AmanziMesh::Mesh& mesh, int c)
{
  int d = mesh.space_dimension();
  AmanziGeometry::Point v(d), xg(d);

  Kokkos::View<AmanziMesh::Entity_ID*> nodes;
  mesh.cell_get_nodes(c, nodes);
  int nnodes = nodes.extent(0);

  for (int i = 0; i < nnodes; ++i) {
    mesh.node_get_coordinates(nodes(i), &v);
    xg += v;
  }
  xg /= nnodes;

  return xg;
}

} // namespace WhetStone
} // namespace Amanzi

#endif
