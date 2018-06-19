/*
  WhetStone, version 2.1
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

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Return centroid weights for 2D and 3D polygons. We take list of 
* nodes as input since it is not cached by the mesh.
* NOTE: polygon must be star-shaped w.r.t. its geometric center.
****************************************************************** */
inline
void PolygonCentroidWeights(
    const AmanziMesh::Mesh& mesh, const AmanziMesh::Entity_ID_List& nodes,
    double area, std::vector<double>& weights)
{
  int d = mesh.space_dimension();
  int nnodes = nodes.size();

  weights.assign(nnodes, 1.0 / (3 * nnodes));
  AmanziGeometry::Point p1(d), p2(d), p3(d), xg(d);

  // geometric center
  for (int i = 0; i < nnodes; ++i) {
    mesh.node_get_coordinates(nodes[i], &p1);
    xg += p1;
  }
  xg /= nnodes;

  // corner volume contributions
  for (int i1 = 0; i1 < nnodes; ++i1) {
    int i2 = (i1 + 1) % nnodes;  
    mesh.node_get_coordinates(nodes[i1], &p1);
    mesh.node_get_coordinates(nodes[i2], &p2);

    p3 = (p1 - xg)^(p2 - xg);
    double tmp = norm(p3) / (6 * area);

    weights[i1] += tmp;
    weights[i2] += tmp;
  }
}


/* ******************************************************************
* Geometric center of a mesh cell 
****************************************************************** */
inline
AmanziGeometry::Point cell_geometric_center(const AmanziMesh::Mesh& mesh, int c)
{
  int d = mesh.space_dimension();
  AmanziGeometry::Point v(d), xg(d);

  AmanziMesh::Entity_ID_List nodes;
  mesh.cell_get_nodes(c, &nodes);
  int nnodes = nodes.size();

  for (int i = 0; i < nnodes; ++i) {
    mesh.node_get_coordinates(nodes[i], &v);
    xg += v;
  } 
  xg /= nnodes;
  
  return xg;
}

}  // namespace WhetStone
}  // namespace Amanzi

#endif

