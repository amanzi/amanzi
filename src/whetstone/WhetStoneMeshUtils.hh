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
inline
void PolygonCentroidWeights(
    const AmanziMesh::Mesh& mesh, 
    const AmanziMesh::Entity_ID_List& nodes,
    double area, std::vector<double>& weights)
{
  int d = mesh.getSpaceDimension();
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
// Faces of ptype of cell c that are connected to node v.
****************************************************************** */
inline
void node_get_cell_faces(const AmanziMesh::Mesh& mesh,
                         const AmanziMesh::Entity_ID v, 
                         const AmanziMesh::Entity_ID c,
                         const AmanziMesh::Parallel_type ptype,
                         AmanziMesh::Entity_ID_List *faces) 
{
  Entity_ID_List nodes;

  faces->clear();  
  int nfaces_owned = mesh.getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_type::OWNED);

  const auto& faces_tmp = mesh.getCellFaces(c);
  int nfaces = faces_tmp.size();

  for (int i = 0; i < nfaces; ++i) {
    int f = faces_tmp[i];
    if (ptype == AmanziMesh::Parallel_type::OWNED && f >= nfaces_owned) continue;

    mesh.getFaceNodes(f, nodes);
    int nnodes = nodes.size();
    for (int k = 0; k < nnodes; ++k) {
      if (nodes[k] == v) {
        faces->push_back(f);
        break;
      }
    }
  }
}


/* ******************************************************************
* Extension of Mesh API: entities in cell
****************************************************************** */
inline
void cell_get_entities(const AmanziMesh::Mesh& mesh, int c, 
                       const AmanziMesh::Entity_kind kind, 
                       AmanziMesh::Entity_ID_List* entities)
{
  if (kind == AmanziMesh::Entity_kind::FACE) {
    mesh.cell_get_faces(c, entities);
  } else if (kind == AmanziMesh::Entity_kind::EDGE) {
    mesh.cell_get_edges(c, entities);
  } else if (kind == AmanziMesh::Entity_kind::NODE) {
    mesh.cell_get_nodes(c, entities);
  } else if (kind == AmanziMesh::Entity_kind::CELL) {
    entities->clear();
    entities->push_back(c);
  } else {
    entities->clear();
  }
}


/* ******************************************************************
* Extension of Mesh API: neighboor of a cell 
****************************************************************** */
inline
int cell_get_face_adj_cell(const AmanziMesh::Mesh& mesh, int c, int f)
{
  AmanziMesh::Entity_ID_List cells;
  mesh.getFaceCells(f, Parallel_type::ALL, cells);

  if (cells.size() == 2)
    return cells[0] + cells[1] - c;

  return -1;
}


/* ******************************************************************
* Exterior boundary normal: dir = 0 for internal face
****************************************************************** */
inline
AmanziGeometry::Point face_normal_exterior(const AmanziMesh::Mesh& mesh, int f, int* dir)
{
  Amanzi::AmanziMesh::Entity_ID_List cells;
  mesh.getFaceCells(f, Amanzi::AmanziMesh::Parallel_type::ALL, cells);

  auto normal = mesh.getFaceNormal(f,  cells[0], dir);
  if (cells.size() > 1) *dir = 0;

  return normal;
}


/* ******************************************************************
* Geometric center of a mesh cell 
****************************************************************** */
inline
AmanziGeometry::Point cell_geometric_center(const AmanziMesh::Mesh& mesh, int c)
{
  int d = mesh.getSpaceDimension();
  AmanziGeometry::Point v(d), xg(d);

  AmanziMesh::Entity_ID_List nodes;
  mesh.getCellNodes(c, nodes);
  int nnodes = nodes.size();

  for (int i = 0; i < nnodes; ++i) {
    mesh.node_get_coordinates(nodes[i], &v);
    xg += v;
  } 
  xg /= nnodes;
  
  return xg;
}


/* ******************************************************************
* Geometric center of a mesh face
****************************************************************** */
inline
AmanziGeometry::Point face_geometric_center(const AmanziMesh::Mesh& mesh, int f)
{
  int d = mesh.getSpaceDimension();
  AmanziGeometry::Point v(d), xg(d);

  AmanziMesh::Entity_ID_List nodes;
  mesh.getFaceNodes(f, nodes);
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

