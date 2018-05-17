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
#include "Point.hh"

namespace Amanzi {
namespace WhetStone {

inline
AmanziGeometry::Point cell_geometric_center(const AmanziMesh::Mesh& mesh, int c)
{
  int d = mesh.space_dimension();
  AmanziGeometry::Point v(d), xg(d);

  Entity_ID_List nodes;
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

