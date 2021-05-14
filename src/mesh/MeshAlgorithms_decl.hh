/*
  Copyright 2010-201x held jointly by LANL, ORNL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Rao Garimella
      Ethan Coon (coonet@ornl.gov)
*/

#pragma once

#include "Point.hh"
#include "MeshDefs.hh"

namespace Amanzi {
namespace AmanziMesh {
namespace MeshAlgorithms {

template<class Mesh_type>
Cell_type getCellType(const Mesh_type& mesh, const Entity_ID c);

template<class Mesh_type>
int getFaceDirectionInCell(const Mesh_type& mesh, const Entity_ID f, const Entity_ID c);

template<class Mesh_type>
std::pair<double,AmanziGeometry::Point>
computeCellGeometry(const Mesh_type& mesh, const Entity_ID c);

template<class Mesh_type>
std::tuple<double,AmanziGeometry::Point,Point_List>
computeFaceGeometry(const Mesh_type& mesh, const Entity_ID f);

template<class Mesh_type>
std::pair<AmanziGeometry::Point,AmanziGeometry::Point>
computeEdgeGeometry(const Mesh_type& mesh, const Entity_ID e);

template<class Mesh_type>
void computeBisectors(const Mesh_type& mesh, const Entity_ID c,
        const Entity_ID_List& faces, Point_List& bisectors);

template<class Mesh_type>
void debugCell(const Mesh_type& mesh, const Entity_ID c);


} // namespace MeshAlgorithms
} // namespace AmanziMesh
} // namspace Amanzi
