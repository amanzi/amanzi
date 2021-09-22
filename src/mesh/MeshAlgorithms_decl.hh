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

class MeshCache;

namespace MeshAlgorithms {

template<class Mesh_type>
Cell_type getCellType(const Mesh_type& mesh, const Entity_ID c);

template<class Mesh_type>
int getFaceDirectionInCell(const Mesh_type& mesh, const Entity_ID f, const Entity_ID c);

//
// topology algorithms
//
template<class Mesh_type>
Entity_ID_List
computeCellEdges(const Mesh_type& mesh, const Entity_ID c);

template<class Mesh_type>
Entity_ID_List
computeCellNodes(const Mesh_type& mesh, const Entity_ID c);


//
// Geometry algorithms
//
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

template<class Mesh_type>
Point_List getEdgeCoordinates(const Mesh_type& mesh, const Entity_ID e);

template<class Mesh_type>
Point_List getFaceCoordinates(const Mesh_type& mesh, const Entity_ID f);

template<class Mesh_type>
Point_List getCellCoordinates(const Mesh_type& mesh, const Entity_ID c);

//
// Deformation algorithms
//

namespace Impl {

// Basic nodal deformation.
//
// NOTE: The user is responsible for ensuring consistency of ghost node
// coordinates -- no communication is done here, so if a node is moved on one
// process, it must be moved on all processes.
//
// NOTE: No geometric consistency is enforced here, so tangling is possible.
//
// NOTE: regions move with the mesh -- deforming a mesh does not cause
// geometric regions to be recalculated.
//
// NOTE: prefer to use deform, not this!
//
template<class Mesh_type>
int
setNodeCoordinates(Mesh_type& mesh,
                   const Entity_ID_List& nodeids,
                   const Point_List& newpos);
}

template<class Mesh_type>
int
deform(Mesh_type& mesh,
       const Entity_ID_List& nodeids,
       const Point_List& newpos);

int
deform(MeshCache& mesh,
       const Entity_ID_List& nodeids,
       const Point_List& newpos);

} // namespace MeshAlgorithms
} // namespace AmanziMesh
} // namspace Amanzi
