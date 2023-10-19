/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Rao Garimella
           Ethan Coon (coonet@ornl.gov)
*/

// A collection of helper functions used to compute basic geometric operations on Mesh objects.
/*
  Note, these are templated on Mesh_type because they may be used by either
  MeshFramework objects or MeshCache objects.

  These provide default calculations of fundamental topological and geometric
  concepts.  While they can be called by users of the Mesh library, typically
  they are used _inside_ of the library, and should not be called by users.

  For the most part, these are used by member functions of the MeshFramework
  and MeshCache classes to provide default implementations.  They are
  implemented as external, nonmember functions because, 1, they use the common
  API of multiple classes, and so are useful for multiple classes
  (MeshFramework and MeshCache), and 2, they use the API, and need no "private"
  data, so they improve encapsulation of the class.

  Note that, if these become commonly used by users of the mesh library, they
  should get moved into Mesh_Algorithms, which provides functions that are
  commonly used by clients _outside_ of this library.
*/


#pragma once

#include "Point.hh"
#include "MeshDefs.hh"

namespace Amanzi {
namespace AmanziMesh {
namespace MeshAlgorithms {

//
// topology algorithms
//
template <class Mesh_type>
KOKKOS_INLINE_FUNCTION Cell_kind
getCellType(const Mesh_type& mesh, const Entity_ID c);

template <class Mesh_type>
int
getFaceDirectionInCell(const Mesh_type& mesh, const Entity_ID f, const Entity_ID c);

template <class Mesh_type>
cEntity_ID_View
computeCellEdges(const Mesh_type& mesh, const Entity_ID c);

template <class Mesh_type>
cEntity_ID_View
computeCellNodes(const Mesh_type& mesh, const Entity_ID c);

template <class Mesh_type>
cEntity_ID_View
computeNodeCells(const Mesh_type& mesh, const Entity_ID n);


//
// Geometry algorithms
//
template <class Mesh_type>
std::pair<double, AmanziGeometry::Point>
computeCellGeometry(const Mesh_type& mesh, const Entity_ID c);

template <class Mesh_type>
std::tuple<double, AmanziGeometry::Point, cPoint_View>
computeFaceGeometry(const Mesh_type& mesh, const Entity_ID f);

template <class Mesh_type>
std::pair<AmanziGeometry::Point, AmanziGeometry::Point>
computeEdgeGeometry(const Mesh_type& mesh, const Entity_ID e);

template <class Mesh_type>
Point_View
computeBisectors(const Mesh_type& mesh, const Entity_ID c, const Entity_ID_View& faces);

template <class Mesh_type>
void
debugCell(const Mesh_type& mesh, const Entity_ID c);

template <class Mesh_type>
cPoint_View
getEdgeCoordinates(const Mesh_type& mesh, const Entity_ID e);

template <class Mesh_type>
cPoint_View
getFaceCoordinates(const Mesh_type& mesh, const Entity_ID f);

template <class Mesh_type>
cPoint_View
getCellCoordinates(const Mesh_type& mesh, const Entity_ID c);

template <class Mesh_type>
std::size_t
getMaxCellNumNodes(const Mesh_type& mesh);

template <class Mesh_type>
std::size_t
getMaxCellNumFaces(const Mesh_type& mesh);

template <class Mesh_type>
std::size_t
getMaxCellNumEdges(const Mesh_type& mesh);

template <class Mesh_type>
Entity_ID_View
getCellFaceAdjacentCells(const Mesh_type& mesh, Entity_ID c, Parallel_kind ptype);

template <class Mesh_type>
int
getFaceAdjacentCell(const Mesh_type& mesh, Entity_ID c, Entity_ID f);

template <class Mesh_type>
KOKKOS_INLINE_FUNCTION AmanziGeometry::Point
getFaceCentroid(const Mesh_type& mesh, const Entity_ID f);

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
template <class Mesh_type>
int
setNodeCoordinates(Mesh_type& mesh, const Entity_ID_View& nodeids, const Point_View& newpos);
} // namespace Impl

template <class Mesh_type>
int
deform(Mesh_type& mesh, const Entity_ID_View& nodeids, const Point_View& newpos);

} // namespace MeshAlgorithms
} // namespace AmanziMesh
} // namespace Amanzi
