/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/
/*

Internally used functions that can be used by both MeshCache and MeshFramework.
These are implemented using the remainder of the Mesh API.

*/

#pragma once

#include "MeshDefs.hh"

namespace Amanzi {
namespace AmanziMesh {
namespace Impl {

//
// Internal-only, non-virtual algorithms -- not for use by client code?
//

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
std::vector<int>
mapFaceToCellEdges(const Mesh_type& mesh, const Entity_ID f, const Entity_ID c);

template <class Mesh_type>
typename Mesh_type::Entity_ID_View
computeCellEdges(const Mesh_type& mesh, const Entity_ID c);

template <class Mesh_type>
typename Mesh_type::Entity_ID_View
computeCellNodes(const Mesh_type& mesh, const Entity_ID c);

template <class Mesh_type>
typename Mesh_type::Entity_ID_View
computeNodeCells(const Mesh_type& mesh, const Entity_ID n);


//
// Geometry algorithms
//
template <class Mesh_type>
std::pair<double, AmanziGeometry::Point>
computeCellGeometry(const Mesh_type& mesh, const Entity_ID c);

template <class Mesh_type>
std::tuple<double, AmanziGeometry::Point, typename Mesh_type::Point_View>
computeFaceGeometry(const Mesh_type& mesh, const Entity_ID f);

template <class Mesh_type>
std::pair<AmanziGeometry::Point, AmanziGeometry::Point>
computeEdgeGeometry(const Mesh_type& mesh, const Entity_ID e);

template <class Mesh_type>
typename Mesh_type::Point_View
computeBisectors(const Mesh_type& mesh,
                 const Entity_ID c,
                 const typename Mesh_type::cEntity_ID_View& faces);

template <class Mesh_type>
void
debugCell(const Mesh_type& mesh, const Entity_ID c);

template <class Mesh_type>
typename Mesh_type::Point_View
computeEdgeCoordinates(const Mesh_type& mesh, const Entity_ID e);

template <class Mesh_type>
typename Mesh_type::Point_View
computeFaceCoordinates(const Mesh_type& mesh, const Entity_ID f);

template <class Mesh_type>
typename Mesh_type::Point_View
computeCellCoordinates(const Mesh_type& mesh, const Entity_ID c);

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
KOKKOS_INLINE_FUNCTION AmanziGeometry::Point
getFaceCentroid(const Mesh_type& mesh, const Entity_ID f);

} // namespace Impl
} // namespace AmanziMesh
} // namespace Amanzi
