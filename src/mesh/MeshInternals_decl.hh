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
#include "MeshCache_decl.hh"

namespace Amanzi {
namespace AmanziMesh {
namespace Impl {

//
// Internal-only, non-virtual algorithms -- not for use by client code?
//

//
// topology algorithms
//
DISABLE_CUDA_WARNING
template <class Mesh_type>
KOKKOS_INLINE_FUNCTION Cell_kind
getCellType(const Mesh_type& mesh, const Entity_ID c)
{
  auto faces = mesh.getCellFaces(c);
  if (mesh.getManifoldDimension() == 2) {
    switch (faces.size()) {
    case 3:
      return Cell_kind::TRI;
      break;
    case 4:
      return Cell_kind::QUAD;
      break;
    default:
      return Cell_kind::POLYGON;
    }
  } else if (mesh.getManifoldDimension() == 3) {
    int nquads = 0;
    for (const auto& f : faces) {
      typename Mesh_type::cEntity_ID_View fnodes;
      mesh.getFaceNodes(f, fnodes);
      if (fnodes.size() == 4) nquads++;
    }

    switch (faces.size()) {
    case 4:
      if (nquads == 0)
        return Cell_kind::TET;
      else
        return Cell_kind::POLYHED;
      break;
    case 5:
      if (nquads == 1)
        return Cell_kind::PYRAMID;
      else if (nquads == 3)
        return Cell_kind::PRISM;
      else
        return Cell_kind::POLYHED;
      break;
    case 6:
      if (nquads == 6)
        return Cell_kind::HEX;
      else
        return Cell_kind::POLYHED;
      break;
    default:
      return Cell_kind::POLYHED;
    }
  } else {
    if (!std::is_same_v<Mesh_type, AmanziMesh::MeshCacheDevice>) {
      assert(false && "Mesh not supported"); 
      //Errors::Message msg;
      //msg << "Mesh of manifold_dimension = " << mesh.getManifoldDimension() << " not supported";
      //Exceptions::amanzi_throw(msg);
    } else {
      assert(false); // "Invalid mesh manifold dimension");
    }
  }
  return Cell_kind::UNKNOWN;
}


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

template <class Mesh_type>
typename Mesh_type::Entity_ID_View
computeEdgeCells(const Mesh_type& mesh, const Entity_ID n);

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
