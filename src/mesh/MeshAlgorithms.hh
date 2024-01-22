/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/
/*

MeshAlgorithms are algorithms that are specific to the MeshFramework
implementation.  They extend the MeshFramework API and provide the COMPUTE
algorithm for the MeshCache implementaiton.  They are kept separate from
the MeshFramework API because then the MeshFramework API can be deleted and
these can be used with the MeshCache instead.

Note that the "virtual" algorithms are kept in a separate class to enable
overriding them for special meshes (MeshLogical).

*/

#pragma once

#include "MeshInternals_decl.hh"
#include "MeshCache_decl.hh"

namespace Amanzi {
namespace AmanziMesh {

//
// This class provides virtual algorithms, with defaults for computing
// geometric quantities given nodal coordinates and topological information.
//
struct MeshAlgorithms {
  virtual ~MeshAlgorithms() = default;

  // lumped things for more efficient calculation
  // KOKKOS_INLINE_FUNCTION
  // virtual std::pair<double, AmanziGeometry::Point>
  // computeCellGeometry(const Mesh& mesh, const Entity_ID c) const;

  virtual std::pair<double, AmanziGeometry::Point>
  computeCellGeometry(const MeshHost& mesh, const Entity_ID c) const {
    return Impl::computeCellGeometry(mesh, c);
  }

  // KOKKOS_INLINE_FUNCTION
  // virtual std::tuple<double, AmanziGeometry::Point, typename Mesh::cPoint_View>
  // computeFaceGeometry(const Mesh& mesh, const Entity_ID f) const;

  virtual std::tuple<double, AmanziGeometry::Point, typename Mesh::cPoint_View>
  computeFaceGeometry(const MeshHost& mesh, const Entity_ID f) const {
    return Impl::computeFaceGeometry(mesh, f);
  }

  // KOKKOS_INLINE_FUNCTION
  // virtual std::pair<AmanziGeometry::Point, AmanziGeometry::Point>
  // computeEdgeGeometry(const Mesh& mesh, const Entity_ID e) const;

  virtual std::pair<AmanziGeometry::Point, AmanziGeometry::Point>
  computeEdgeGeometry(const MeshHost& mesh, const Entity_ID e) const {
    return Impl::computeEdgeGeometry(mesh, e);
  }

  // KOKKOS_INLINE_FUNCTION
  // virtual double computeCellVolume(const Mesh& mesh, const Entity_ID c) const;

  virtual double computeCellVolume(const MeshHost& mesh, const Entity_ID c) const {
    return Impl::computeCellGeometry(mesh, c).first;
  }

  // KOKKOS_INLINE_FUNCTION
  // virtual AmanziGeometry::Point computeCellCentroid(const Mesh& mesh, const Entity_ID c) const;

  virtual AmanziGeometry::Point
  computeCellCentroid(const MeshHost& mesh, const Entity_ID c) const {
    return Impl::computeCellGeometry(mesh, c).second;
  }

  // KOKKOS_INLINE_FUNCTION
  // virtual double computeFaceArea(const Mesh& mesh, const Entity_ID f) const;

  virtual double computeFaceArea(const MeshHost& mesh, const Entity_ID f) const {
    return std::get<0>(Impl::computeFaceGeometry(mesh, f));
  }

  // KOKKOS_INLINE_FUNCTION
  // virtual AmanziGeometry::Point computeFaceCentroid(const Mesh& mesh, const Entity_ID f) const;

  virtual AmanziGeometry::Point
  computeFaceCentroid(const MeshHost& mesh, const Entity_ID f) const {
    return std::get<1>(Impl::computeFaceGeometry(mesh, f));
  }

  // KOKKOS_INLINE_FUNCTION
  // virtual AmanziGeometry::Point computeFaceNormal(const Mesh& mesh, const Entity_ID f,
  //         const Entity_ID c, int * const orientation) const;

  virtual AmanziGeometry::Point computeFaceNormal(const MeshHost& mesh,
                                                         const Entity_ID f,
                                                         const Entity_ID c,
                                                         int* const orientation) const;

  // KOKKOS_INLINE_FUNCTION
  // virtual double computeEdgeLength(const Mesh& mesh, const Entity_ID e) const;

  virtual double computeEdgeLength(const MeshHost& mesh, const Entity_ID e) const {
    return AmanziGeometry::norm(Impl::computeEdgeGeometry(mesh, e).first);
  }

  // KOKKOS_INLINE_FUNCTION
  // virtual AmanziGeometry::Point
  // computeEdgeVector(const Mesh& mesh, const Entity_ID e, const Entity_ID n, int * const orientation) const;

  virtual AmanziGeometry::Point computeEdgeVector(const MeshHost& mesh,
                                                         const Entity_ID e,
                                                         const Entity_ID n,
                                                         int* const orientation) const;

  // KOKKOS_INLINE_FUNCTION
  // virtual AmanziGeometry::Point
  // computeEdgeCentroid(const Mesh& mesh, const Entity_ID e) const;

  virtual AmanziGeometry::Point
  computeEdgeCentroid(const MeshHost& mesh, const Entity_ID e) const {
    return Impl::computeEdgeGeometry(mesh, e).second;
  }

  // KOKKOS_INLINE_FUNCTION
  // virtual void computeCellFacesAndBisectors(const Mesh& mesh, const Entity_ID cellid,
  //         typename Mesh::cEntity_ID_View& faceids, typename Mesh::cPoint_View * const bisectors) const;

  virtual void
  computeCellFacesAndBisectors(const MeshHost& mesh,
                               const Entity_ID cellid,
                               typename MeshHost::cEntity_ID_View& faceids,
                               typename MeshHost::cPoint_View* const bisectors) const;
};



} // namespace AmanziMesh
} // namespace Amanzi
