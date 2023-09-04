/*
  Copyright 2010-201x held jointly by LANL, ORNL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/
// Helper functions for resolving regions on meshes.
/*

  Here, we use the word "resolve" to mean figuring out which entities are in
  the set.  The word is not particularly clear, but it is used in the sense of:
  "resolving whether an entity is inside or outside of the set."

*/

#pragma once

#include <map>

#include "Region.hh"
#include "RegionEnumerated.hh"
#include "RegionLabeledSet.hh"
#include "RegionLogical.hh"
#include "RegionPoint.hh"
#include "MeshDefs.hh"
//#include "Mesh_Algorithms.hh"
#include "MeshMaps.hh"

namespace Amanzi {
namespace AmanziMesh {

template <MemSpace_kind>
struct MeshCache;

//
// This casts as needed to find the right implementation.  Note that some
// require casting, some do not, so we avoid if possible.
//
View_type<const Entity_ID, MemSpace_kind::HOST>
resolveMeshSet(const AmanziGeometry::Region& region,
               const Entity_kind kind,
               const Parallel_kind ptype,
               const MeshCache<MemSpace_kind::HOST>& mesh);

View_type<const Entity_ID, MemSpace_kind::HOST>
resolveMeshSetVolumeFractions(const AmanziGeometry::Region& region,
                              const Entity_kind kind,
                              const Parallel_kind ptype,
                              Double_View& vol_fracs,
                              const MeshCache<MemSpace_kind::HOST>& mesh);


namespace Impl {

View_type<const Entity_ID, MemSpace_kind::HOST>
resolveMeshSet_(const AmanziGeometry::Region& region,
                const Entity_kind kind,
                const Parallel_kind ptype,
                const MeshCache<MemSpace_kind::HOST>& mesh);

View_type<const Entity_ID, MemSpace_kind::HOST>
resolveBoundaryEntityMeshSet(const AmanziGeometry::Region& region,
                             const Entity_kind kind,
                             const Parallel_kind ptype,
                             const MeshCache<MemSpace_kind::HOST>& parent_mesh);

View_type<const Entity_ID, MemSpace_kind::HOST>
resolveIDMeshSetFromParent(const AmanziGeometry::Region& region,
                           const Entity_kind kind,
                           const Parallel_kind ptype,
                           const MeshCache<MemSpace_kind::HOST>& mesh,
                           const MeshCache<MemSpace_kind::HOST>& parent_mesh);

View_type<const Entity_ID, MemSpace_kind::HOST>
resolveGeometricMeshSetFromParent(const AmanziGeometry::Region& region,
                                  const Entity_kind kind,
                                  const Parallel_kind ptype,
                                  const MeshCache<MemSpace_kind::HOST>& mesh,
                                  const MeshCache<MemSpace_kind::HOST>& parent_mesh);

View_type<const Entity_ID, MemSpace_kind::HOST>
resolveMeshSetAll(const AmanziGeometry::Region& region,
                  const Entity_kind kind,
                  const Parallel_kind ptype,
                  const MeshCache<MemSpace_kind::HOST>& mesh);

View_type<const Entity_ID, MemSpace_kind::HOST>
resolveMeshSetBoundary(const AmanziGeometry::Region& region,
                       const Entity_kind kind,
                       const Parallel_kind ptype,
                       const MeshCache<MemSpace_kind::HOST>& mesh);

View_type<const Entity_ID, MemSpace_kind::HOST>
resolveMeshSetEnumerated(const AmanziGeometry::RegionEnumerated& region,
                         const Entity_kind kind,
                         const Parallel_kind ptype,
                         const MeshCache<MemSpace_kind::HOST>& mesh);

View_type<const Entity_ID, MemSpace_kind::HOST>
resolveMeshSetGeometric(const AmanziGeometry::Region& region,
                        const Entity_kind kind,
                        const Parallel_kind ptype,
                        const MeshCache<MemSpace_kind::HOST>& mesh);

View_type<const Entity_ID, MemSpace_kind::HOST>
resolveMeshSetLabeledSet(const AmanziGeometry::RegionLabeledSet& region,
                         const Entity_kind kind,
                         const Parallel_kind ptype,
                         const MeshCache<MemSpace_kind::HOST>& mesh);

View_type<const Entity_ID, MemSpace_kind::HOST>
resolveMeshSetLogical(const AmanziGeometry::RegionLogical& region,
                      const Entity_kind kind,
                      const Parallel_kind ptype,
                      const MeshCache<MemSpace_kind::HOST>& mesh);

cEntity_ID_View
resolveMeshSetPoint(const AmanziGeometry::RegionPoint& region,
                    const Entity_kind kind,
                    const Parallel_kind ptype,
                    const MeshCache<MemSpace_kind::HOST>& mesh);

View_type<const Entity_ID, MemSpace_kind::HOST>
filterParentEntities(const MeshCache<MemSpace_kind::HOST>& mesh,
                     Entity_kind kind,
                     Parallel_kind ptype,
                     const View_type<const Entity_ID, MemSpace_kind::HOST>& parent_entities);

View_type<const Entity_ID, MemSpace_kind::HOST>
filterParentEntities_SurfaceCellToCell(
  const MeshCache<MemSpace_kind::HOST>& mesh,
  Parallel_kind ptype,
  const View_type<const Entity_ID, MemSpace_kind::HOST>& parent_entities);

View_type<const Entity_ID, MemSpace_kind::HOST>
filterParentEntities_SurfaceFaceToFace(
  const MeshCache<MemSpace_kind::HOST>& mesh,
  Parallel_kind ptype,
  const View_type<const Entity_ID, MemSpace_kind::HOST>& parent_entities);

} // namespace Impl
} // namespace AmanziMesh
} // namespace Amanzi
