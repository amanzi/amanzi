/*
  Copyright 2010-201x held jointly by LANL, ORNL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/
// Helper functions for resolving regions on meshes.

#pragma once

#include <map>

#include "Region.hh"
#include "RegionEnumerated.hh"
#include "RegionLabeledSet.hh"
#include "RegionLogical.hh"
#include "MeshDefs.hh"

namespace Amanzi {
namespace AmanziMesh {

using MeshSets = std::map<std::tuple<std::string,Entity_kind,Parallel_type>, Entity_ID_View>;
using MeshSetVolumeFractions = std::map<std::tuple<std::string,Entity_kind,Parallel_type>, Double_View>;


class MeshCache;

//
// This casts as needed to find the right implementation.  Note that some
// require casting, some do not, so we avoid if possible.
//
Entity_ID_List resolveMeshSet(const AmanziGeometry::Region& region,
        const Entity_kind kind,
        const Parallel_type ptype,
        const MeshCache& mesh);

Entity_ID_List resolveMeshSetVolumeFractions(
  const AmanziGeometry::Region& region,
  const Entity_kind kind,
  const Parallel_type ptype,
  Double_View& vol_fracs,
  const MeshCache& mesh);


namespace Impl {

Entity_ID_List resolveMeshSet_(const AmanziGeometry::Region& region,
        const Entity_kind kind,
        const Parallel_type ptype,
        const MeshCache& mesh);

Entity_ID_List
resolveBoundaryEntityMeshSet(const AmanziGeometry::Region& region,
                             const Entity_kind kind,
                             const Parallel_type ptype,
                             const MeshCache& parent_mesh);

Entity_ID_List
resolveIDMeshSetFromParent(const AmanziGeometry::Region& region,
                           const Entity_kind kind,
                           const Parallel_type ptype,
                           const MeshCache& mesh,
                           const MeshCache& parent_mesh);

Entity_ID_List
resolveGeometricMeshSetFromParent(const AmanziGeometry::Region& region,
        const Entity_kind kind,
        const Parallel_type ptype,
        const MeshCache& mesh,
        const MeshCache& parent_mesh);

Entity_ID_List resolveMeshSetAll(const AmanziGeometry::Region& region,
        const Entity_kind kind,
        const Parallel_type ptype,
        const MeshCache& mesh);

Entity_ID_List resolveMeshSetBoundary(const AmanziGeometry::Region& region,
        const Entity_kind kind,
        const Parallel_type ptype,
        const MeshCache& mesh);

Entity_ID_List resolveMeshSetEnumerated(const AmanziGeometry::RegionEnumerated& region,
        const Entity_kind kind,
        const Parallel_type ptype,
        const MeshCache& mesh);

Entity_ID_List resolveMeshSetGeometric(const AmanziGeometry::Region& region,
        const Entity_kind kind,
        const Parallel_type ptype,
        const MeshCache& mesh);

Entity_ID_List resolveMeshSetLabeledSet(const AmanziGeometry::RegionLabeledSet& region,
        const Entity_kind kind,
        const Parallel_type ptype,
        const MeshCache& mesh);

Entity_ID_List resolveMeshSetLogical(const AmanziGeometry::RegionLogical& region,
        const Entity_kind kind,
        const Parallel_type ptype,
        const MeshCache& mesh);


Entity_ID_List filterParentEntities(const MeshCache& mesh,
        Entity_kind kind,
        Parallel_type ptype,
        const Entity_ID_List& parent_entities);

Entity_ID_List filterParentEntities_SurfaceCellToCell(const MeshCache& mesh,
        Parallel_type ptype,
        const Entity_ID_List& parent_entities);

Entity_ID_List filterParentEntities_SurfaceFaceToFace(const MeshCache& mesh,
        Parallel_type ptype,
        const Entity_ID_List& parent_entities);

} // namespace Impl
} // namespace AmanziMesh
} // namespace Amanzi


