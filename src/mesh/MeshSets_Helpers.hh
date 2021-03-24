/*
  Copyright 2010-201x held jointly by LANL, ORNL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/
// Helper functions for resolving regions on meshes.

#pragma once

#include "RegionLogical.hh"
#include "RegionEnumerated.hh"

namespace Amanzi {
namespace AmanziMesh {

class MeshCache;

//
// This casts as needed to find the right implementation.
//
Entity_ID_List resolveMeshSet(const AmanziGeometry::Region& region,
        const Entity_kind kind,
        const Parallel_type ptype,
        const MeshCache& mesh);

//
// The default implementation assumes that region is a geometric region.
//
Entity_ID_List resolveMeshSetGeometric(const AmanziGeometry::Region& region,
        const Entity_kind kind,
        const Parallel_type ptype,
        const MeshCache& mesh);

//
// Implements RegionLogical
//
Entity_ID_List resolveMeshSet(const AmanziGeometry::RegionLogical& region,
        const Entity_kind kind,
        const Parallel_type ptype,
        const MeshCache& mesh);

//
// Implements RegionEnumerated
//
Entity_ID_List resolveMeshSet(const AmanziGeometry::RegionEnumerated& region,
        const Entity_kind kind,
        const Parallel_type ptype,
        const MeshCache& mesh);


} // namespace AmanziMesh
} // namespace Amanzi
