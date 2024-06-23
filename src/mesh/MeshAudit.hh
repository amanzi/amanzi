/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#pragma once

#include "Mesh.hh"
#include "MeshCache.hh"
#include "MeshAudit_decl.hh"
#include "MeshAudit_impl.hh"

namespace Amanzi {
namespace AmanziMesh {

using MeshAuditCache = MeshAudit_<Mesh, MeshCache, Impl::MeshAudit_Sets>;
using MeshAudit = MeshAudit_<Mesh, Mesh const *, Impl::MeshAudit_Sets>;

} // namespace AmanziMesh
} // namespace Amanzi
