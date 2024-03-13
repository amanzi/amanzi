/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#pragma once

#include "MeshFramework.hh"
#include "MeshAudit_decl.hh"
#include "MeshAudit_impl.hh"

namespace Amanzi {
namespace AmanziMesh {

using MeshFrameworkAudit = MeshAudit_<MeshFramework, Impl::MeshAudit_Geometry>;

}
} // namespace Amanzi
