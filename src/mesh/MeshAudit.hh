/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.
*/

#pragma once

#include "Mesh.hh"
#include "MeshAudit_decl.hh"
#include "MeshAudit_impl.hh"

namespace Amanzi {
namespace AmanziMesh {

using MeshAudit = MeshAudit_<Mesh, Impl::MeshAudit_Sets>;

}
}
