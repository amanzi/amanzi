/*
  Copyright 2010-201x held jointly by LANL, ORNL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

// Generates MeshCache objects for use in testing.

#pragma once

#include "MeshFrameworkTraits.hh"
#include "MeshCache.hh"

using namespace Amanzi;
using namespace AmanziMesh;

inline
Teuchos::RCP<MeshCache>
createMeshFromFramework(const Teuchos::RCP<MeshFramework>& mesh) {
  auto mc = Teuchos::rcp(new MeshCache(mesh));
  Utils::cacheAll(*mc);
  return mc;
}
