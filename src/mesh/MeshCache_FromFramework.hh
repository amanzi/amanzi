/*
  Copyright 2010-201x held jointly by LANL, ORNL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
           Julien Loiseau (jloiseau@lanl.gov)
*/
//! Builds a MeshCache from a MeshFramework

#pragma once

#include "MeshFramework.hh"
#include "MeshCache.hh"

namespace Amanzi {
namespace AmanziMesh {

class MeshCache_FromFramework : public MeshCache {
 public:
  MeshCache_FromFramework(const Teuchos::RCP<MeshFramework>& framework_mesh)
    : framework_mesh_(framework_mesh) {
    Initialize();
  }

  void Finalize() {
    framework_mesh_ = Teuchos::null;
  }

  void Initialize();

}


} // namespace AmanziMesh
} // namespace Amanzi
