/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
           Julien Loiseau (jloiseau@lanl.gov)
*/

#pragma once

#include "MeshCache.hh"

namespace Amanzi {
namespace AmanziMesh {

using Mesh = MeshCache<MemSpace_kind::HOST>;

} // namespace AmanziMesh
} // namespace Amanzi
