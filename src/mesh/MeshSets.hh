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

#include "MeshDefs.hh"

namespace Amanzi {
namespace AmanziMesh {

using MeshSets = std::map<std::tuple<std::string,Entity_kind,Parallel_type>, Entity_ID_View>;
using MeshSetVolumeFractions = std::map<std::tuple<std::string,Entity_kind,Parallel_type>, Double_View>;

}
}

