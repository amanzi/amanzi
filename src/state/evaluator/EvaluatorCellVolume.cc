/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

//!

/*
  State

  License: BSD
  Author: Ethan Coon

*/

#include "EvaluatorCellVolume.hh"

namespace Amanzi {

// ---------------------------------------------------------------------------
// Does the actual work to update the value in the state.
// ---------------------------------------------------------------------------
void
EvaluatorCellVolume::Update_(State& S)
{
  auto& vec = S.GetW<CompositeVector>(my_key_, my_tag_, my_key_);
  const AmanziMesh::Mesh* mesh = vec.Mesh().get();
  for (const auto& comp : vec) {
    if (comp == "cell") {
      auto vec_c = vec.ViewComponent("cell", false);
      Impl::copyCellVolume(mesh, vec_c);
    } else if (comp == "face") {
      auto vec_c = vec.ViewComponent("face", false);
      Impl::copyFaceArea(mesh, vec_c);
    } else {
      Errors::Message message;
      message << "EvaluatorCellVolume for \"" << my_key_
              << "\" requires component \"" << comp
              << "\" which is not \"face\" or \"cell\"";
      throw message;
    }
  }
}

} // namespace Amanzi
