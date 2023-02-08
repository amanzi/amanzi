/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*
  State

*/

#include "EvaluatorCellVolume.hh"

namespace Amanzi {


void
EvaluatorCellVolume::EnsureCompatibility(State& S)
{
  EvaluatorIndependent<CompositeVector, CompositeVectorSpace>::EnsureCompatibility(S);
  Key domain = Keys::getDomain(my_key_);
  S.Require<CompositeVector, CompositeVectorSpace>(my_key_, my_tag_)
    .SetMesh(S.GetMesh(domain))
    ->SetGhosted()
    ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
}


// ---------------------------------------------------------------------------
// Does the actual work to update the value in the state.
// ---------------------------------------------------------------------------
void
EvaluatorCellVolume::Update_(State& S)
{
  auto& vec = S.GetW<CompositeVector>(my_key_, my_tag_, my_key_);
  for (const auto& comp : vec) {
    if (comp == "cell") {
      int ncells = vec.Mesh()->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
      auto& vec_c = *vec.ViewComponent("cell", false);
      for (int c = 0; c != ncells; ++c) { vec_c[0][c] = vec.Mesh()->getCellVolume(c); }
    } else if (comp == "face") {
      int nfaces = vec.Mesh()->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED);
      auto& vec_f = *vec.ViewComponent("face", false);
      for (int f = 0; f != nfaces; ++f) { vec_f[0][f] = vec.Mesh()->getFaceArea(f); }
    } else {
      Errors::Message message;
      message << "EvaluatorCellVolume for \"" << my_key_ << "\" requires component \"" << comp
              << "\" which is not \"face\" or \"cell\"";
      throw message;
    }
  }
}

} // namespace Amanzi
