/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

#include "Patch.hh"
#include "EvaluatorSecondaryVectorAsPatch.hh"

namespace Amanzi {

const std::string EvaluatorSecondaryVectorAsPatch::eval_type = "vector as patch";

void
EvaluatorSecondaryVectorAsPatch::EnsureCompatibility(State& S)
{
  auto my_key_tag = my_keys_.front();
  S.Require<MultiPatch<double>,MultiPatchSpace>(my_key_tag.first, my_key_tag.second, my_key_tag.first);

  auto dep = dependencies_.front();
  S.Require<CompositeVector,CompositeVectorSpace>(dep.first, dep.second);
}


void
EvaluatorSecondaryVectorAsPatch::Update_(State& S)
{
  S.GetW<MultiPatch<double>>(my_keys_.front().first, my_keys_.front().second, my_keys_.front().first)[0].data =
    S.GetW<CompositeVector>(dependencies_.front().first, dependencies_.front().second, dependencies_.front().first)
    .viewComponent("cell", false);
}

} // namespace Amanzi
