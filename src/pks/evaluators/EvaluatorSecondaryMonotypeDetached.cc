/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov
*/

#include "EvaluatorSecondaryMonotypeDetached.hh"

namespace Amanzi {

void
EvaluatorSecondaryMonotypeDetached::UpdateDerivative_(State& S,
                                                      const Key& wrt_key,
                                                      const Tag& wrt_tag)
{
  std::vector<CompositeVector*> results(my_keys_.size());
  int j = 0;
  for (const auto& keytag : my_keys_) {
    results[j] = &S.GetDerivativeW<CompositeVector>(
      keytag.first, keytag.second, wrt_key, wrt_tag, keytag.first);
    results[j]->PutScalarMasterAndGhosted(0.0);
    ++j;
  }

  // if provides key, then the result is 1
  if (ProvidesKey(wrt_key, wrt_tag)) {
    auto keytag = std::make_pair(wrt_key, wrt_tag);
    int i = std::find(my_keys_.begin(), my_keys_.end(), keytag) - my_keys_.begin();
    results[i]->PutScalar(1.0);
    return;
  }

  // dF/dx = sum_(deps) partial F / partial dep * ddep/dx + partial F/partial x
  for (auto& dep : dependencies_) {
    if (std::find(detached_dependencies_.begin(), detached_dependencies_.end(), dep) != detached_dependencies_.end()) continue;

    if (wrt_key == dep.first && wrt_tag == dep.second) {
      // partial F / partial x
      std::vector<CompositeVector> tmp_data(my_keys_.size(), *results[0]);
      std::vector<CompositeVector*> tmp(my_keys_.size());
      for (int i = 0; i != my_keys_.size(); ++i) {
        tmp[i] = &tmp_data[i];
      }
      EvaluatePartialDerivative_(S, wrt_key, wrt_tag, tmp);
      for (int i = 0; i != my_keys_.size() ; ++i) results[i]->Update(1., tmp_data[i], 1.);

    } else if (!S.GetEvaluator(dep.first, dep.second).ProvidesKey(wrt_key, wrt_tag) &&
               S.GetEvaluator(dep.first, dep.second).IsDifferentiableWRT(S, wrt_key, wrt_tag)) {
      // partial F / partial dep * ddep/dx
      // -- ddep/dx
      const auto& ddep = S.GetDerivative<CompositeVector>(dep.first, dep.second, wrt_key, wrt_tag);

      // -- partial F / partial dep
      std::vector<CompositeVector> tmp_data(my_keys_.size(), *results[0]);
      std::vector<CompositeVector*> tmp(my_keys_.size());
      for (int i = 0; i != my_keys_.size(); ++i) {
        tmp[i] = &tmp_data[i];
      }
      EvaluatePartialDerivative_(S, dep.first, dep.second, tmp);

      // sum
      for (int i = 0; i != my_keys_.size() ; ++i) results[i]->Multiply(1., ddep, tmp_data[i], 1.);
    }
  }
}

} // namespace Amanzi
