/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  MultiPhase PK

*/

#include "SaturationEvaluator.hh"

namespace Amanzi {
namespace Multiphase {

/* ******************************************************************
* Simple constructor
****************************************************************** */
SaturationEvaluator::SaturationEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(plist)
{
  if (my_keys_.size() == 0) {
    my_keys_.push_back(std::make_pair(plist_.get<std::string>("my key"), Tags::DEFAULT));
  }
  saturation_complement_key_ = plist_.get<std::string>("saturation complement key");
  dependencies_.insert(std::make_pair(saturation_complement_key_, Tags::DEFAULT));
}


/* ******************************************************************
* Copy constructor.
****************************************************************** */
Teuchos::RCP<Evaluator>
SaturationEvaluator::Clone() const
{
  return Teuchos::rcp(new SaturationEvaluator(*this));
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void
SaturationEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& results)
{
  const auto& sl = *S.Get<CompositeVector>(saturation_complement_key_).ViewComponent("cell");
  auto& result_c = *results[0]->ViewComponent("cell");

  int ncells = result_c.MyLength();
  for (int c = 0; c != ncells; ++c) {
    result_c[0][c] = 1.0 - sl[0][c];
  }
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void
SaturationEvaluator::EvaluatePartialDerivative_(const State& S,
                                                const Key& wrt_key,
                                                const Tag& wrt_tag,
                                                const std::vector<CompositeVector*>& results)
{
  auto& result_c = *results[0]->ViewComponent("cell");

  int ncells = result_c.MyLength();
  if (wrt_key == saturation_complement_key_) {
    for (int c = 0; c != ncells; ++c) result_c[0][c] = -1.0;
  }
}

} // namespace Multiphase
} // namespace Amanzi
