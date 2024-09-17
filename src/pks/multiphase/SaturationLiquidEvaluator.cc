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

#include "SaturationLiquidEvaluator.hh"

namespace Amanzi {
namespace Multiphase {

/* ******************************************************************
* Simple constructor
****************************************************************** */
SaturationLiquidEvaluator::SaturationLiquidEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(plist)
{
  if (my_keys_.size() == 0) {
    my_keys_.push_back(std::make_pair(plist_.get<std::string>("my key"), Tags::DEFAULT));
  }
  saturation_gas_key_ = plist_.get<std::string>("saturation gas key");
  dependencies_.insert(std::make_pair(saturation_gas_key_, Tags::DEFAULT));
}


/* ******************************************************************
* Copy constructor.
****************************************************************** */
Teuchos::RCP<Evaluator>
SaturationLiquidEvaluator::Clone() const
{
  return Teuchos::rcp(new SaturationLiquidEvaluator(*this));
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void
SaturationLiquidEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& results)
{
  const auto& sg = *S.Get<CompositeVector>(saturation_gas_key_).ViewComponent("cell");
  auto& result_c = *results[0]->ViewComponent("cell");

  int ncells = result_c.MyLength();
  for (int c = 0; c != ncells; ++c) { result_c[0][c] = 1.0 - sg[0][c]; }
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void
SaturationLiquidEvaluator::EvaluatePartialDerivative_(const State& S,
                                                      const Key& wrt_key,
                                                      const Tag& wrt_tag,
                                                      const std::vector<CompositeVector*>& results)
{
  auto& result_c = *results[0]->ViewComponent("cell");

  int ncells = result_c.MyLength();
  if (wrt_key == saturation_gas_key_) {
    for (int c = 0; c != ncells; ++c) result_c[0][c] = -1.0;
  }
}

} // namespace Multiphase
} // namespace Amanzi
