/*
  Shallow Water PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include "errors.hh"

#include "DischargeEvaluator.hh"

namespace Amanzi {
namespace ShallowWater {

/* ******************************************************************
* Simple constructor
****************************************************************** */
DischargeEvaluator::DischargeEvaluator(Teuchos::ParameterList& plist)
    : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(plist)
{
  if (my_keys_.size() == 0) {
    my_keys_.push_back(std::make_pair(plist_.get<std::string>("my key", "surface-discharge"), Tags::DEFAULT));
  }
  std::string domain = Keys::getDomain(my_keys_[0].first);

  ponded_depth_key_ = plist_.get<std::string>("ponded depth key", Keys::getKey(domain, "ponded_depth"));
  velocity_key_ = plist_.get<std::string>("velocity key", Keys::getKey(domain, "velocity"));

  dependencies_.insert(std::make_pair(ponded_depth_key_, Tags::DEFAULT));
  dependencies_.insert(std::make_pair(velocity_key_, Tags::DEFAULT));
}


/* ******************************************************************
* Copy constructor.
****************************************************************** */
Teuchos::RCP<Evaluator> DischargeEvaluator::Clone() const {
  return Teuchos::rcp(new DischargeEvaluator(*this));
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void DischargeEvaluator::Evaluate_(
    const State& S, const std::vector<CompositeVector*>& results)
{
  const auto& h_c = *S.Get<CompositeVector>(ponded_depth_key_).ViewComponent("cell");
  const auto& u_c = *S.Get<CompositeVector>(velocity_key_).ViewComponent("cell");
  auto& result_c = *results[0]->ViewComponent("cell");

  int ncells = result_c.MyLength();
  for (int c = 0; c != ncells; ++c) {
    for (int i = 0; i < 2; ++i) {
      result_c[i][c] = h_c[0][c] * u_c[i][c];
    }
  }
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void DischargeEvaluator::EvaluatePartialDerivative_(
    const State& S, const Key& wrt_key, const Tag& wrt_tag,
    const std::vector<CompositeVector*>& results) 
{
  const auto& h_c = *S.Get<CompositeVector>(ponded_depth_key_).ViewComponent("cell");
  const auto& u_c = *S.Get<CompositeVector>(velocity_key_).ViewComponent("cell");
  auto& result_c = *results[0]->ViewComponent("cell");

  int ncells = result_c.MyLength();
  if (wrt_key == ponded_depth_key_) {
    for (int c = 0; c != ncells; ++c) {
      for (int i = 0; i < 2; ++i) {
        result_c[i][c] = u_c[i][c];
      }
    }
  } else {
    AMANZI_ASSERT(false);
  }
}

}  // namespace ShallowWater
}  // namespace Amanzi

