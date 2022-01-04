/*
  Multiphase PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
           Konstantin Lipnikov (lipnikov@lanl.gov)

  Relative permeability as a function of liquid saturation.
*/

#include "MultiphaseDefs.hh"
#include "MultiphaseTypeDefs.hh"
#include "RelPermEvaluator.hh"

namespace Amanzi {
namespace Multiphase {

/* ******************************************************************
* Two constructors.
****************************************************************** */
RelPermEvaluator::RelPermEvaluator(Teuchos::ParameterList& plist,
                                   const Teuchos::RCP<WRMmpPartition>& wrm)
    : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(plist),
      wrm_(wrm)
{
  if (my_keys_.size() == 0) {
    my_keys_.push_back(std::make_pair(plist_.get<std::string>("my key"), Tags::DEFAULT));
  }
  std::string domain = Keys::getDomain(my_keys_[0].first);
  saturation_liquid_key_ = plist.get<std::string>("saturation key", Keys::getKey(domain, "saturation_liquid"));

  std::string name = plist.get<std::string>("phase name");
  if (name == "liquid")
    phase_ = MULTIPHASE_PHASE_LIQUID;
  else if (name == "gas")
    phase_ = MULTIPHASE_PHASE_GAS;

  dependencies_.push_back(std::make_pair(saturation_liquid_key_, Tags::DEFAULT));
}


RelPermEvaluator::RelPermEvaluator(const RelPermEvaluator& other)
    : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(other),
      wrm_(other.wrm_) {};


/* ******************************************************************
* Copy constructor.
****************************************************************** */
Teuchos::RCP<Evaluator> RelPermEvaluator::Clone() const {
  return Teuchos::rcp(new RelPermEvaluator(*this));
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void RelPermEvaluator::Evaluate_(
    const State& S, const std::vector<CompositeVector*>& results)
{
  const auto& sat_c = *S.Get<CompositeVector>(saturation_liquid_key_).ViewComponent("cell");
  auto& result_c = *results[0]->ViewComponent("cell");

  int ncells = result_c.MyLength();
  for (int c = 0; c != ncells; ++c) {
    result_c[0][c] = wrm_->second[(*wrm_->first)[c]]->k_relative(sat_c[0][c], phase_);
  }
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void RelPermEvaluator::EvaluatePartialDerivative_(
    const State& S, const Key& wrt_key, const Tag& wrt_tag,
    const std::vector<CompositeVector*>& results)
{
  const auto& sat_c = *S.Get<CompositeVector>(saturation_liquid_key_).ViewComponent("cell");
  auto& result_c = *results[0]->ViewComponent("cell");

  int ncells = result_c.MyLength();
  for (int c = 0; c != ncells; ++c) {
    result_c[0][c] = wrm_->second[(*wrm_->first)[c]]->dKdS(sat_c[0][c], phase_);
  }
}

}  // namespace Flow
}  // namespace Amanzi
