/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
      Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Flow PK

  The WRM Evaluator simply calls the WRM with the correct arguments.
*/

#include "FlowDefs.hh"
#include "WRMEvaluator.hh"
#include "WRMPartition.hh"

namespace Amanzi {
namespace Flow {

/* ******************************************************************
* Two constructors.
****************************************************************** */
WRMEvaluator::WRMEvaluator(Teuchos::ParameterList& plist, const Teuchos::RCP<WRMPartition>& wrm)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(plist), wrm_(wrm)
{
  InitializeFromPlist_();
}


WRMEvaluator::WRMEvaluator(const WRMEvaluator& other)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(other),
    wrm_(other.wrm_),
    pressure_key_(other.pressure_key_){};


/* ******************************************************************
* Copy constructor.
****************************************************************** */
Teuchos::RCP<Evaluator>
WRMEvaluator::Clone() const
{
  return Teuchos::rcp(new WRMEvaluator(*this));
}


/* ******************************************************************
* Initialization.
****************************************************************** */
void
WRMEvaluator::InitializeFromPlist_()
{
  if (my_keys_.size() == 0) {
    my_keys_.push_back(std::make_pair(
      plist_.get<std::string>("saturation key", "saturation_liquid"), Tags::DEFAULT));
  }

  // my dependency is pressure.
  pressure_key_ = plist_.get<std::string>("pressure key", "pressure");
  dependencies_.insert(std::make_pair(pressure_key_, Tags::DEFAULT));
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void
WRMEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& results)
{
  auto& sat_c = *results[0]->ViewComponent("cell");
  const auto& pres_c = *S.Get<CompositeVector>(pressure_key_).ViewComponent("cell");
  double patm = S.Get<double>("atmospheric_pressure");

  int ncells = sat_c.MyLength();
  for (int c = 0; c != ncells; ++c) {
    sat_c[0][c] = wrm_->second[(*wrm_->first)[c]]->saturation(patm - pres_c[0][c]);
  }
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void
WRMEvaluator::EvaluatePartialDerivative_(const State& S,
                                         const Key& wrt_key,
                                         const Tag& wrt_tag,
                                         const std::vector<CompositeVector*>& results)
{
  auto& sat_c = *results[0]->ViewComponent("cell");
  const auto& pres_c = *S.Get<CompositeVector>(pressure_key_).ViewComponent("cell");
  double patm = S.Get<double>("atmospheric_pressure");

  int ncells = sat_c.MyLength();
  for (int c = 0; c != ncells; ++c) {
    // Negative sign indicates that dSdP = -dSdPc.
    sat_c[0][c] = -wrm_->second[(*wrm_->first)[c]]->dSdPc(patm - pres_c[0][c]);
  }
}

} // namespace Flow
} // namespace Amanzi
