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
*/

// Amanzi
#include "PDE_HelperDiscretization.hh"

// Amanzi::Flow
#include "FlowDefs.hh"
#include "ApertureDarcyEvaluator.hh"

namespace Amanzi {
namespace Flow {

/* ******************************************************************
* Two constructors.
****************************************************************** */
ApertureDarcyEvaluator::ApertureDarcyEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(plist)
{
  if (my_keys_.size() == 0) {
    my_keys_.push_back(std::make_pair(plist_.get<std::string>("aperture key"), Tags::DEFAULT));
  }

  // dependencies
  ref_pressure_key_ = plist_.get<std::string>("reference pressure key");
  ref_aperture_key_ = plist_.get<std::string>("reference aperture key");

  pressure_key_ = plist_.get<std::string>("pressure key");
  compliance_key_ = plist_.get<std::string>("compliance key");

  dependencies_.insert(std::make_pair(pressure_key_, Tags::DEFAULT));
  dependencies_.insert(std::make_pair(ref_aperture_key_, Tags::DEFAULT));
}


ApertureDarcyEvaluator::ApertureDarcyEvaluator(const ApertureDarcyEvaluator& other)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(other),
    ref_pressure_key_(other.ref_pressure_key_),
    pressure_key_(other.pressure_key_),
    compliance_key_(other.compliance_key_){};


/* ******************************************************************
* Copy constructor.
****************************************************************** */
Teuchos::RCP<Evaluator>
ApertureDarcyEvaluator::Clone() const
{
  return Teuchos::rcp(new ApertureDarcyEvaluator(*this));
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void
ApertureDarcyEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& results)
{
  const auto& compliance_c = *S.Get<CompositeVector>(compliance_key_).ViewComponent("cell");
  const auto& a0_c = *S.Get<CompositeVector>(ref_aperture_key_).ViewComponent("cell");
  const auto& p0_c = *S.Get<CompositeVector>(ref_pressure_key_).ViewComponent("cell");
  const auto& p_c = *S.Get<CompositeVector>(pressure_key_).ViewComponent("cell");

  auto& a_c = *results[0]->ViewComponent("cell");
  int ncells = a_c.MyLength();

  for (int c = 0; c != ncells; ++c) {
    a_c[0][c] = a0_c[0][c] + compliance_c[0][c] * (p_c[0][c] - p0_c[0][c]);
  }
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void
ApertureDarcyEvaluator::EvaluatePartialDerivative_(const State& S,
                                                   const Key& wrt_key,
                                                   const Tag& wrt_tag,
                                                   const std::vector<CompositeVector*>& results)
{
  auto& result_c = *results[0]->ViewComponent("cell");
  const auto& compliance_c = *S.Get<CompositeVector>(compliance_key_).ViewComponent("cell");

  result_c = compliance_c;
}

} // namespace Flow
} // namespace Amanzi
