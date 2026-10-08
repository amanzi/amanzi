/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Mechanics PK

  The evaluator calls the small strain models with the correct arguments.
*/

#include "SSMEvaluator.hh"
#include "SSMPartition.hh"

namespace Amanzi {
namespace Mechanics {

/* ******************************************************************
* Two constructors.
****************************************************************** */
SSMEvaluator::SSMEvaluator(Teuchos::ParameterList& plist, const Teuchos::RCP<SSMPartition>& ssm)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(plist), ssm_(ssm)
{
  if (my_keys_.size() == 0) {
    my_keys_.push_back(std::make_pair(
      plist_.get<std::string>("shear modulus key", "shear_modulus_key"), Tags::DEFAULT));
  }

  // my dependencies
  shear_strain_key_ = plist_.get<std::string>("shear strain key", "shear_strain");
  dependencies_.insert(std::make_pair(shear_strain_key_, Tags::DEFAULT));
}


SSMEvaluator::SSMEvaluator(const SSMEvaluator& other)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(other),
    ssm_(other.ssm_),
    shear_strain_key_(other.shear_strain_key_) {};


/* ******************************************************************
* Copy constructor.
****************************************************************** */
Teuchos::RCP<Evaluator>
SSMEvaluator::Clone() const
{
  return Teuchos::rcp(new SSMEvaluator(*this));
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void
SSMEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& results)
{
  const auto& gamma_c = *S.Get<CompositeVector>(shear_strain_key_).ViewComponent("cell");
  auto& result_c = *results[0]->ViewComponent("cell");
  int ncells = result_c.MyLength();

  for (int c = 0; c != ncells; ++c) {
    result_c[0][c] = ssm_->second[(*ssm_->first)[c]]->ShearStress(gamma_c[0][c]);
  }
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void
SSMEvaluator::EvaluatePartialDerivative_(const State& S,
                                         const Key& wrt_key,
                                         const Tag& wrt_tag,
                                         const std::vector<CompositeVector*>& results)
{}

} // namespace Mechanics
} // namespace Amanzi
