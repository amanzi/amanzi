/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Mechanics PK

  Evaluator for computing hydrostatic stress.
*/

#include "PDE_Elasticity.hh"

#include "HydrostaticStressEvaluator.hh"

namespace Amanzi {
namespace Mechanics {

/* ******************************************************************
* Constructor
****************************************************************** */
HydrostaticStressEvaluator::HydrostaticStressEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(plist)
{
  displacement_key_ = plist.get<std::string>("displacement key", "displacement");
  dependencies_.insert(std::make_pair(displacement_key_, Tags::DEFAULT));
}


/* ******************************************************************
* A copy constructor.
****************************************************************** */
HydrostaticStressEvaluator::HydrostaticStressEvaluator(const HydrostaticStressEvaluator& other)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(other),
    displacement_key_(other.displacement_key_),
    op_(other.op_)
{}


/* ******************************************************************
* Clone with unclear yet purpose.
****************************************************************** */
Teuchos::RCP<Evaluator>
HydrostaticStressEvaluator::Clone() const
{
  return Teuchos::rcp(new HydrostaticStressEvaluator(*this));
}


/* ******************************************************************
* Evaluator
****************************************************************** */
void
HydrostaticStressEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& results)
{
  // operator may not be available during initialization time
  if (op_.get()) {
    auto u = S.Get<CompositeVector>(displacement_key_, Tags::DEFAULT);
    op_->ComputeHydrostaticStress(u, *results[0]);
  } else {
    results[0]->PutScalar(0.0);
  }
}

} // namespace Mechanics
} // namespace Amanzi
