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

#include "VolumetricStrainEvaluator.hh"

namespace Amanzi {
namespace Mechanics {

/* ******************************************************************
* Constructor
****************************************************************** */
VolumetricStrainEvaluator::VolumetricStrainEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(plist)
{
  displacement_key_ = plist.get<std::string>("displacement key", "displacement");
  dependencies_.insert(std::make_pair(displacement_key_, Tags::DEFAULT));
}


/* ******************************************************************
* A copy constructor.
****************************************************************** */
VolumetricStrainEvaluator::VolumetricStrainEvaluator(const VolumetricStrainEvaluator& other)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(other),
    displacement_key_(other.displacement_key_),
    op_(other.op_)
{}


/* ******************************************************************
* Clone with unclear yet purpose.
****************************************************************** */
Teuchos::RCP<Evaluator>
VolumetricStrainEvaluator::Clone() const
{
  return Teuchos::rcp(new VolumetricStrainEvaluator(*this));
}


/* ******************************************************************
* Evaluator
****************************************************************** */
void
VolumetricStrainEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& results)
{
  // operator may not be available during initialization time
  if (op_.get()) {
    auto u = S.Get<CompositeVector>(displacement_key_, Tags::DEFAULT);
    op_->ComputeVolumetricStrain(u, *results[0]);
  } else {
    results[0]->PutScalar(0.0);
  }
}

} // namespace Mechanics
} // namespace Amanzi
