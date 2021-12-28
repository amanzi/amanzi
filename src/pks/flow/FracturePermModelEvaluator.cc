/*
  Flow PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
           Konstantin Lipnikov (lipnikov@lanl.gov)

  The fracture permeability model evaluator simply calls the 
  permeability model with the correct arguments.
*/

#include "FlowDefs.hh"
#include "FracturePermModel_CubicLaw.hh"
#include "FracturePermModel_Linear.hh"
#include "FracturePermModelEvaluator.hh"
#include "FracturePermModelPartition.hh"

namespace Amanzi {
namespace Flow {

/* ******************************************************************
* Two constructors.
****************************************************************** */
FracturePermModelEvaluator::FracturePermModelEvaluator(
    Teuchos::ParameterList& plist, Teuchos::RCP<FracturePermModelPartition> fpm)
    : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(plist),
      fpm_(fpm)
{
  if (my_keys_.size() == 0) {
    my_keys_.push_back(std::make_pair(plist.get<std::string>("permeability key"), Tags::DEFAULT));
  }

  aperture_key_ = plist.get<std::string>("aperture key");
  dependencies_.push_back(std::make_pair(aperture_key_, Tags::DEFAULT));
}


FracturePermModelEvaluator::FracturePermModelEvaluator(const FracturePermModelEvaluator& other)
    : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(other),
      aperture_key_(other.aperture_key_),
      fpm_(other.fpm_) {};


/* ******************************************************************
* Copy constructor.
****************************************************************** */
Teuchos::RCP<Evaluator> FracturePermModelEvaluator::Clone() const {
  return Teuchos::rcp(new FracturePermModelEvaluator(*this));
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void FracturePermModelEvaluator::Evaluate_(
    const State& S, const std::vector<CompositeVector*>& results)
{
  auto& perm_c = *results[0]->ViewComponent("cell");
  const auto& aperture_c = *S.Get<CompositeVector>(aperture_key_).ViewComponent("cell");

  int ncells = perm_c.MyLength();
  for (int c = 0; c != ncells; ++c) {
    perm_c[0][c] = fpm_->second[(*fpm_->first)[c]]->Permeability(aperture_c[0][c]);
  }
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void FracturePermModelEvaluator::EvaluatePartialDerivative_(
    const State& S, const Key& wrt_key, const Tag& wrt_tag,
    const std::vector<CompositeVector*>& results)
{
  results[0]->PutScalar(0.0);
}

}  // namespace Flow
}  // namespace Amanzi
