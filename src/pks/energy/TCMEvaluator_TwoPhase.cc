/*
  Energy

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)

  Interface for a thermal conductivity model with two phases.
*/

#include "dbc.hh"
#include "TCMFactory_TwoPhase.hh"
#include "TCMEvaluator_TwoPhase.hh"

namespace Amanzi {
namespace Energy {

/* ******************************************************************
* Constructor.
****************************************************************** */
TCMEvaluator_TwoPhase::TCMEvaluator_TwoPhase(Teuchos::ParameterList& plist)
    : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(plist)
{
  if (my_keys_.size() == 0) {
    my_keys_.push_back(make_pair(plist_.get<std::string>("thermal conductivity key"), Tags::DEFAULT));
  }
  auto prefix = Keys::getDomainPrefix(my_keys_[0].first);

  porosity_key_ = plist_.get<std::string>("porosity key", prefix + "porosity");
  dependencies_.insert(std::make_pair(porosity_key_, Tags::DEFAULT));

  saturation_key_ = plist_.get<std::string>("saturation key", prefix + "saturation_liquid");
  dependencies_.insert(std::make_pair(saturation_key_, Tags::DEFAULT));

  AMANZI_ASSERT(plist_.isSublist("thermal conductivity parameters"));
  Teuchos::ParameterList sublist = plist_.sublist("thermal conductivity parameters");
  TCMFactory_TwoPhase fac;
  tc_ = fac.CreateTCM(sublist);
}


/* ******************************************************************
* Copy constructor.
****************************************************************** */
TCMEvaluator_TwoPhase::TCMEvaluator_TwoPhase(const TCMEvaluator_TwoPhase& other)
    : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(other),
      tc_(other.tc_),
      porosity_key_(other.porosity_key_),
      saturation_key_(other.saturation_key_) {}


/* ******************************************************************
* TBW.
****************************************************************** */
Teuchos::RCP<Evaluator> TCMEvaluator_TwoPhase::Clone() const {
  return Teuchos::rcp(new TCMEvaluator_TwoPhase(*this));
}


/* ******************************************************************
* Evaluator body.
****************************************************************** */
void TCMEvaluator_TwoPhase::Evaluate_(
    const State& S, const std::vector<CompositeVector*>& results)
{
  // pull out the dependencies
  auto poro = S.GetPtr<CompositeVector>(porosity_key_);
  auto sat = S.GetPtr<CompositeVector>(saturation_key_);

  for (auto comp = results[0]->begin(); comp != results[0]->end(); ++comp) {
    // much more efficient to pull out vectors first
    const Epetra_MultiVector& poro_v = *poro->ViewComponent(*comp);
    const Epetra_MultiVector& sat_v = *sat->ViewComponent(*comp);
    Epetra_MultiVector& result_v = *results[0]->ViewComponent(*comp);

    int ncomp = results[0]->size(*comp);
    for (int i = 0; i != ncomp; ++i) {
      result_v[0][i] = tc_->ThermalConductivity(poro_v[0][i], sat_v[0][i]);
    }
  }
}


/* ******************************************************************
* Evaluator of derivarives.
****************************************************************** */
void TCMEvaluator_TwoPhase::EvaluatePartialDerivative_(
    const State& S, const Key& wrt_key, const Tag& wrt_tag,
    const std::vector<CompositeVector*>& results)
{
  AMANZI_ASSERT(0);
}

}  // namespace Energy
}  // namespace Amanzi
