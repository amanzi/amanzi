/*
  EOS
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)

  Determining the molar fraction of a gas component within a gas mixture.
*/

#include "MolarFractionGasEvaluator.hh"

#include "EOSFactory.hh"

namespace Amanzi {
namespace AmanziEOS {

/* ******************************************************************
* Initialize various PK objects
****************************************************************** */
MolarFractionGasEvaluator::MolarFractionGasEvaluator(Teuchos::ParameterList& plist)
    : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(plist)
{
  // set up the actual model
  AMANZI_ASSERT(plist_.isSublist("vapor pressure model parameters"));
  EOSFactory<EOS_SaturatedVaporPressure> svp_fac;
  svp_model_ = svp_fac.CreateEOS(plist_.sublist("vapor pressure model parameters"));

  // process the list for my provided field.
  if (my_keys_.size() == 0)
    my_keys_.push_back(std::make_pair(plist_.get<std::string>("viscosity key", "viscosity_liquid"), Tags::DEFAULT));

  // set up dependencies
  std::string domain = Keys::getDomain(my_keys_[0].first);
  Tag tag = make_tag("");
  temp_key_ = plist_.get<std::string>("temperature key", Keys::getKey(domain, "temperature"));
  dependencies_.insert(std::make_pair(temp_key_, tag));
}


MolarFractionGasEvaluator::MolarFractionGasEvaluator(const MolarFractionGasEvaluator& other)
    : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(other),
      svp_model_(other.svp_model_),
      temp_key_(other.temp_key_) {};


Teuchos::RCP<Evaluator> MolarFractionGasEvaluator::Clone() const {
  return Teuchos::rcp(new MolarFractionGasEvaluator(*this));
}


void MolarFractionGasEvaluator::Evaluate_(
    const State& S, const std::vector<CompositeVector*>& results)
{
  // Pull dependencies out of state.
  auto temp = S.GetPtr<CompositeVector>(temp_key_);
  double p_atm = S.Get<double>("atmospheric_pressure");

  // evaluate p_s / p_atm
  for (auto comp = results[0]->begin(); comp != results[0]->end(); ++comp) {
    const Epetra_MultiVector& temp_v = *temp->ViewComponent(*comp);
    Epetra_MultiVector& result_v = *results[0]->ViewComponent(*comp);

    int count = results[0]->size(*comp);
    for (int i = 0; i != count; ++i) {
      result_v[0][i] = svp_model_->Pressure(temp_v[0][i]) / p_atm;
    }
  }
}


void MolarFractionGasEvaluator::EvaluatePartialDerivative_(
    const State& S, const Key& wrt_key, const Tag& wrt_tag,
    const std::vector<CompositeVector*>& results)
{
  AMANZI_ASSERT(wrt_key == temp_key_);

  auto temp = S.GetPtr<CompositeVector>(temp_key_);
  double p_atm = S.Get<double>("atmospheric_pressure");

  for (auto comp = results[0]->begin(); comp != results[0]->end(); ++comp) {
    const Epetra_MultiVector& temp_v = *temp->ViewComponent(*comp);
    Epetra_MultiVector& result_v = *results[0]->ViewComponent(*comp);

    int count = results[0]->size(*comp);
    for (int i = 0; i != count; ++i) {
      result_v[0][i] = svp_model_->DPressureDT(temp_v[0][i]) / p_atm;
    }
  }
}

}  // namespace AmanziEOS
}  // namespace Amanzi

