/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  EffectivePressureEvaluator evaluates p_eff = max(p_atm, p_liquid), which is used for EOS.

  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "factory.hh"
#include "effective_pressure_evaluator.hh"

namespace Amanzi {
namespace Relations {

// registry of method
Utils::RegisteredFactory<FieldEvaluator,EffectivePressureEvaluator> EffectivePressureEvaluator::factory_("effective_pressure");

EffectivePressureEvaluator::EffectivePressureEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist) {
  if (my_key_ == std::string("")) {
    my_key_ = ep_plist_.get<std::string>("effective pressure key", "effective_pressure");
  }

  setLinePrefix(my_key_+std::string(" evaluator"));

  pres_key_ = ep_plist_.get<std::string>("pressure key", "pressure");
  dependencies_.insert(pres_key_);
}


EffectivePressureEvaluator::EffectivePressureEvaluator(
        const EffectivePressureEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    pres_key_(other.pres_key_) {}


Teuchos::RCP<FieldEvaluator> EffectivePressureEvaluator::Clone() const {
  return Teuchos::rcp(new EffectivePressureEvaluator(*this));
}


void EffectivePressureEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
                         const Teuchos::Ptr<CompositeVector>& result) {
  // Pull dependencies out of state.
  Teuchos::RCP<const CompositeVector> pres = S->GetFieldData(pres_key_);
  const double& p_atm = *(S->GetScalarData("atmospheric_pressure"));

  // evaluate effective pressure as max(pres, p_atm)
  for (CompositeVector::name_iterator comp=result->begin();
       comp!=result->end(); ++comp) {
    const Epetra_MultiVector& pres_v = *(pres->ViewComponent(*comp,false));
    Epetra_MultiVector& result_v = *(result->ViewComponent(*comp,false));

    int count = result->size(*comp);
    for (int id=0; id!=count; ++id) {
      result_v[0][id] = std::max<double>(pres_v[0][id], p_atm);
    }
  }
}

void EffectivePressureEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {
  // Pull dependencies out of state.
  Teuchos::RCP<const CompositeVector> pres = S->GetFieldData(pres_key_);
  const double& p_atm = *(S->GetScalarData("atmospheric_pressure"));

  ASSERT(wrt_key == pres_key_);
  // pressure is max(pres, p_atm), so derivative is 1 or 0
  for (CompositeVector::name_iterator comp=result->begin();
       comp!=result->end(); ++comp) {
    const Epetra_MultiVector& pres_v = *(pres->ViewComponent(*comp,false));
    Epetra_MultiVector& result_v = *(result->ViewComponent(*comp,false));

    int count = result->size(*comp);
    for (int id=0; id!=count; ++id) {
      result_v[0][id] = pres_v[0][id] > p_atm ? 1.0 : 0.0;
    }
  }
};

} // namespace
} // namespace


