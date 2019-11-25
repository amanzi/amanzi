/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  EffectiveConcentrationEvaluator evaluates p_eff = max(p_atm, p_liquid), which is used for EOS.

  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "Factory.hh"
#include "effective_concentration_evaluator.hh"

namespace Amanzi {
namespace Relations {

EffectiveConcentrationEvaluator::EffectiveConcentrationEvaluator(Teuchos::ParameterList& plist) :
  SecondaryVariableFieldEvaluator(plist), index_id_(0) {
  if (my_key_ == std::string("")) {
    my_key_ = ep_plist_.get<std::string>("effective concentration key", "effective_concentration");
  }

  Key domain_name = Keys::getDomain(my_key_);

  // -- concentration
  conc_key_ = plist_.get<std::string>("concentration key",
          Keys::getKey(domain_name, "total_component_concentration"));
  dependencies_.insert(conc_key_);

  // -- logging
  if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
    Teuchos::OSTab tab = vo_->getOSTab();
    for (KeySet::const_iterator dep=dependencies_.begin();
         dep!=dependencies_.end(); ++dep) {
      *vo_->os() << " dep: " << *dep << std::endl;
    }
  }

}


EffectiveConcentrationEvaluator::EffectiveConcentrationEvaluator(
        const EffectiveConcentrationEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    conc_key_(other.conc_key_) {}


Teuchos::RCP<FieldEvaluator> EffectiveConcentrationEvaluator::Clone() const {
  return Teuchos::rcp(new EffectiveConcentrationEvaluator(*this));
}


void EffectiveConcentrationEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
                         const Teuchos::Ptr<CompositeVector>& result) {
  // Pull dependencies out of state.
  Teuchos::RCP<const CompositeVector> conc = S->GetFieldData(conc_key_);

  for (CompositeVector::name_iterator comp=result->begin();
       comp!=result->end(); ++comp) {
    const Epetra_MultiVector& conc_v = *(conc->ViewComponent(*comp,false));
    Epetra_MultiVector& result_v = *(result->ViewComponent(*comp,false));

    int count = result->size(*comp);
    for (int id=0; id!=count; ++id) {
      result_v[0][id] = conc_v[index_id_][id];
    }
  }
}

void EffectiveConcentrationEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {
  // Pull dependencies out of state.
  Teuchos::RCP<const CompositeVector> conc = S->GetFieldData(conc_key_);


  AMANZI_ASSERT(wrt_key == conc_key_);
  for (CompositeVector::name_iterator comp=result->begin();
       comp!=result->end(); ++comp) {
    const Epetra_MultiVector& conc_v = *(conc->ViewComponent(*comp,false));
    Epetra_MultiVector& result_v = *(result->ViewComponent(*comp,false));
    result_v.PutScalar(1.0);
  }
};

} // namespace
} // namespace


