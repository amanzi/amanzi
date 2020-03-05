/*
  EOS
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)

  EffectivePressureEvaluator evaluates p_eff = max(p_atm, p_liquid).
*/

#include "EffectivePressureEvaluator.hh"

namespace Amanzi {
namespace AmanziEOS {

/* *******************************************************************
* Constructor.
******************************************************************* */
EffectivePressureEvaluator::EffectivePressureEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist) {
  if (my_key_ == std::string("")) {
    my_key_ = ep_plist_.get<std::string>("effective pressure key", "effective_pressure");
  }

  std::string domain("");
  auto end = my_key_.find_first_of("-");
  if (end != std::string::npos) domain = my_key_.substr(0, end);

  // -- pressure
  pres_key_ = plist_.get<std::string>("pressure key", Keys::getKey(domain, "pressure"));
  dependencies_.insert(pres_key_);

  // -- logging
  if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
    Teuchos::OSTab tab = vo_->getOSTab();
    for (auto dep = dependencies_.begin(); dep != dependencies_.end(); ++dep) {
      *vo_->os() << " dep: " << *dep << std::endl;
    }
  }
}


/* *******************************************************************
* Copy constructor.
******************************************************************* */
EffectivePressureEvaluator::EffectivePressureEvaluator(
        const EffectivePressureEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    pres_key_(other.pres_key_) {};


Teuchos::RCP<FieldEvaluator> EffectivePressureEvaluator::Clone() const {
  return Teuchos::rcp(new EffectivePressureEvaluator(*this));
}


/* *******************************************************************
* Main evaluator.
******************************************************************* */
void EffectivePressureEvaluator::EvaluateField_(
    const Teuchos::Ptr<State>& S, const Teuchos::Ptr<CompositeVector>& result)
{
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


/* *******************************************************************
* Main evaluator of derivatives.
******************************************************************* */
void EffectivePressureEvaluator::EvaluateFieldPartialDerivative_(
    const Teuchos::Ptr<State>& S,
    Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{
  // Pull dependencies out of state.
  Teuchos::RCP<const CompositeVector> pres = S->GetFieldData(pres_key_);
  const double& p_atm = *(S->GetScalarData("atmospheric_pressure"));

  AMANZI_ASSERT(wrt_key == pres_key_);
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
}

}  // namespace AmanziEOS
}  // namespace Amanzi


