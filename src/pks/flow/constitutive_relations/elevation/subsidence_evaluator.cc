/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/*
  License: see $ATS_DIR/COPYRIGHT
  Authors: Ahmad Jan (jana@ornl.gov)
*/

//! Determine the change in subsidence (trough subsidence) from surface elevation.
/*!

* `"elevation key`" ``[string]`` **DOMAIN-elevation**
         The name of del_max, the max microtopography value.
*/

#include "subsidence_evaluator.hh"

namespace Amanzi {
namespace Flow {

SubsidenceEvaluator::SubsidenceEvaluator(Teuchos::ParameterList& plist) :
     SecondaryVariableFieldEvaluator(plist)
{
  Key domain = Keys::getDomain(my_key_);

  my_key_ = Keys::getKey(domain, "subsidence");
  // dependencies
  elev_key_ = Keys::readKey(plist_, domain, "elevation key", "elevation");
  dependencies_.insert(elev_key_);
  elev_init_key_ = Keys::readKey(plist_, domain, "elevation init key", "elevation_init");
  dependencies_.insert(elev_init_key_);
}

SubsidenceEvaluator::SubsidenceEvaluator(const SubsidenceEvaluator& other) :
  SecondaryVariableFieldEvaluator(other),
  elev_key_(other.elev_key_),
  elev_init_key_(other.elev_init_key_) {}

void
SubsidenceEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result)
{
  auto& res_c = *result->ViewComponent("cell",false);
  const auto& elev_init_c = *S->GetFieldData(elev_init_key_)->ViewComponent("cell",false);
  const auto& elev_c = *S->GetFieldData(elev_key_)->ViewComponent("cell",false);
  
  for (int c=0; c!=res_c.MyLength(); ++c){
      res_c[0][c] = elev_init_c[0][c] - elev_c[0][c];
  }
}


void
SubsidenceEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {}


} //namespace
} //namespace
