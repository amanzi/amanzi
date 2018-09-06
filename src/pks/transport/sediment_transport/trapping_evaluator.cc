/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Determining the molar fraction of a gas component within a gas mixture.

  License: BSD
  Authors:
*/

#include "trapping_evaluator.hh"

namespace Amanzi {

TrappingRateEvaluator :: TrappingRateEvaluator(Teuchos::ParameterList& plist) :
  SecondaryVariableFieldEvaluator(plist) {

  Key domain_name = "surface";
  
  velocity_key_ = plist_.get<std::string>("velocity key",
                                     Keys::getKey(domain_name,"velocity"));

  tau_e_ = plist_.get<double>("critical shear stress");
  Qe_0_ = plist_.get<double>("empirical coefficient");
  gamma_ = plist_.get<double>("specific weight of water");
  lambda_ = plist_.get<double>("bottom friction coefficient");
    
  dependencies_.insert("surface-effective_pressure");
    
}

  
TrappingRateEvaluator ::TrappingRateEvaluator (const TrappingRateEvaluator & other) :
  SecondaryVariableFieldEvaluator(other),
  velocity_key_(other.velocity_key_) {

  tau_e_ = other.tau_e_;
  Qe_0_ = other.Qe_0_;
  gamma_ = other.gamma_;
  lambda_ = other.lambda_;
} 


Teuchos::RCP<FieldEvaluator> TrappingRateEvaluator ::Clone() const {
  return Teuchos::rcp(new TrappingRateEvaluator (*this));
}


void TrappingRateEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result) {

  result -> PutScalar(0.);

}

void TrappingRateEvaluator::EvaluateFieldPartialDerivative_ (const Teuchos::Ptr<State>& S,
                                                            Key wrt_key,
                                                            const Teuchos::Ptr<CompositeVector>& result) {
   AMANZI_ASSERT(0); 
}
  
  
} // namespace
