/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/*
  License: see $ATS_DIR/COPYRIGHT
  Author: Ahmad Jan (jana@ornl.gov)
*/

//! Evaluates snow melt 
//! Source: USDA - Natural Resources Conservation Service, 
//! National Engineering Handbook (NEH) part 630, Chapter 11
/*!

Requires the following dependencies:

* `"air temperature key`" ``[string]`` **DOMAIN-air_temperature**
* `"precipitation snow key`" ``[string]`` **DOMAIN-precipitation_snow**
*/

#ifndef AMANZI_FLOW_RELATIONS_SMR_EVALUATOR_HH_
#define AMANZI_FLOW_RELATIONS_SMR_EVALUATOR_HH_

#include "Factory.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {

class SnowMeltRateEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  explicit
  SnowMeltRateEvaluator(Teuchos::ParameterList& plist);
  SnowMeltRateEvaluator(const SnowMeltRateEvaluator& other) = default;

  virtual Teuchos::RCP<FieldEvaluator> Clone() const override {
    return Teuchos::rcp(new SnowMeltRateEvaluator(*this));
  }

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result) override;
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) override;

 protected:
  
  Key at_key_, snow_key_;
  double melt_rate_;
  double snow_transition_depth_;

  Key domain_, domain_surf_;
  bool compatibility_checked_;
  
 private:
  static Utils::RegisteredFactory<FieldEvaluator,SnowMeltRateEvaluator> reg_;

};

} //namespace
} //namespace
} //namespace
#endif
