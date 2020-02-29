/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  PCLiquidEvaluator is the interface between state/data and the model, an EOS.

  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_RELATIONS_PC_LIQUID_EVALUATOR_HH_
#define AMANZI_RELATIONS_PC_LIQUID_EVALUATOR_HH_

#include "secondary_variable_field_evaluator.hh"
#include "Factory.hh"

namespace Amanzi {
namespace Flow {

class PCLiqAtm;

class PCLiquidEvaluator : public SecondaryVariableFieldEvaluator {

 public:

  // constructor format for all derived classes
  explicit
  PCLiquidEvaluator(Teuchos::ParameterList& plist);
  PCLiquidEvaluator(const PCLiquidEvaluator& other);
  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);

  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

  virtual void EnsureCompatibility(const Teuchos::Ptr<State>& S) {
    S->RequireScalar("atmospheric_pressure");
    SecondaryVariableFieldEvaluator::EnsureCompatibility(S);
  }

  
  Teuchos::RCP<PCLiqAtm> get_PCLiqAtm() { return model_; }

 protected:
  // the actual model
  Teuchos::RCP<PCLiqAtm> model_;

  // Keys for fields
  // dependencies
  Key pres_key_;
  Key p_atm_key_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,PCLiquidEvaluator> factory_;

};

} // namespace
} // namespace

#endif
