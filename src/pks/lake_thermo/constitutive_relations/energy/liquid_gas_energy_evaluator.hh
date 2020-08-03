/*
  The liquid+gas energy evaluator is an algebraic evaluator of a given model.

  Generated via evaluator_generator with:
Energy for a two-phase, liquid+water vapor evaluator.
    
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_ENERGY_LIQUID_GAS_ENERGY_EVALUATOR_HH_
#define AMANZI_ENERGY_LIQUID_GAS_ENERGY_EVALUATOR_HH_

#include "Factory.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace Energy {
namespace Relations {

class LiquidGasEnergyModel;

class LiquidGasEnergyEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  explicit
  LiquidGasEnergyEvaluator(Teuchos::ParameterList& plist);
  LiquidGasEnergyEvaluator(const LiquidGasEnergyEvaluator& other);

  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

  Teuchos::RCP<LiquidGasEnergyModel> get_model() { return model_; }

 protected:
  void InitializeFromPlist_();

  Key phi_key_;
  Key phi0_key_;
  Key sl_key_;
  Key nl_key_;
  Key ul_key_;
  Key sg_key_;
  Key ng_key_;
  Key ug_key_;
  Key rho_r_key_;
  Key ur_key_;
  Key cv_key_;

  Teuchos::RCP<LiquidGasEnergyModel> model_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,LiquidGasEnergyEvaluator> reg_;

};

} //namespace
} //namespace
} //namespace

#endif