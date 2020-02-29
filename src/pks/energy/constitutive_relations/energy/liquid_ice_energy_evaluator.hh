/*
  The liquid+ice energy evaluator is an algebraic evaluator of a given model.

  Generated via evaluator_generator with:
Energy for a two-phase, liquid+ice evaluator.
    
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_ENERGY_LIQUID_ICE_ENERGY_EVALUATOR_HH_
#define AMANZI_ENERGY_LIQUID_ICE_ENERGY_EVALUATOR_HH_

#include "Factory.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace Energy {
namespace Relations {

class LiquidIceEnergyModel;

class LiquidIceEnergyEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  explicit
  LiquidIceEnergyEvaluator(Teuchos::ParameterList& plist);
  LiquidIceEnergyEvaluator(const LiquidIceEnergyEvaluator& other);

  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

  Teuchos::RCP<LiquidIceEnergyModel> get_model() { return model_; }

 protected:
  void InitializeFromPlist_();

  Key phi_key_;
  Key phi0_key_;
  Key sl_key_;
  Key nl_key_;
  Key ul_key_;
  Key si_key_;
  Key ni_key_;
  Key ui_key_;
  Key rho_r_key_;
  Key ur_key_;
  Key cv_key_;

  Teuchos::RCP<LiquidIceEnergyModel> model_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,LiquidIceEnergyEvaluator> reg_;

};

} //namespace
} //namespace
} //namespace

#endif