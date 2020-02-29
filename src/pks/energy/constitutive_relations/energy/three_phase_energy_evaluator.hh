/*
  The three phase energy evaluator is an algebraic evaluator of a given model.

  Generated via evaluator_generator with:
Energy for a three-phase, gas+liquid+ice evaluator.
    
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_ENERGY_THREE_PHASE_ENERGY_EVALUATOR_HH_
#define AMANZI_ENERGY_THREE_PHASE_ENERGY_EVALUATOR_HH_

#include "Factory.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace Energy {
namespace Relations {

class ThreePhaseEnergyModel;

class ThreePhaseEnergyEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  explicit
  ThreePhaseEnergyEvaluator(Teuchos::ParameterList& plist);
  ThreePhaseEnergyEvaluator(const ThreePhaseEnergyEvaluator& other);

  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

  Teuchos::RCP<ThreePhaseEnergyModel> get_model() { return model_; }

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
  Key sg_key_;
  Key ng_key_;
  Key ug_key_;
  Key rho_r_key_;
  Key ur_key_;
  Key cv_key_;

  Teuchos::RCP<ThreePhaseEnergyModel> model_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,ThreePhaseEnergyEvaluator> reg_;

};

} //namespace
} //namespace
} //namespace

#endif