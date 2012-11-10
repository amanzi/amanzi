/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  EffectivePressureEvaluator evaluates p_eff = max(p_atm, p_liquid), which is used for EOS.

  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_RELATIONS_EOS_EVALUATOR_HH_
#define AMANZI_RELATIONS_EOS_EVALUATOR_HH_

#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace Relations {

class EffectivePressureEvaluator : public SecondaryVariableFieldEvaluator {

 public:

  // constructor format for all derived classes
  explicit
  EffectivePressureEvaluator(Teuchos::ParameterList& ep_plist);

  EffectivePressureEvaluator(const EffectivePressureEvaluator& other);
  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

 protected:

  // PList
  Teuchos::ParameterList ep_plist_;

  // Keys for fields
  // dependencies
  Key pres_key_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,EffectivePressureEvaluator> factory_;
};

} // namespace
} // namespace

#endif
