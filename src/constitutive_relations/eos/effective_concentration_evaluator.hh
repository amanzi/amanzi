/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  EffectiveConcentrationEvaluator evaluates p_eff = max(p_atm, p_liquid), which is used for EOS.

  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_EFFECTIVE_CONCENTRATION_EVALUATOR_HH_
#define AMANZI_EFFECTIVE_CONCENTRATION_EVALUATOR_HH_

#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace Relations {

class EffectiveConcentrationEvaluator : public SecondaryVariableFieldEvaluator {

 public:

  // constructor format for all derived classes
  explicit
  EffectiveConcentrationEvaluator(Teuchos::ParameterList& ep_plist);

  EffectiveConcentrationEvaluator(const EffectiveConcentrationEvaluator& other);
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
  Key conc_key_;
  int index_id_;
  

 private:
  static Utils::RegisteredFactory<FieldEvaluator,EffectiveConcentrationEvaluator> factory_;
};

} // namespace
} // namespace

#endif
