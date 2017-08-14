/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

A field evaluator base for secondary variables.  A FieldEvaluator is a node in the
Phalanx-like dependency tree.

Secondary variable evaluators, such as equations of state, water retention evaluators,
internal energy evaluators, etc should inherit this class, implementing the
missing UpdateField_() and UpdateFieldDerivative_() methods.

------------------------------------------------------------------------- */

#ifndef STATE_SECONDARY_VARIABLE_FIELD_EVALUATOR_HH_
#define STATE_SECONDARY_VARIABLE_FIELD_EVALUATOR_HH_

#include <string>
#include <vector>

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

#include "State.hh"
#include "FieldEvaluator.hh"

namespace Amanzi {

class SecondaryVariableFieldEvaluator : public FieldEvaluator {

 public:
  explicit
  SecondaryVariableFieldEvaluator(Teuchos::ParameterList& plist);

  SecondaryVariableFieldEvaluator(const SecondaryVariableFieldEvaluator& other);

  // virtual destructor
  virtual ~SecondaryVariableFieldEvaluator() {}

  virtual void operator=(const FieldEvaluator& other);

  // ---------------------------------------------------------------------------
  // Answers the question, has this Field changed since it was last requested
  // for Field Key reqest.  Updates the field if needed.
  // ---------------------------------------------------------------------------
  virtual bool HasFieldChanged(const Teuchos::Ptr<State>& S, Key request);

  // ---------------------------------------------------------------------------
  // Answers the question, Has This Field's derivative with respect to Key
  // wrt_key changed since it was last requested for Field Key reqest.
  // Updates the derivative if needed.
  // ---------------------------------------------------------------------------
  virtual bool HasFieldDerivativeChanged(const Teuchos::Ptr<State>& S,
                                         Key request, Key wrt_key);

  virtual bool IsDependency(const Teuchos::Ptr<State>&S, Key key) const;
  virtual bool ProvidesKey(Key key) const;

  virtual void EnsureCompatibility(const Teuchos::Ptr<State>& S);

  virtual std::string WriteToString() const;

 protected:
  // These do the actual work
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
                             const Teuchos::Ptr<CompositeVector>& result) = 0;
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) = 0;

  virtual void UpdateField_(const Teuchos::Ptr<State>& S);
  virtual void UpdateFieldDerivative_(const Teuchos::Ptr<State>& S, Key wrt_key);

  virtual void CheckDerivative_(const Teuchos::Ptr<State>& S, Key wrt_key);

 protected:
  Key my_key_;
  KeySet requests_;
  KeyPairSet deriv_requests_;
  KeySet dependencies_;
  bool check_derivative_;
  bool nonlocal_dependencies_;

}; // class FieldEvaluator

} // namespace Amanzi

#endif
