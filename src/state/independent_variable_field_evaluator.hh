/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

A field evaluator with no dependencies specified by a function.

------------------------------------------------------------------------- */

#ifndef AMANZI_INDEPENDENT_FIELD_EVALUATOR_
#define AMANZI_INDEPENDENT_FIELD_EVALUATOR_

#include "composite_vector_function.hh"
#include "FieldEvaluator.hh"
#include "FieldEvaluator_Factory.hh"

namespace Amanzi {

class IndependentVariableFieldEvaluator : public FieldEvaluator {

 public:

  // ---------------------------------------------------------------------------
  // Constructors
  // ---------------------------------------------------------------------------
  explicit
  IndependentVariableFieldEvaluator(Teuchos::ParameterList& plist);
  IndependentVariableFieldEvaluator(const IndependentVariableFieldEvaluator& other);

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

  virtual bool IsDependency(const Teuchos::Ptr<State>& S, Key key) const;
  virtual bool ProvidesKey(Key key) const;

  virtual void EnsureCompatibility(const Teuchos::Ptr<State>& S);

  virtual std::string WriteToString() const;

 protected:
  // ---------------------------------------------------------------------------
  // Update the value in the state.
  // ---------------------------------------------------------------------------
  virtual void UpdateField_(const Teuchos::Ptr<State>& S) = 0;

 protected:
  Key my_key_;
  Teuchos::RCP<Functions::CompositeVectorFunction> func_;

  double time_;
  bool temporally_variable_;
  bool computed_once_;

  KeySet requests_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,IndependentVariableFieldEvaluator> fac_;
};

} // namespace


#endif
