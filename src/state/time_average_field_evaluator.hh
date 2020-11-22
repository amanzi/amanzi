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

#ifndef STATE_TIME_AVERAGE_FIELD_EVALUATOR_HH_
#define STATE_TIME_AVERAGE_FIELD_EVALUATOR_HH_

#include <string>
#include <vector>

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

#include "State.hh"
#include "FieldEvaluator.hh"

namespace Amanzi {

class TimeAverageFieldEvaluator : public FieldEvaluator {

 public:
  explicit
  TimeAverageFieldEvaluator(Teuchos::ParameterList& plist);
  TimeAverageFieldEvaluator(const TimeAverageFieldEvaluator& other);

  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

  // virtual destructor
  virtual ~TimeAverageFieldEvaluator() {}

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
                                         Key request, Key wrt_key){return false;};

  virtual bool IsDependency(const Teuchos::Ptr<State>&S, Key key) const;
  virtual bool ProvidesKey(Key key) const;
  virtual void EnsureCompatibility(const Teuchos::Ptr<State>& S);
  virtual std::string WriteToString() const;

 protected:
  // These do the actual work
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
                             const Teuchos::Ptr<CompositeVector>& result);
  // virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
  //                                              Key wrt_key, const Teuchos::Ptr<CompositeVector>& result){};
  virtual void UpdateField_(const Teuchos::Ptr<State>& S);
  //  virtual void UpdateFieldDerivative_(const Teuchos::Ptr<State>& S, Key wrt_key){};


 protected:
  Key my_key_, base_key_;
  KeySet requests_;
  KeySet dependencies_;
  int number_bins_, current_bin_;
  double bin_period_, update_time_;
  std::vector<double> bin_time_;  
  Teuchos::RCP<CompositeVector> field_time_bins_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,TimeAverageFieldEvaluator> fac_;
  
}; // class FieldEvaluator

} // namespace Amanzi

#endif
