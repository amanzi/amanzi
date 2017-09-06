/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

A field evaluator base for secondary variables.  A Evaluator is a node in the
dependency graph.

Secondary variable evaluators, such as equations of state, water retention evaluators,
internal energy evaluators, etc should inherit this class, implementing the
missing Update_() and UpdateFieldDerivative_() methods.

This handles multiple fields that are related and/or better calculated as a
unit.  NOT to be confused wiht EvaluatorSecondary (note plural!)

------------------------------------------------------------------------- */

#ifndef STATE_EVALUATOR_SECONDARIES_HH_
#define STATE_EVALUATOR_SECONDARIES_HH_

#include <string>
#include <vector>

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

#include "State.hh"
#include "Evaluator.hh"

namespace Amanzi {

class EvaluatorSecondaries : public Evaluator {

 public:
  explicit
  EvaluatorSecondaries(Teuchos::ParameterList& plist);

  EvaluatorSecondaries(const EvaluatorSecondaries& other) = default;

  virtual Evaluator& operator=(const Evaluator& other) override;
  EvaluatorSecondaries& operator=(const EvaluatorSecondaries& other);

  // ---------------------------------------------------------------------------
  // Answers the question, has this Field changed since it was last requested
  // for Field Key reqest.  Updates the field if needed.
  // ---------------------------------------------------------------------------
  virtual bool Update(State& S, const Key& request) override;

  // ---------------------------------------------------------------------------
  // Answers the question, Has This Field's derivative with respect to Key
  // wrt_key changed since it was last requested for Field Key reqest.
  // Updates the derivative if needed.
  // ---------------------------------------------------------------------------
  virtual bool UpdateDerivative(State& S,
          const Key& request, const Key& wrt_key) override;

  virtual bool IsDependency(const State& S, const Key& key) const override;
  virtual bool ProvidesKey(const Key& key) const override;

  virtual void EnsureCompatibility(State& S) override;

  virtual std::string WriteToString() const override;

 protected:
  // These do the actual work
  virtual void EvaluateField_(State& S,
          const std::vector<Teuchos::Ptr<CompositeVector> >& results) = 0;
  virtual void EvaluateFieldPartialDerivative_(State& S, const Key& wrt_key,
          const std::vector<Teuchos::Ptr<CompositeVector> > & results) = 0;

  virtual void Update_(State& S);
  virtual void UpdateFieldDerivative_(State& S, const Key& wrt_key);
  virtual void CheckDerivative_(State& S, const Key& wrt_key);

 protected:
  std::vector<Key> my_keys_;
  Key my_tag_;

  KeySet requests_;
  KeyPairSet deriv_requests_;
  KeySet dependencies_;
  bool check_derivative_;

  VerboseObject vo_;
  Teuchos::ParameterList plist_;
  
}; // class Evaluator

} // namespace Amanzi

#endif
