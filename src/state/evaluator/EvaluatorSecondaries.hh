/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

A field evaluator base for algebraic variables.  A Evaluator is a node in the
dependency graph.

Algebraic variable evaluators, such as equations of state, water retention evaluators,
internal energy evaluators, etc should inherit this class, implementing the
missing Update_() and UpdateFieldDerivative_() methods.

This handles multiple fields that are related and/or better calculated as a
unit.  NOT to be confused wiht EvaluatorAlgebraic (note plural!)

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
  
//template<typename Data_t, typename DataFactory_t=NullFactory>
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
  virtual bool Update(State& S, const Key& request) override final;

  // ---------------------------------------------------------------------------
  // Answers the question, Has This Field's derivative with respect to Key
  // wrt_key changed since it was last requested for Field Key reqest.
  // Updates the derivative if needed.
  // ---------------------------------------------------------------------------
  virtual bool UpdateDerivative(State& S,
          const Key& request, const Key& wrt_key, const Key& wrt_tag) override final;

  virtual bool IsDependency(const State& S, const Key& key, const Key& tag) const override final;
  virtual bool ProvidesKey(const Key& key, const Key& tag) const override final;
  virtual bool IsDifferentiableWRT(const State& S, const Key& wrt_key, const Key& wrt_tag) const override {
    return IsDependency(S, wrt_key, wrt_tag);
  }
  
  virtual void EnsureCompatibility(State& S) override;

  virtual std::string WriteToString() const override;

 protected:
  // These do the actual work
  virtual void Evaluate_(const State& S,
                         const std::vector<Teuchos::Ptr<CompositeVector> > & results) = 0;
  virtual void EvaluatePartialDerivative_(const State& S,
          const Key& wrt_key, const Key& wrt_tag,
          const std::vector<Teuchos::Ptr<CompositeVector> > & results) = 0;
  
  virtual void Update_(State& S);
  virtual void UpdateDerivative_(State& S, const Key& wrt_key, const Key& wrt_tag);
  //virtual void CheckDerivative_(State& S, const Key& wrt_key, const Key& wrt_tag);
  
 protected:
  KeyPairVector my_keys_;

  KeySet requests_;
  KeyTripleSet deriv_requests_;
  KeyPairVector dependencies_;

  VerboseObject vo_;
  Teuchos::ParameterList plist_;
  
}; // class Evaluator

} // namespace Amanzi

#endif
