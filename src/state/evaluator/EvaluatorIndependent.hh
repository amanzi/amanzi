/*
  State

  License: BSD
  Author: Ethan Coon

  An evaluator with no dependencies specified by a function.
*/

#ifndef AMANZI_INDEPENDENT_EVALUATOR_
#define AMANZI_INDEPENDENT_EVALUATOR_

#include "CompositeVectorFunction.hh"
#include "Evaluator.hh"
#include "Evaluator_Factory.hh"

namespace Amanzi {

class EvaluatorIndependent : public Evaluator {

 public:

  // ---------------------------------------------------------------------------
  // Constructors, assignement operators, etc
  // ---------------------------------------------------------------------------
  explicit
  EvaluatorIndependent(Teuchos::ParameterList& plist);
  EvaluatorIndependent(const EvaluatorIndependent& other) = default;

  // ---------------------------------------------------------------------------
  // Lazy evaluation of the evaluator.
  //
  // Updates the data, if needed.  Returns true if the value of the data has
  // changed since the last request for an update.
  // ---------------------------------------------------------------------------
  virtual bool Update(State& S, const Key& request) override;

  // ---------------------------------------------------------------------------
  // Lazy evaluation of derivatives of evaluator.
  //
  // Updates the derivative, if needed.  Returns true if the value of the
  // derivative with respect to wrt_key has changed since the last request for
  // an update.
  // ---------------------------------------------------------------------------
  virtual bool UpdateDerivative(State& S,
          const Key& request, const Key& wrt_key) override;

  virtual bool IsDependency(const State& S, const Key& key) const override;
  virtual bool ProvidesKey(const Key& key) const override;

  virtual void EnsureCompatibility(State& S) override;

  virtual std::string WriteToString() const override;

 protected:
  // ---------------------------------------------------------------------------
  // Does the actual work to update the value in the state.
  // ---------------------------------------------------------------------------
  virtual void Update_(State& S) = 0;

 protected:
  Key my_key_;
  Key my_tag_;

  double time_;
  bool temporally_variable_;
  bool computed_once_;

  KeySet requests_;
  Teuchos::ParameterList& plist_;
  VerboseObject vo_;

 private:
  static Utils::RegisteredFactory<Evaluator,EvaluatorIndependent> fac_;
};

} // namespace


#endif
