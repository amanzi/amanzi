/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*
  State

*/

/*!

A non-leaf node in the tree, secondary variables have dependencies on other
evaluators.  They may provide one or more than one variable, and may depend on
a collection of things with arbitrary names, tags, and data types.

*/

#ifndef STATE_EVALUATOR_SECONDARY_HH_
#define STATE_EVALUATOR_SECONDARY_HH_

#include <string>
#include <vector>

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

#include "Evaluator.hh"
#include "State.hh"

namespace Amanzi {

class EvaluatorSecondary : public Evaluator {
 public:
  explicit EvaluatorSecondary(Teuchos::ParameterList& plist);

  EvaluatorSecondary(const EvaluatorSecondary& other) = default;

  virtual Evaluator& operator=(const Evaluator& other) override;
  EvaluatorSecondary& operator=(const EvaluatorSecondary& other);

  // Requires evaluators for the full dependency graph.
  virtual void EnsureEvaluators(State& S) override;

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
  virtual bool
  UpdateDerivative(State& S, const Key& request, const Key& wrt_key, const Tag& wrt_tag) override;

  virtual bool IsDependency(const State& S, const Key& key, const Tag& tag) const override;
  virtual bool ProvidesKey(const Key& key, const Tag& tag) const override;
  virtual bool
  IsDifferentiableWRT(const State& S, const Key& wrt_key, const Tag& wrt_tag) const override
  {
    // note, provides key means the value is 1, and there may be times we have to use this value...
    return ProvidesKey(wrt_key, wrt_tag) || IsDependency(S, wrt_key, wrt_tag);
  }

  virtual std::string WriteToString() const override;

 protected:
  // These do the actual work
  virtual void Update_(State& S) = 0;
  virtual void UpdateDerivative_(State& S, const Key& wrt_key, const Tag& wrt_tag) = 0;

  // may need to check for dependencies on other ranks
  virtual Teuchos::Ptr<const Comm_type> get_comm_(const State& s) const { return Teuchos::null; }

  // some may override this to force checkpointing
  virtual void EnsureCompatibility_Flags_(State& S);

 protected:
  KeyTagVector my_keys_;
  KeyTagSet dependencies_;
  bool nonlocal_dependencies_;

  KeySet requests_;
  DerivativeTripleSet deriv_requests_;

  bool updated_once_;

  Teuchos::ParameterList plist_;
  VerboseObject vo_;

}; // class Evaluator

} // namespace Amanzi

#endif
