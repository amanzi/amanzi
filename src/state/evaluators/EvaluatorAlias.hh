/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*!

An evaluator that aliases to another key/tag pair, doing no work itself.

Note this may alias to any type of evaluator -- primary, secondary, or
independent.  These are typically not specified by the user in the input spec,
but instead manually constructed by a PK or MPC.

*/

#pragma once

#include "Evaluator.hh"

namespace Amanzi {

class EvaluatorAlias : public Evaluator {
 public:
  EvaluatorAlias(const Key& key, const Tag& alias_tag, const Tag& target_tag) :
    Evaluator(EvaluatorType::ALIAS),
    my_key_(key),
    my_tag_(alias_tag),
    target_key_(key),
    target_tag_(target_tag)
  {}

  virtual Teuchos::RCP<Evaluator> Clone() const override;

  EvaluatorAlias(const EvaluatorAlias& other) = default;
  EvaluatorAlias& operator=(const EvaluatorAlias& other);

  virtual Evaluator& operator=(const Evaluator& other) override;

  KeyTag get_target() const;

  bool Update(State& S, const Key& request) override;

  bool
  UpdateDerivative(State& S, const Key& requester, const Key& wrt_key, const Tag& wrt_tag) override;

  bool IsDependency(const State& S, const Key& key, const Tag& tag) const override;

  // Is this evaluator differentiable with respect to the primary variable in
  // wrt_key:wrt_tag?
  //
  // Searches the dependency graph to see if this evaluator depends upon the
  // evaluator named key.
  bool
  IsDifferentiableWRT(const State& S, const Key& wrt_key, const Tag& wrt_tag) const override;

  // Does this provide key?
  // Returns true if key is a field owned by this evaluator, false otherwise.
  virtual bool ProvidesKey(const Key& key, const Tag& tag) const override;

  // Requires evaluators for the full dependency graph.
  virtual void EnsureEvaluators(State& S) override;

  // Checks that all data requirements on dependencies of this evaluator are
  // satisfied by other evaluators in the dependency graph.
  virtual void EnsureCompatibility(State& S) override;

  // Virtual method for debugging, plotting the dependency graph, etc.
  virtual std::string WriteToString() const override;

  void SetChanged(State& S);

 protected:
  Key my_key_, target_key_;
  Tag my_tag_, target_tag_;
};

} // namespace Amanzi

