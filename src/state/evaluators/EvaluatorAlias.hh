/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*

An evaluator that aliases to another key/tag pair, doing no work itself.

Note this may be any type of evaluator -- primary, secondary, or independent.
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

  virtual Teuchos::RCP<Evaluator> Clone() const override {
    return Teuchos::rcp(new EvaluatorAlias(*this));
  }

  EvaluatorAlias(const EvaluatorAlias& other) = default;
  EvaluatorAlias& operator=(const EvaluatorAlias& other) = default;

  Evaluator& operator=(const Evaluator& other) override
  {
    if (this != &other) {
      const EvaluatorAlias* other_p = dynamic_cast<const EvaluatorAlias*>(&other);
      AMANZI_ASSERT(other_p != NULL);
      *this = *other_p;
    }
    return *this;
  }

  KeyTag get_target() const {
    return KeyTag{target_key_, target_tag_};
  }

  bool Update(State& S, const Key& request) override {
    return S.GetEvaluator(target_key_, target_tag_).Update(S, request);
  }
  bool
  UpdateDerivative(State& S, const Key& requester, const Key& wrt_key, const Tag& wrt_tag) override {
    return S.GetEvaluator(target_key_, target_tag_).UpdateDerivative(S, requester, wrt_key, wrt_tag);
  }

  bool IsDependency(const State& S, const Key& key, const Tag& tag) const override {
    if (key == target_key_ && tag == target_tag_) return true;
    if (S.HasEvaluator(target_key_, target_tag_) &&
        S.GetEvaluator(target_key_, target_tag_).IsDependency(S, key, tag)) {
      return true;
    }
    return false;
  }

  // Is this evaluator differentiable with respect to the primary variable in
  // wrt_key:wrt_tag?
  //
  // Searches the dependency graph to see if this evaluator depends upon the
  // evaluator named key.
  bool
  IsDifferentiableWRT(const State& S, const Key& wrt_key, const Tag& wrt_tag) const override {
    return (wrt_key == my_key_ && wrt_tag == my_tag_) ||
      S.GetEvaluator(target_key_, target_tag_).IsDifferentiableWRT(S, wrt_key, wrt_tag);
  }

  // Does this provide key?
  // Returns true if key is a field owned by this evaluator, false otherwise.
  virtual bool ProvidesKey(const Key& key, const Tag& tag) const override {
    return (key == my_key_ && tag == my_tag_);
  }

  // Requires evaluators for the full dependency graph.
  virtual void EnsureEvaluators(State& S) override {
    S.RequireEvaluator(target_key_, target_tag_).EnsureEvaluators(S);
  }

  // Checks that all data requirements on dependencies of this evaluator are
  // satisfied by other evaluators in the dependency graph.
  virtual void EnsureCompatibility(State& S) override {
    S.GetEvaluator(target_key_, target_tag_).EnsureCompatibility(S);
    if (S.HasRecord(target_key_, target_tag_)) {
      S.GetRecordSetW(my_key_).AliasRecord(target_tag_, my_tag_);
    }
  }

  // Virtual method for debugging, plotting the dependency graph, etc.
  virtual std::string WriteToString() const override {
    std::string ret("Alias ");
    ret = ret + Keys::getKey(my_key_, my_tag_) + " --> " + Keys::getKey(target_key_, target_tag_);
    return ret;
  }

protected:
  Key my_key_, target_key_;
  Tag my_tag_, target_tag_;
};

} // namespace Amanzi

