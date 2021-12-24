/*
  State

  Copyright 2010-201x held jointly, see COPYRIGHT.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon
*/

//! EvaluatorSecondaryMonotype is algebraic in data from the same time tag.

/*!

Algebraic evaluators are secondary evaluators that read only dependencies of
the same type as they calculate.  This allows requirements placed on the
calculated variable to be pushed down to the dependencies, checking
consistency, and also allows derivatives to be calculated automatically.

Algebraic variable evaluators, such as equations of state, water retention
evaluators, internal energy evaluators, etc should inherit this class,
implementing the missing Update_() and UpdateFieldDerivative_() methods.

*/

#ifndef AMANZI_STATE_EVALUATOR_ALGEBRAIC_HH_
#define AMANZI_STATE_EVALUATOR_ALGEBRAIC_HH_

#include <string>
#include <vector>

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

#include "exceptions.hh"

#include "EvaluatorSecondary.hh"

namespace Amanzi {

// By default, this class adds nothing on top of EvaluatorSecondary.
// Specializations can do useful things though.
template <typename Data_t, typename DataFactory_t = NullFactory>
class EvaluatorSecondaryMonotype : public EvaluatorSecondary {
public:
  using EvaluatorSecondary::EvaluatorSecondary;

  virtual void EnsureCompatibility(State& S) override;

protected:
  virtual void Update_(State& S) override;

  virtual void UpdateDerivative_(State& S, const Key& wrt_key,
          const Key& wrt_tag) override;
  
  virtual void EvaluatePartialDerivative_(const State& S,
          const Key& wrt_key, const Key& wrt_tag, const std::vector<Data_t*>& results) = 0;

  virtual void Evaluate_(const State& S, const std::vector<Data_t*>& results) = 0;
};


// implement generic versions
template <typename Data_t, typename DataFactory_t>
inline void
EvaluatorSecondaryMonotype<Data_t,DataFactory_t>::EnsureCompatibility(State& S) {
  // claim ownership, declare type
  for (auto keytag : my_keys_) {
    const auto& fac = S.Require<Data_t,DataFactory_t>(keytag.first, keytag.second, keytag.first);
  }

  // grab a factory for this
  auto akeytag = my_keys_[0];
  const auto& fac = S.Require<Data_t,DataFactory_t>(akeytag.first, akeytag.second);

  bool has_derivs = false;
  for (auto keytag : my_keys_)
    has_derivs |= S.HasDerivativeSet(keytag.first, keytag.second);
  if (has_derivs) {
    for (const auto& keytag : my_keys_) {
      for (const auto& deriv : S.GetDerivativeSet(keytag.first, keytag.second)) {
        auto wrt = Keys::splitKeyTag(deriv.first);
        auto& dfac = S.RequireDerivative<Data_t,DataFactory_t>(keytag.first, keytag.second,
                                                               wrt.first, wrt.second, keytag.first);
        dfac = fac; // derivatives are of the same type -- pointwise
      }
    }
  }

  // check plist for vis or checkpointing control
  EnsureCompatibility_Flags_(S);

  // requirements on dependencies as doubles
  for (auto& dep : dependencies_) {
    auto& depfac = S.Require<Data_t,DataFactory_t>(dep.first, dep.second);
    depfac = fac;
  }

  // require evaluators for dependencies and push down derivative info
  for (auto& dep : dependencies_) {
    auto& eval = S.RequireEvaluator(dep.first, dep.second);
    if (has_derivs) {
      for (const auto& deriv : S.GetDerivativeSet(akeytag.first, akeytag.second)) {
        auto wrt = Keys::splitKeyTag(deriv.first);
        if (eval.IsDifferentiableWRT(S, wrt.first, wrt.second)) {
          auto& dfac = S.RequireDerivative<Data_t,DataFactory_t>(dep.first, dep.second, wrt.first, wrt.second);
          dfac = fac;
        }
      }
    }
    S.GetEvaluator(dep.first, dep.second).EnsureCompatibility(S);
  }
}


// ---------------------------------------------------------------------------
// Updates the field value in state S.
// ---------------------------------------------------------------------------
template <typename Data_t, typename DataFactory_t>
inline void
EvaluatorSecondaryMonotype<Data_t,DataFactory_t>::Update_(State& S) {
  // vector of pointers to results
  std::vector<Data_t*> results;
  for (const auto& keytag : my_keys_) {
    results.push_back(&S.GetW<Data_t>(keytag.first, keytag.second, keytag.first));
  }

  // call the evaluate method
  Evaluate_(S, results);
}


// ---------------------------------------------------------------------------
// Updates the derivative
// ---------------------------------------------------------------------------
template <typename Data_t, typename DataFactory_t>
inline void
EvaluatorSecondaryMonotype<Data_t,DataFactory_t>::UpdateDerivative_(
    State& S, const Key& wrt_key, const Key& wrt_tag) {
  Errors::Message message;
  message << "EvaluatorSecondaryMonotype: "
          << my_keys_[0].first << "," << my_keys_[0].second
          << " has no implemented UpdateDerivative_()";
  throw(message);
}


// declare specializations
template <>
void EvaluatorSecondaryMonotype<double>::UpdateDerivative_(
    State& S, const Key& wrt_key, const Key& wrt_tag);

template <>
void EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>::UpdateDerivative_(
    State& S, const Key& wrt_key, const Key& wrt_tag);

template <>
void EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>::EnsureCompatibility(State& S);

}  // namespace Amanzi

#endif
