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

  virtual void EnsureCompatibility(State& S) override final;

 protected:
  virtual void Update_(State& S) override;

  virtual void UpdateDerivative_(State& S, const Key& wrt_key, const Tag& wrt_tag) override;

  virtual void EvaluatePartialDerivative_(const State& S,
          const Key& wrt_key, const Tag& wrt_tag, const std::vector<Data_t*>& results) = 0;

  virtual void Evaluate_(const State& S, const std::vector<Data_t*>& results) = 0;

 protected:
  // NOTE: EnsureCompatiblity() is often overridden by derived classes because
  // the factories are nearly, but not quite, identical between either various
  // keys evaluated by this evaluator, or from my key to dependency.  These
  // helper functions can be used by derived classes to do parts of, but not
  // all of, the default EnsureCompatibility().
  //
  // helper function that does the basics that likely all
  // EvaluatorSecondaryMonotype will use.
  virtual void EnsureCompatibility_Basics_(State& S);

  // helper function -- calls Require and claims ownership of all of my_keys_.
  // Called in Basics
  virtual void EnsureCompatibility_ClaimOwnership_(State& S);

  // helper function -- calls Require and claims ownership of all of my derivs.
  // Called in Basics
  virtual void EnsureCompatibility_ClaimDerivs_(State& S);

  // helper function -- calls Require on all dependencies
  // Called in Basics
  virtual void EnsureCompatibility_Deps_(State& S);

  // helper function -- requires deps to have derivatives
  // Called in Basics
  virtual void EnsureCompatibility_DepDerivs_(State& S);

  // helper function -- calls EnsureCompatibility on all dependency evaluators
  virtual void EnsureCompatibility_DepEnsureCompatibility_(State& S);

  // NOTE: the following helper functions are only implemented for
  // CompositeVectors.  They deal with the mesh and vector structure, and how
  // those relate between my_keys and dependencies.
  //
  // ensures a given structure of all of my_keys
  virtual void EnsureCompatibility_Structure_(State& S) {}

  // ensures all of my_keys have the same structure
  virtual void EnsureCompatibility_StructureSame_(State& S) {}

  // push metadata from my_keys to my_derivs
  virtual void EnsureCompatibility_DerivStructure_(State& S) {}

  // my_keys take their metadata from dependencies
  virtual void EnsureCompatibility_FromDeps_(State& S, const std::string& policy) {}

  // my_keys push requirements down to dependencies
  //
  // Calls the below method with fac = my_key's fac
  virtual void EnsureCompatibility_ToDeps_(State& S) {}

  // my_keys push requirements down to dependencies, using this factory, along
  // with a potential option for changing the dependency's domain name based on
  // the dependency's domain prefix.
  virtual void EnsureCompatibility_ToDeps_(State& S, const DataFactory_t& fac) {}

  // my_keys push requirements down to dependencies, each dependency has its
  // own mesh based on domain name but a common structure.
  virtual void EnsureCompatibility_ToDeps_(State& S,
          const std::vector<std::string>& names,
          const std::vector<AmanziMesh::Entity_kind>& locations,
          const std::vector<int>& num_dofs) {}

};


// implement generic versions
template <typename Data_t, typename DataFactory_t>
inline void
EvaluatorSecondaryMonotype<Data_t,DataFactory_t>::EnsureCompatibility(State& S)
{
  EnsureCompatibility_Basics_(S);
  EnsureCompatibility_Structure_(S);

  std::string consistency_policy =
      plist_.get<std::string>("consistency policy", "give to child");

  if (consistency_policy == "none") {
    // Requirements must be set by requiring code.
    //
    // Set requirements on my derivatives.
    EnsureCompatibility_DerivStructure_(S);

    // Ensure compatibility of children
    EnsureCompatibility_DepEnsureCompatibility_(S);

  } else if (consistency_policy == "give to child") {
    // Requirements set on my_keys must be provided by my dependencies.  This
    // is the most common, and also the most commonly customized with minor
    // tweaks to EnsureCompatibility_ToDeps_ implemented by specific
    // evaluators.

    // Set requirements on my derivatives.
    EnsureCompatibility_DerivStructure_(S);

    // Give my requirements to my children
    EnsureCompatibility_ToDeps_(S);

    // Now that the dependency requirements are set, call EnsureCompatibility
    EnsureCompatibility_DepEnsureCompatibility_(S);

  } else if (Keys::starts_with(consistency_policy, "take from child")) {
    // my_keys structure is set by the structure of my dependencies.  This is
    // most commonly used for evaluators where the structure is based on an
    // unknown, runtime-provided discretization.
    //
    // First call the dependencies EnsureCompatibility, which will set the
    // structure of my dependencies.
    EnsureCompatibility_DepEnsureCompatibility_(S);

    // Then, take requirements from children
    EnsureCompatibility_FromDeps_(S, consistency_policy);

    // Finally, push that into derivatives as well
    EnsureCompatibility_DerivStructure_(S);
  } else {
    Errors::Message msg;
    msg << "EvaluatorSecondaryMonotype for " << my_keys_.front().first << "@"
        << my_keys_.front().second.get() << ": invalid consistency policy \""
        << consistency_policy << "\"";
    Exceptions::amanzi_throw(msg);
  }
}


// ---------------------------------------------------------------------------
// Helper function that does the basics of ensuring requirements of evaluators
// are met.  Likely this is called by ALL implementations of
// EnsureCompatibility, whether there is a factory or not.
// ---------------------------------------------------------------------------
template <typename Data_t, typename DataFactory_t>
inline void
EvaluatorSecondaryMonotype<Data_t,DataFactory_t>::EnsureCompatibility_Basics_(State& S)
{
  EnsureCompatibility_ClaimOwnership_(S);
  EnsureCompatibility_Flags_(S);
  EnsureCompatibility_Deps_(S);
  EnsureCompatibility_ClaimDerivs_(S);
  EnsureCompatibility_DepDerivs_(S);
}

// ---------------------------------------------------------------------------
// Helper function that claims ownership of all of my_keys and sets the type.
// ---------------------------------------------------------------------------
template <typename Data_t, typename DataFactory_t>
inline void
EvaluatorSecondaryMonotype<Data_t,DataFactory_t>::EnsureCompatibility_ClaimOwnership_(
  State& S)
{
  // claim ownership, declare type
  for (auto keytag : my_keys_) {
    S.Require<Data_t,DataFactory_t>(keytag.first, keytag.second, keytag.first);
  }
}


// ---------------------------------------------------------------------------
// Helper function that sets requirements on derivatives.
//
// This assumes that all derivatives have the same structure as their primary
// variable's structure.
// ---------------------------------------------------------------------------
template <typename Data_t, typename DataFactory_t>
inline void
EvaluatorSecondaryMonotype<Data_t,DataFactory_t>::EnsureCompatibility_ClaimDerivs_(
  State& S)
{
  // require existence of derivatives
  bool has_derivs = false;
  for (auto keytag : my_keys_)
    has_derivs |= S.HasDerivativeSet(keytag.first, keytag.second);

  if (has_derivs) {
    // a derivative on any of my_keys_ means we need it on all other of my_keys_
    for (const auto& keytag : my_keys_) {
      if (S.HasDerivativeSet(keytag.first, keytag.second)) {
        for (const auto& deriv : S.GetDerivativeSet(keytag.first, keytag.second)) {
          auto wrt = Keys::splitKeyTag(deriv.first.get());
          for (const auto& other_keytag : my_keys_) {
            S.RequireDerivative<Data_t, DataFactory_t>(
              other_keytag.first, other_keytag.second,
              wrt.first, wrt.second, other_keytag.first);
          }
        }
      }
    }
  }
}


// ---------------------------------------------------------------------------
// Helper function that simply requires the dependencies to exist
// ---------------------------------------------------------------------------
template <typename Data_t, typename DataFactory_t>
inline void
EvaluatorSecondaryMonotype<Data_t,DataFactory_t>::EnsureCompatibility_Deps_(
  State& S)
{
  for (auto& dep : dependencies_) {
    auto& depfac = S.Require<Data_t,DataFactory_t>(dep.first, dep.second);
  }
}


// ---------------------------------------------------------------------------
// Helper function that requires the derivatives of dependencies needed to
// perform the chain rule.
//
// Note this assumes that EnsureCompatibility_ClaimDerivs_() and
// EnsureCompatibility_Deps_() have already been called.
// ---------------------------------------------------------------------------
template <typename Data_t, typename DataFactory_t>
inline void
EvaluatorSecondaryMonotype<Data_t,DataFactory_t>::EnsureCompatibility_DepDerivs_(State& S)
{
  auto akeytag = my_keys_.front();
  if (S.HasDerivativeSet(akeytag.first, akeytag.second)) {
    for (auto& dep : dependencies_) {
      auto& dep_eval = S.GetEvaluator(dep.first, dep.second);
      for (const auto& deriv : S.GetDerivativeSet(akeytag.first, akeytag.second)) {
        auto wrt = Keys::splitKeyTag(deriv.first.get());
        if (wrt != dep && dep_eval.IsDifferentiableWRT(S, wrt.first, wrt.second)) {
          S.RequireDerivative<Data_t,DataFactory_t>(dep.first, dep.second,
                  wrt.first, wrt.second);
        }
      }
    }
  }
}


// ---------------------------------------------------------------------------
// Helper function that recurses, calling EnsureCompatibility of dependencies.
// ---------------------------------------------------------------------------
template <typename Data_t, typename DataFactory_t>
inline void
EvaluatorSecondaryMonotype<Data_t,DataFactory_t>::EnsureCompatibility_DepEnsureCompatibility_(
  State& S)
{
  for (auto& dep : dependencies_) {
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
// Updates the derivative -- no valid derivative implementation for a generic
// type.
// ---------------------------------------------------------------------------
template <typename Data_t, typename DataFactory_t>
inline void
EvaluatorSecondaryMonotype<Data_t,DataFactory_t>::UpdateDerivative_(
    State& S, const Key& wrt_key, const Tag& wrt_tag) {
  Errors::Message msg;
  msg << "EvaluatorSecondaryMonotype: "
      << my_keys_[0].first << "," << my_keys_[0].second.get()
      << " has no implemented UpdateDerivative_()";
  throw(msg);
}


//
// SPECIALIZATIONS for double and CompositeVector
// ---------------------------------------------------------------------------


// ---------------------------------------------------------------------------
// Updates the derivative for doubles, doing the chain rule
// ---------------------------------------------------------------------------
template <>
void EvaluatorSecondaryMonotype<double>::UpdateDerivative_(
    State& S, const Key& wrt_key, const Tag& wrt_tag);

// ---------------------------------------------------------------------------
// Updates the derivative for CompositeVector, doing the chain rule pointwise
// ---------------------------------------------------------------------------
template <>
void EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>::UpdateDerivative_(
    State& S, const Key& wrt_key, const Tag& wrt_tag);

// ---------------------------------------------------------------------------
// The default EnsureCompatibility for CompositeVectors deals with how metadata
// relates between my_keys and dependencies.
// ---------------------------------------------------------------------------
template <>
void
EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>::EnsureCompatibility_StructureSame_(
  State& S);

template <>
void
EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>::EnsureCompatibility_DerivStructure_(
  State& S);

template <>
void
EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>::EnsureCompatibility_FromDeps_(
  State& S, const std::string& policy);

template <>
void
EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>::EnsureCompatibility_ToDeps_(
  State& S);

template <>
void
EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>::EnsureCompatibility_ToDeps_(
  State& S, const CompositeVectorSpace& fac);

template <>
void
EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>::EnsureCompatibility_ToDeps_(
  State& S,
  const std::vector<std::string>& names,
  const std::vector<AmanziMesh::Entity_kind>& locations,
  const std::vector<int>& num_dofs);


using EvaluatorSecondaryMonotypeCV = EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>;

}  // namespace Amanzi

#endif
