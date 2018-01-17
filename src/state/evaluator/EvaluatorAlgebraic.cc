/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

A field evaluator base for secondary variables.  A Evaluator is a node in the
Phalanx-like dependency tree.

Secondary variable evaluators, such as equations of state, water retention
evaluators, internal energy evaluators, etc should inherit this class,
implementing the missing Update_() and UpdateFieldDerivative_() methods.

Secondary secondary evaluator where all dependencies and this are
doubles.

------------------------------------------------------------------------- */

#include "EvaluatorAlgebraic.hh"

namespace Amanzi {

// ---------------------------------------------------------------------------
// Ensures that dependencies provide doubles
// ---------------------------------------------------------------------------
template <> void EvaluatorAlgebraic<double>::EnsureCompatibility(State &S) {
  // Ensure my field exists.  Requirements should be already set.
  ASSERT(my_key_ != std::string(""));

  // claim ownership, declare type
  S.Require<double>(my_key_, my_tag_, my_key_);

  // check plist for vis or checkpointing control
  bool io_my_key = plist_.get<bool>(std::string("visualize ") + my_key_, true);
  S.GetRecordW(my_key_, my_tag_, my_key_).set_io_vis(io_my_key);
  bool checkpoint_my_key =
      plist_.get<bool>(std::string("checkpoint ") + my_key_, false);
  S.GetRecordW(my_key_, my_tag_, my_key_).set_io_checkpoint(checkpoint_my_key);

  // requirements on dependencies as doubles
  for (auto &dep : dependencies_) {
    S.Require<double>(dep.first, dep.second);
  }

  // Recurse into the tree to propagate info to leaves.
  for (auto &dep : dependencies_) {
    S.RequireEvaluator(dep.first, dep.second)->EnsureCompatibility(S);
  }
}

// ---------------------------------------------------------------------------
// Updates the derivative for doubles
// ---------------------------------------------------------------------------
template <>
void EvaluatorAlgebraic<double>::UpdateDerivative_(State &S, const Key &wrt_key,
                                                   const Key &wrt_tag) {
  if (!S.HasDerivativeData(my_key_, my_tag_, wrt_key, wrt_tag)) {
    // Create the data structure
    S.RequireDerivative(my_key_, my_tag_, wrt_key, wrt_tag, my_key_);
  }

  double &dmy =
      S.GetDerivativeW<double>(my_key_, my_tag_, wrt_key, wrt_tag, my_key_);
  dmy = 0.;

  // dF/dx = sum_(deps) partial F/ partial dep * ddep/dx + partial F/partial x
  for (auto &dep : dependencies_) {

    if (wrt_key == dep.first && wrt_tag == dep.second) {
      // partial F / partial x
      double tmp;
      EvaluatePartialDerivative_(S, wrt_key, wrt_tag, tmp);
      dmy += tmp;

    } else if (S.GetEvaluator(dep.first, dep.second)
                   ->IsDependency(S, wrt_key, wrt_tag)) {
      // partial F / partial dep * ddep/dx
      // note this has already been Updated in the public version of this
      // function
      // -- ddep/dx
      const auto &ddep =
          S.GetDerivative<double>(dep.first, dep.second, wrt_key, wrt_tag);

      // -- partial F / partial dep
      double tmp;
      EvaluatePartialDerivative_(S, dep.first, dep.second, tmp);
      dmy += ddep * tmp;
    }
  }
}

// ---------------------------------------------------------------------------
// Ensures that dependencies provide the vector structure we need for this.
// ---------------------------------------------------------------------------
template <>
void EvaluatorAlgebraic<CompositeVector,
                        CompositeVectorSpace>::EnsureCompatibility(State &S) {
  // Ensure my field exists.  Requirements should be already set.
  ASSERT(my_key_ != std::string(""));

  // claim ownership
  CompositeVectorSpace &my_fac =
      S.Require<CompositeVector, CompositeVectorSpace>(my_key_, my_tag_,
                                                       my_key_);

  // check plist for vis or checkpointing control
  bool io_my_key = plist_.get<bool>(std::string("visualize ") + my_key_, true);
  S.GetRecordW(my_key_, my_tag_, my_key_).set_io_vis(io_my_key);
  bool checkpoint_my_key =
      plist_.get<bool>(std::string("checkpoint ") + my_key_, false);
  S.GetRecordW(my_key_, my_tag_, my_key_).set_io_checkpoint(checkpoint_my_key);

  // If my requirements have not yet been set, we'll have to hope they
  // get set by someone later.  For now just defer.
  std::string consistency_policy =
      plist_.get<std::string>("consistency policy", "give to child");
  if (consistency_policy == "none") {
    // no consistency required
    return;
  } else if (consistency_policy == "give to child") {
    // give my requirements to my children
    CompositeVectorSpace my_fac_copy(my_fac);
    my_fac_copy.SetOwned(false);

    // Loop over my dependencies, ensuring they meet the requirements.
    for (const auto &dep : dependencies_) {
      auto &fac = S.Require<CompositeVector, CompositeVectorSpace>(dep.first,
                                                                   dep.second);
      fac.Update(my_fac_copy);
    }

    // Recurse into the tree to propagate info to leaves.
    for (const auto &dep : dependencies_) {
      S.RequireEvaluator(dep.first, dep.second)->EnsureCompatibility(S);
    }
  } else if (consistency_policy.substr(0, 15) == "take from child") {
    // first call children to set their info
    for (const auto &dep : dependencies_) {
      S.RequireEvaluator(dep.first, dep.second)->EnsureCompatibility(S);
    }

    // take my requirements as the intersection of my children
    CompositeVectorSpace my_space;
    for (const auto &dep : dependencies_) {
      auto fac = S.Require<CompositeVector, CompositeVectorSpace>(dep.first,
                                                                  dep.second);
      if (my_space.size() == 0) {
        my_space = fac;
      } else {
        if (!my_space.SameAs(fac)) {
          if ((consistency_policy == "take from child: intersection" &&
               fac.SubsetOf(my_space)) ||
              (consistency_policy == "take from child: union" &&
               my_space.SubsetOf(fac))) {
            my_space = fac;
          }
        }
      }
    }
    my_fac.Update(my_space);
  }
}

// ---------------------------------------------------------------------------
// Updates the derivative for CompositeVectors
// ---------------------------------------------------------------------------
template <>
void EvaluatorAlgebraic<CompositeVector, CompositeVectorSpace>::
    UpdateDerivative_(State &S, const Key &wrt_key, const Key &wrt_tag) {
  if (!S.HasDerivativeData(my_key_, my_tag_, wrt_key, wrt_tag)) {
    // Create the data structure
    S.RequireDerivative(my_key_, my_tag_, wrt_key, wrt_tag, my_key_);
  }

  CompositeVector &dmy = S.GetDerivativeW<CompositeVector>(
      my_key_, my_tag_, wrt_key, wrt_tag, my_key_);
  dmy.PutScalarMasterAndGhosted(0.);

  // dF/dx = sum_(deps) partial F/ partial dep * ddep/dx + partial F/partial x
  for (auto &dep : dependencies_) {

    if (wrt_key == dep.first && wrt_tag == dep.second) {
      // partial F / partial x
      CompositeVector tmp(dmy);
      EvaluatePartialDerivative_(S, wrt_key, wrt_tag, tmp);
      dmy.Update(1., tmp, 1.);

    } else if (S.GetEvaluator(dep.first, dep.second)
                   ->IsDependency(S, wrt_key, wrt_tag)) {
      // partial F / partial dep * ddep/dx
      // note this has already been Updated in the public version of this
      // function
      // -- ddep/dx
      const auto &ddep = S.GetDerivative<CompositeVector>(dep.first, dep.second,
                                                          wrt_key, wrt_tag);

      // -- partial F / partial dep
      CompositeVector tmp(dmy);
      EvaluatePartialDerivative_(S, dep.first, dep.second, tmp);
      dmy.Multiply(1., ddep, tmp, 1.);
    }
  }
}

} // namespace Amanzi
