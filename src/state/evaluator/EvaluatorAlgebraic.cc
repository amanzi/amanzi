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
  bool has_derivs = S.HasDerivativeSet(my_key_, my_tag_);
  if (has_derivs) {
    for (const auto& deriv : S.GetDerivativeSet(my_key_, my_tag_)) {
      auto wrt = Keys::splitKeyTag(deriv.first);
      S.RequireDerivative<double>(my_key_, my_tag_, wrt.first, wrt.second, my_key_);
    }
  }

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

  // require evaluators for dependencies and push down derivative info
  for (auto &dep : dependencies_) {
    auto& eval = S.RequireEvaluator(dep.first, dep.second);

    if (has_derivs) {
      for (const auto& deriv : S.GetDerivativeSet(my_key_, my_tag_)) {
        auto wrt = Keys::splitKeyTag(deriv.first);
        if (eval.IsDifferentiableWRT(S, wrt.first, wrt.second)) {
          S.RequireDerivative<double>(dep.first, dep.second, wrt.first, wrt.second);
        }
      }
    }
    S.GetEvaluator(dep.first, dep.second).EnsureCompatibility(S);
  }
}

// ---------------------------------------------------------------------------
// Updates the derivative for doubles
// ---------------------------------------------------------------------------
template <>
void EvaluatorAlgebraic<double>::UpdateDerivative_(State &S, const Key &wrt_key,
                                                   const Key &wrt_tag) {
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
                   .IsDependency(S, wrt_key, wrt_tag)) {
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
      S.Require<CompositeVector, CompositeVectorSpace>(my_key_, my_tag_, my_key_);
  bool has_derivs = S.HasDerivativeSet(my_key_, my_tag_);
  if (has_derivs) {
    for (const auto& deriv : S.GetDerivativeSet(my_key_, my_tag_)) {
      auto wrt = Keys::splitKeyTag(deriv.first);
      S.RequireDerivative<CompositeVector,CompositeVectorSpace>(my_key_, my_tag_,
              wrt.first, wrt.second, my_key_);
    }
  }

  // check plist for vis or checkpointing control
  bool io_my_key = plist_.get<bool>(std::string("visualize ") + my_key_, true);
  S.GetRecordW(my_key_, my_tag_, my_key_).set_io_vis(io_my_key);
  bool checkpoint_my_key =
      plist_.get<bool>(std::string("checkpoint ") + my_key_, false);
  S.GetRecordW(my_key_, my_tag_, my_key_).set_io_checkpoint(checkpoint_my_key);

  // require evaluators for dependencies
  for (auto &dep : dependencies_) S.RequireEvaluator(dep.first, dep.second);

  // set requirements on myself, my derivatives, my dependencies, and their derivatives
  std::string consistency_policy =
      plist_.get<std::string>("consistency policy", "give to child");
  if (consistency_policy == "none") {
    // I must have requirements since I won't get them here.
    // Set requirements on my derivatives.
    if (has_derivs) {
      for (const auto& deriv : S.GetDerivativeSet(my_key_, my_tag_)) {
        auto wrt = Keys::splitKeyTag(deriv.first);
        S.RequireDerivative<CompositeVector,CompositeVectorSpace>(my_key_, my_tag_,
                wrt.first, wrt.second, my_key_).Update(my_fac);
      }
    }

  } else if (consistency_policy == "give to child" &&
             my_fac.Mesh().get()) {
    // set requirements on my derivatives
    if (has_derivs) {
      for (const auto& deriv : S.GetDerivativeSet(my_key_, my_tag_)) {
        auto wrt = Keys::splitKeyTag(deriv.first);
        S.RequireDerivative<CompositeVector,CompositeVectorSpace>(my_key_, my_tag_,
                wrt.first, wrt.second, my_key_).Update(my_fac);
      }
    }

    // give my requirements to my children
    CompositeVectorSpace my_fac_copy(my_fac);
    my_fac_copy.SetOwned(false);

    for (const auto &dep : dependencies_) {
      auto &fac = S.Require<CompositeVector, CompositeVectorSpace>(dep.first, dep.second);
      fac.Update(my_fac_copy);

      // Set requirements on derivatives too
      if (has_derivs) {
        for (const auto& deriv : S.GetDerivativeSet(my_key_, my_tag_)) {
          auto wrt = Keys::splitKeyTag(deriv.first);
          if (S.GetEvaluator(dep.first, dep.second).IsDifferentiableWRT(S, wrt.first, wrt.second)) {
            auto &fac_d = S.RequireDerivative<CompositeVector,CompositeVectorSpace>(dep.first, dep.second, wrt.first, wrt.second);
            fac_d.Update(my_fac_copy);
          }
        }
      }

      // call ensure compatibility, now that dep requirements are set
      S.GetEvaluator(dep.first, dep.second).EnsureCompatibility(S);
    }

  } else if (consistency_policy.substr(0, 15) == "take from child") {
    // require derivatives of my children
    if (has_derivs) {
      for (const auto& deriv : S.GetDerivativeSet(my_key_, my_tag_)) {
        auto wrt = Keys::splitKeyTag(deriv.first);
        for (const auto& dep : dependencies_) {
          S.RequireDerivative<CompositeVector,CompositeVectorSpace>(dep.first, dep.second,
                  wrt.first, wrt.second);
        }
      }
    }

    // then call children to set their info, which will also set up deriv info
    for (const auto &dep : dependencies_) {
      S.GetEvaluator(dep.first, dep.second).EnsureCompatibility(S);
    }

    // take my requirements as the intersection or union of my children
    CompositeVectorSpace my_space;
    for (const auto &dep : dependencies_) {
      const auto& fac = S.Require<CompositeVector, CompositeVectorSpace>(dep.first, dep.second);
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

    // now push that into my derivative as well
    if (has_derivs) {
      for (const auto& deriv : S.GetDerivativeSet(my_key_, my_tag_)) {
        auto wrt = Keys::splitKeyTag(deriv.first);
        S.RequireDerivative<CompositeVector,CompositeVectorSpace>(my_key_, my_tag_,
                wrt.first, wrt.second, my_key_).Update(my_fac);
      }
    }

  }
}

// ---------------------------------------------------------------------------
// Updates the derivative for CompositeVectors
// ---------------------------------------------------------------------------
template <>
void EvaluatorAlgebraic<CompositeVector, CompositeVectorSpace>::
    UpdateDerivative_(State &S, const Key &wrt_key, const Key &wrt_tag) {
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

    } else if (S.GetEvaluator(dep.first, dep.second).IsDependency(S, wrt_key, wrt_tag)) {
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
