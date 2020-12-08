/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
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

#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {


// ---------------------------------------------------------------------------
// Updates the derivative for doubles
// ---------------------------------------------------------------------------
template <>
void
EvaluatorSecondaryMonotype<double>::UpdateDerivative_(State& S,
                                                      const Key& wrt_key,
                                                      const Key& wrt_tag)
{
  // zero out derivatives
  std::vector<double*> results(my_keys_.size());
  int j = 0;
  for (const auto& keytag : my_keys_) {
    results[j] = &S.GetDerivativeW<double>(
      keytag.first, keytag.second, wrt_key, wrt_tag, keytag.first);
    *results[j] = 0.;
    ++j;
  }

  // if provides key, then the result is 1
  if (ProvidesKey(wrt_key, wrt_tag)) {
    int i = std::find(my_keys_.begin(), my_keys_.end(), std::make_pair(wrt_key, wrt_tag))
            - my_keys_.begin();
    AMANZI_ASSERT(i < my_keys_.size()); // ensured by IsDifferentiableWRT() check previously
    *results[i] = 1.0;
    return;
  }
  
  // dF/dx = sum_(deps) partial F/ partial dep * ddep/dx + partial F/partial x
  for (auto& dep : dependencies_) {
    if (wrt_key == dep.first && wrt_tag == dep.second) {
      // partial F / partial x
      std::vector<double> tmp_data(my_keys_.size(), 0.);
      std::vector<double*> tmp(my_keys_.size());
      for (int i = 0; i != my_keys_.size(); ++i) { tmp[i] = &tmp_data[i]; }
      EvaluatePartialDerivative_(S, wrt_key, wrt_tag, tmp);
      for (int i = 0; i != my_keys_.size(); ++i) (*results[i]) += tmp_data[i];


    } else if (S.GetEvaluator(dep.first, dep.second)
                 .IsDependency(S, wrt_key, wrt_tag)) {
      // partial F / partial dep * ddep/dx
      // note this has already been Updated in the public version of this
      // function

      // -- ddep/dx
      const auto& ddep =
        S.GetDerivative<double>(dep.first, dep.second, wrt_key, wrt_tag);

      // -- partial F / partial dep
      std::vector<double> tmp_data(my_keys_.size(), 0.);
      std::vector<double*> tmp(my_keys_.size());
      for (int i = 0; i != my_keys_.size(); ++i) { tmp[i] = &tmp_data[i]; }
      EvaluatePartialDerivative_(S, dep.first, dep.second, tmp);

      // sum
      for (int i = 0; i != my_keys_.size(); ++i)
        (*results[i]) += ddep * tmp_data[i];
    }
  }
}


// ---------------------------------------------------------------------------
// Updates the derivative for CompositeVectors
// ---------------------------------------------------------------------------
template <>
void
EvaluatorSecondaryMonotype<
  CompositeVector, CompositeVectorSpace>::UpdateDerivative_(State& S,
                                                            const Key& wrt_key,
                                                            const Key& wrt_tag)
{
  // zero out values
  std::vector<CompositeVector*> results(my_keys_.size());
  int j = 0;
  for (const auto& keytag : my_keys_) {
    results[j] = &S.GetDerivativeW<CompositeVector>(
      keytag.first, keytag.second, wrt_key, wrt_tag, keytag.first);
    results[j]->putScalarMasterAndGhosted(0.);
    ++j;
  }

  // if provides key, then the result is 1
  if (ProvidesKey(wrt_key, wrt_tag)) {
    Errors::Message msg;
    msg << "EvaluatorSecondary (" << my_keys_[0].first << "," << my_keys_[0].second << ") provides key (" << wrt_key << "," << wrt_tag << ") and so should not be differentiated with respect to this key.";
    throw(msg);
  }

  // dF/dx = sum_(deps) partial F/ partial dep * ddep/dx + partial F/partial x
  for (auto& dep : dependencies_) {
    auto& dep_eval = S.GetEvaluator(dep.first, dep.second);
    if (dep_eval.ProvidesKey(wrt_key, wrt_tag)) {
      // partial F / partial x
      std::vector<CompositeVector> tmp_data(my_keys_.size(), *results[0]);
      std::vector<CompositeVector*> tmp(my_keys_.size());
      for (int i = 0; i != my_keys_.size(); ++i) { tmp[i] = &tmp_data[i]; }
      EvaluatePartialDerivative_(S, wrt_key, wrt_tag, tmp);
      for (int i = 0; i != my_keys_.size(); ++i)
        results[i]->update(1., tmp_data[i], 1.);

    } else if (S.GetEvaluator(dep.first, dep.second)
                 .IsDifferentiableWRT(S, wrt_key, wrt_tag)) {
      // partial F / partial dep * ddep/dx
      // note this has already been Updated in the public version of this
      // function

      // -- ddep/dx
      const auto& ddep = S.GetDerivative<CompositeVector>(
        dep.first, dep.second, wrt_key, wrt_tag);

      // -- partial F / partial dep
      std::vector<CompositeVector> tmp_data(my_keys_.size(), *results[0]);
      std::vector<CompositeVector*> tmp(my_keys_.size());
      for (int i = 0; i != my_keys_.size(); ++i) { tmp[i] = &tmp_data[i]; }
      EvaluatePartialDerivative_(S, dep.first, dep.second, tmp);

      // sum
      for (int i = 0; i != my_keys_.size(); ++i)
        results[i]->elementWiseMultiply(1., ddep, tmp_data[i], 1.);
    }
  }
}


// ---------------------------------------------------------------------------
// Ensures that dependencies provide the vector structure we need for this.
// ---------------------------------------------------------------------------
template <>
void
EvaluatorSecondaryMonotype<CompositeVector,
                           CompositeVectorSpace>::EnsureCompatibility(State& S)
{
  // claim ownership, declare type
  for (auto keytag : my_keys_) {
    const auto& fac = S.Require<CompositeVector, CompositeVectorSpace>(
      keytag.first, keytag.second, keytag.first);
  }

  // grab my factory
  auto akeytag = my_keys_[0];
  auto& my_fac = S.Require<CompositeVector, CompositeVectorSpace>(
    akeytag.first, akeytag.second, akeytag.first);

  // require existence of derivatives
  bool has_derivs = false;
  for (auto keytag : my_keys_)
    has_derivs |= S.HasDerivativeSet(keytag.first, keytag.second);
  if (has_derivs) {
    for (const auto& keytag : my_keys_) {
      for (const auto& deriv :
           S.GetDerivativeSet(keytag.first, keytag.second)) {
        auto wrt = Keys::splitKeyTag(deriv.first);
        S.RequireDerivative<CompositeVector, CompositeVectorSpace>(
          keytag.first, keytag.second, wrt.first, wrt.second, keytag.first);
      }
    }
  }

  // check plist for vis or checkpointing control
  EnsureCompatibility_Flags_(S);
  if (my_fac.Mesh().get() && !db_.get()) {
    db_ = Teuchos::rcp(new Debugger(my_fac.Mesh(), my_keys_[0].first, plist_));
  }

  // require evaluators for dependencies
  // NOTE: unnecessary now, as require put in State::SetEvaluator()
  //  for (auto& dep : dependencies_) S.RequireEvaluator(dep.first, dep.second);

  // set requirements on myself, my derivatives, my dependencies, and their
  // derivatives
  std::string consistency_policy =
    plist_.get<std::string>("consistency policy", "give to child");
  if (consistency_policy == "none") {
    // I must have requirements since I won't get them here.
    // Set requirements on my derivatives.
    if (has_derivs) {
      for (const auto& keytag : my_keys_) {
        for (const auto& deriv :
             S.GetDerivativeSet(keytag.first, keytag.second)) {
          auto wrt = Keys::splitKeyTag(deriv.first);
          S.RequireDerivative<CompositeVector, CompositeVectorSpace>(
             keytag.first, keytag.second, wrt.first, wrt.second, keytag.first)
            .Update(my_fac);
        }
      }
    }

  } else if (consistency_policy == "give to child" && my_fac.Mesh().get()) {
    // set requirements on my derivatives
    if (has_derivs) {
      for (const auto& keytag : my_keys_) {
        for (const auto& deriv :
             S.GetDerivativeSet(keytag.first, keytag.second)) {
          auto wrt = Keys::splitKeyTag(deriv.first);
          S.RequireDerivative<CompositeVector, CompositeVectorSpace>(
             keytag.first, keytag.second, wrt.first, wrt.second, keytag.first)
            .Update(my_fac);
        }
      }
    }

    // give my requirements to my children
    for (const auto& dep : dependencies_) {
      auto& fac =
        S.Require<CompositeVector, CompositeVectorSpace>(dep.first, dep.second);
      fac.Update(my_fac);

      // Set requirements on derivatives too
      if (has_derivs) {
        for (const auto& keytag : my_keys_) {
          for (const auto& deriv :
               S.GetDerivativeSet(keytag.first, keytag.second)) {
            auto wrt = Keys::splitKeyTag(deriv.first);
            auto& dep_eval = S.GetEvaluator(dep.first, dep.second);
            if (dep_eval.IsDifferentiableWRT(S, wrt.first, wrt.second) &&
                !dep_eval.ProvidesKey(wrt.first, wrt.second)) {
              S.RequireDerivative<CompositeVector, CompositeVectorSpace>(
                 dep.first, dep.second, wrt.first, wrt.second)
                .Update(my_fac);
            }
          }
        }
      }

      // call ensure compatibility, now that dep requirements are set
      S.GetEvaluator(dep.first, dep.second).EnsureCompatibility(S);
    }

  } else if (consistency_policy.substr(0, 15) == "take from child") {
    // require derivatives of my children
    if (has_derivs) {
      for (const auto& keytag : my_keys_) {
        for (const auto& deriv :
             S.GetDerivativeSet(keytag.first, keytag.second)) {
          auto wrt = Keys::splitKeyTag(deriv.first);
          for (const auto& dep : dependencies_) {
            S.RequireDerivative<CompositeVector, CompositeVectorSpace>(
              dep.first, dep.second, wrt.first, wrt.second);
          }
        }
      }
    }

    // then call children to set their info, which will also set up deriv info
    for (const auto& dep : dependencies_) {
      S.GetEvaluator(dep.first, dep.second).EnsureCompatibility(S);
    }

    // take my requirements as the intersection or union of my children
    CompositeVectorSpace my_space;
    for (const auto& dep : dependencies_) {
      const auto& fac =
        S.Require<CompositeVector, CompositeVectorSpace>(dep.first, dep.second);
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

    // update my facs with this info
    my_fac.Update(my_space);
    for (auto keytag : my_keys_) {
      S.Require<CompositeVector, CompositeVectorSpace>(
         keytag.first, keytag.second, keytag.first)
        .Update(my_fac);
    }

    // now push that into my derivative as well
    if (has_derivs) {
      for (auto keytag : my_keys_) {
        for (const auto& deriv :
             S.GetDerivativeSet(keytag.first, keytag.second)) {
          auto wrt = Keys::splitKeyTag(deriv.first);
          S.RequireDerivative<CompositeVector, CompositeVectorSpace>(
             keytag.first, keytag.second, wrt.first, wrt.second, keytag.first)
            .Update(my_fac);
        }
      }
    }
  }
}



template <>
void
EvaluatorSecondaryMonotype<CompositeVector,
                           CompositeVectorSpace>::Debug_(const State& S)
{
  if (vo_.os_OK(Teuchos::VERB_HIGH)) {
    std::vector<Teuchos::Ptr<const CompositeVector>> my_vecs;
    std::vector<std::string> my_names;
    for (const auto& key_tag : dependencies_) {
      if (key_tag.first != "time") { // this helps make life a bit simpler...
        my_names.push_back(Keys::getKeyTag(key_tag.first, key_tag.second));
        my_vecs.push_back(S.GetPtr<CompositeVector>(key_tag.first, key_tag.second).ptr());
      }
    }
    db_->WriteVectors(my_names, my_vecs);
    db_->WriteDivider();

    my_vecs.clear(); my_names.clear();
    for (const auto& key_tag : my_keys_) {
      my_names.push_back(Keys::getKeyTag(key_tag.first, key_tag.second));
      my_vecs.push_back(S.GetPtr<CompositeVector>(key_tag.first, key_tag.second).ptr());
    }
    db_->WriteVectors(my_names, my_vecs);
  }
}



} // namespace Amanzi
