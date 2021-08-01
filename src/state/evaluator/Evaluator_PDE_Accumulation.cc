/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon (coonet@ornl.gov)
*/

//! Wraps a PDE_Accumulation to be an Evaluator.

/*!

Lots of options here, document me!

*/

#include "Evaluator_PDE_Accumulation.hh"

namespace Amanzi {

Evaluator_PDE_Accumulation::Evaluator_PDE_Accumulation(
  Teuchos::ParameterList& plist)
    : EvaluatorSecondaryMonotype<CompositeVector,CompositeVectorSpace>(plist)
{
  AMANZI_ASSERT(my_keys_.size() == 1);
  auto domain = Keys::getDomain(my_keys_[0].first);

  conserved_key_ = Keys::readKey(plist, domain, "conserved quantity");
  cv_key_ = Keys::readKey(plist, domain, "cell volume", "cell_volume");
  tag_old_ = plist.get<std::string>("tag old");
  tag_new_ = plist.get<std::string>("tag new");

  dependencies_.emplace_back(std::make_pair(conserved_key_, tag_old_));
  dependencies_.emplace_back(std::make_pair(conserved_key_, tag_new_));
  dependencies_.emplace_back(std::make_pair(cv_key_, tag_old_));
  dependencies_.emplace_back(std::make_pair(cv_key_, tag_new_));
  dependencies_.emplace_back(std::make_pair("time", tag_old_));
  dependencies_.emplace_back(std::make_pair("time", tag_new_));
}

void
Evaluator_PDE_Accumulation::EnsureCompatibility(State& S)
{
  // this is a bit of a hack, but makes life so much easier as dealing with
  // derivatives, etc, can be a bit tricky.  Pop the double dependencies in
  // time, then call EnsureCompatibility(), then push them back on.
  dependencies_.pop_back();
  dependencies_.pop_back();
  EvaluatorSecondaryMonotype<CompositeVector,
                             CompositeVectorSpace>::EnsureCompatibility(S);
  dependencies_.emplace_back(std::make_pair("time", tag_old_));
  dependencies_.emplace_back(std::make_pair("time", tag_new_));
  
  S.Require<double>("time", tag_old_);
  S.Require<double>("time", tag_new_);
  S.RequireEvaluator("time", tag_old_);
  S.RequireEvaluator("time", tag_new_);
}


void
Evaluator_PDE_Accumulation::Evaluate_(
  const State& S, const std::vector<CompositeVector*>& results)
{
  Teuchos::OSTab tab = vo_.getOSTab();
  AMANZI_ASSERT(results.size() == 1);
  double dt = S.Get<double>("time", tag_new_) - S.Get<double>("time", tag_old_);
  if (vo_.os_OK(Teuchos::VERB_EXTREME))
    *vo_.os() << "Evaluating PDE_Accumulation with dt = " << dt << std::endl;
  results[0]->elementWiseMultiply(
    1. / dt,
    S.Get<CompositeVector>(conserved_key_, tag_new_),
    S.Get<CompositeVector>(cv_key_, tag_new_),
    0.);
  results[0]->elementWiseMultiply(
    -1. / dt,
    S.Get<CompositeVector>(conserved_key_, tag_old_),
    S.Get<CompositeVector>(cv_key_, tag_old_),
    1.);
  Debug_(S);
}


void
Evaluator_PDE_Accumulation::EvaluatePartialDerivative_(
  const State& S, const Key& wrt_key, const Key& wrt_tag,
  const std::vector<CompositeVector*>& results)
{
  Teuchos::OSTab tab = vo_.getOSTab();

  AMANZI_ASSERT(results.size() == 1);
  double dt = S.Get<double>("time", tag_new_) - S.Get<double>("time", tag_old_);
  if (vo_.os_OK(Teuchos::VERB_EXTREME))
    *vo_.os() << "Evaluating PDE_Accumulation deriv with dt = " << dt << std::endl;

  if (wrt_key == conserved_key_) {
    if (wrt_tag == tag_old_) {
      (*results[0]) = S.Get<CompositeVector>(cv_key_, tag_old_);
      results[0]->scale(-1. / dt);
    } else if (wrt_tag == tag_new_) {
      (*results[0]) = S.Get<CompositeVector>(cv_key_, tag_new_);
      results[0]->scale(1. / dt);
    } else {
      AMANZI_ASSERT(0);
    }
  } else if (wrt_key == cv_key_) {
    if (wrt_tag == tag_old_) {
      (*results[0]) = S.Get<CompositeVector>(conserved_key_, tag_old_);
      results[0]->scale(-1. / dt);
    } else if (wrt_tag == tag_new_) {
      (*results[0]) = S.Get<CompositeVector>(conserved_key_, tag_new_);
      results[0]->scale(1. / dt);
    } else {
      AMANZI_ASSERT(0);
    }
  } else {
    AMANZI_ASSERT(0);
  }

  if (vo_.os_OK(Teuchos::VERB_EXTREME))
    results[0]->print(*vo_.os());
}

} // namespace Amanzi
