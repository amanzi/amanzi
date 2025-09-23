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

#include "EvaluatorPrimary.hh"
#include "EvaluatorAlias.hh"

namespace Amanzi {

Teuchos::RCP<Evaluator>
EvaluatorAlias::Clone() const
{
  return Teuchos::rcp(new EvaluatorAlias(*this));
}

EvaluatorAlias&
EvaluatorAlias::operator=(const EvaluatorAlias& other)
{
  if (this != &other) {
    AMANZI_ASSERT(my_key_ == other.my_key_);
    AMANZI_ASSERT(target_key_ == other.target_key_);
    AMANZI_ASSERT(my_tag_ == other.my_tag_);
    AMANZI_ASSERT(target_tag_ == other.target_tag_);
    // nothing to do?
  }
  return *this;
}

Evaluator&
EvaluatorAlias::operator=(const Evaluator& other)
{
  if (this != &other) {
    const EvaluatorAlias* other_p = dynamic_cast<const EvaluatorAlias*>(&other);
    AMANZI_ASSERT(other_p != NULL);
    *this = *other_p;
  }
  return *this;
}


KeyTag
EvaluatorAlias::get_target() const
{
  return KeyTag{ target_key_, target_tag_ };
}


bool
EvaluatorAlias::Update(State& S, const Key& request)
{
  return S.GetEvaluator(target_key_, target_tag_).Update(S, request);
}

bool
EvaluatorAlias::UpdateDerivative(State& S,
                                 const Key& requester,
                                 const Key& wrt_key,
                                 const Tag& wrt_tag)
{
  return S.GetEvaluator(target_key_, target_tag_).UpdateDerivative(S, requester, wrt_key, wrt_tag);
}


bool
EvaluatorAlias::IsDependency(const State& S, const Key& key, const Tag& tag) const
{
  if (key == target_key_ && tag == target_tag_) return true;
  if (S.HasEvaluator(target_key_, target_tag_) &&
      S.GetEvaluator(target_key_, target_tag_).IsDependency(S, key, tag)) {
    return true;
  }
  return false;
}


bool
EvaluatorAlias::IsDifferentiableWRT(const State& S, const Key& wrt_key, const Tag& wrt_tag) const
{
  return (wrt_key == my_key_ && wrt_tag == my_tag_) ||
         S.GetEvaluator(target_key_, target_tag_).IsDifferentiableWRT(S, wrt_key, wrt_tag);
}


bool
EvaluatorAlias::ProvidesKey(const Key& key, const Tag& tag) const
{
  return (key == my_key_ && tag == my_tag_);
}


void
EvaluatorAlias::EnsureEvaluators(State& S)
{
  S.RequireEvaluator(target_key_, target_tag_).EnsureEvaluators(S);
}


void
EvaluatorAlias::EnsureCompatibility(State& S)
{
  S.GetEvaluator(target_key_, target_tag_).EnsureCompatibility(S);
  if (S.HasRecord(target_key_, target_tag_)) {
    S.GetRecordSetW(my_key_).AliasRecord(target_tag_, my_tag_);
  }
}


std::string
EvaluatorAlias::WriteToString() const
{
  std::string ret("Alias ");
  ret = ret + Keys::getKey(my_key_, my_tag_) + " --> " + Keys::getKey(target_key_, target_tag_);
  return ret;
}


void
EvaluatorAlias::SetChanged(State& S)
{
  Teuchos::RCP<Evaluator> target_eval = S.GetEvaluatorPtr(target_key_, target_tag_);
  auto target_eval_as_primary = Teuchos::rcp_dynamic_cast<EvaluatorPrimary_>(target_eval);
  if (target_eval_as_primary == Teuchos::null) {
    Errors::Message msg(
      "SetChanged called on an aliased evaluator whose target is not a primary evaluator");
    Exceptions::amanzi_throw(msg);
  }
  target_eval_as_primary->SetChanged();
}


} // namespace Amanzi
