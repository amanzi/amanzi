/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*
  State

  A field evaluator with no dependencies specified by a function.
*/

#include <memory>

#include "EvaluatorIndependentPatchFunction.hh"

namespace Amanzi {

// ---------------------------------------------------------------------------
// Constructor
// ---------------------------------------------------------------------------
EvaluatorIndependentPatchFunction::EvaluatorIndependentPatchFunction(Teuchos::ParameterList& plist)
  : EvaluatorIndependent<MultiPatch<double>, MultiPatchSpace>(plist)
{
  function_outer_name_ = plist_.get<std::string>("function list name", "function");
  function_inner_name_ = plist_.get<std::string>("function inner list name", "function");
};


// ---------------------------------------------------------------------------
// Virtual Copy constructor
// ---------------------------------------------------------------------------
Teuchos::RCP<Evaluator>
EvaluatorIndependentPatchFunction::Clone() const
{
  return Teuchos::rcp(new EvaluatorIndependentPatchFunction(*this));
}


// ---------------------------------------------------------------------------
// Operator=
// ---------------------------------------------------------------------------
Evaluator&
EvaluatorIndependentPatchFunction::operator=(const Evaluator& other)
{
  if (this != &other) {
    const EvaluatorIndependentPatchFunction* other_p =
      dynamic_cast<const EvaluatorIndependentPatchFunction*>(&other);
    AMANZI_ASSERT(other_p != NULL);
    *this = *other_p;
  }
  return *this;
}


EvaluatorIndependentPatchFunction&
EvaluatorIndependentPatchFunction::operator=(const EvaluatorIndependentPatchFunction& other)
{
  if (this != &other) {
    AMANZI_ASSERT(my_key_ == other.my_key_);
    requests_ = other.requests_;
  }
  return *this;
}


void
EvaluatorIndependentPatchFunction::EnsureCompatibility(State& S)
{
  // EvaluatorIndependent::EnsureCompatibility sets flags, requires and claims data
  EvaluatorIndependent<MultiPatch<double>, MultiPatchSpace>::EnsureCompatibility(S);

  // but unlike EvaluatorIndependent function, here we need to ensure we have a
  // MultiPatchSpace that matches the functions that were set in the input
  // spec.
  auto& mps = S.Require<MultiPatch<double>, MultiPatchSpace>(my_key_, my_tag_, my_key_);
  if (func_ == Teuchos::null && mps.mesh != Teuchos::null) {
    AMANZI_ASSERT(mps.size() == 0);
    func_ = Teuchos::rcp(new Functions::MeshFunction(plist_.sublist(function_outer_name_),
            mps.mesh,
            function_inner_name_,
            mps.entity_kind,
            mps.flag_type));
    for (const auto& spec : *func_) mps.addPatch(std::get<1>(spec));
    AMANZI_ASSERT(mps.size() == func_->size());
    AMANZI_ASSERT(mps.size() > 0);
  }
}


// ---------------------------------------------------------------------------
// Update the value in the state.
// ---------------------------------------------------------------------------
void
EvaluatorIndependentPatchFunction::Update_(State& S)
{
  auto& mp = S.GetW<MultiPatch<double>>(my_key_, my_tag_, my_key_);

  // NOTE: EvaluatorIndependentPatchFunctions own their own data.
  time_ = S.get_time(my_tag_);
  func_->Compute(time_, mp);
}

} // namespace Amanzi
