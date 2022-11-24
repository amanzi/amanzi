/*
  State

  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

  A field evaluator with no dependencies specified by a function.
*/

#include <memory>

#include "EvaluatorIndependentFunction.hh"
#include "CompositeVectorFunctionFactory.hh"

namespace Amanzi {

// ---------------------------------------------------------------------------
// Constructor
// ---------------------------------------------------------------------------
EvaluatorIndependentFunction::EvaluatorIndependentFunction(Teuchos::ParameterList& plist)
  : EvaluatorIndependent<CompositeVector, CompositeVectorSpace>(plist){};


// ---------------------------------------------------------------------------
// Virtual Copy constructor
// ---------------------------------------------------------------------------
Teuchos::RCP<Evaluator>
EvaluatorIndependentFunction::Clone() const
{
  return Teuchos::rcp(new EvaluatorIndependentFunction(*this));
}


// ---------------------------------------------------------------------------
// Operator=
// ---------------------------------------------------------------------------
Evaluator&
EvaluatorIndependentFunction::operator=(const Evaluator& other)
{
  if (this != &other) {
    const EvaluatorIndependentFunction* other_p =
      dynamic_cast<const EvaluatorIndependentFunction*>(&other);
    AMANZI_ASSERT(other_p != NULL);
    *this = *other_p;
  }
  return *this;
}


EvaluatorIndependentFunction&
EvaluatorIndependentFunction::operator=(const EvaluatorIndependentFunction& other)
{
  if (this != &other) {
    AMANZI_ASSERT(my_key_ == other.my_key_);
    requests_ = other.requests_;
  }
  return *this;
}


// ---------------------------------------------------------------------------
// Update the value in the state.
// ---------------------------------------------------------------------------
void
EvaluatorIndependentFunction::Update_(State& S)
{
  if (!computed_once_) {
    // Create the function.
    auto& cv = S.Get<CompositeVector>(my_key_, my_tag_);
    AMANZI_ASSERT(plist_.isSublist("function"));

    std::vector<std::string> complist;
    func_ =
      Functions::CreateCompositeVectorFunction(plist_.sublist("function"), cv.Map(), complist);
  }

  // NOTE: EvaluatorIndependentFunctions own their own data.
  auto cv = S.GetPtrW<CompositeVector>(my_key_, my_tag_, my_key_);
  time_ = S.get_time(my_tag_);
  func_->Compute(time_, cv.ptr());
}

} // namespace Amanzi
