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

#include "EvaluatorIndependentFunction.hh"

namespace Amanzi {

const std::string EvaluatorIndependentFunction::name = "independent variable";

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
  auto cv = S.GetPtrW<CompositeVector>(my_key_, my_tag_, my_key_);

  if (!computed_once_) {
    // Create the function.
    AMANZI_ASSERT(plist_.isSublist("function"));
    func_ = Teuchos::rcp(new Functions::CompositeVectorFunction(plist_.sublist("function"), cv->getMesh()));
  }

  // NOTE: EvaluatorIndependentFunctions own their own data.
  time_ = S.get_time(my_tag_);
  func_->Compute(time_, *cv);
}

} // namespace Amanzi
