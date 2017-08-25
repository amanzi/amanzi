/*
  State

  Author: Ethan Coon

  A field evaluator with no dependencies specified by a function.
*/

#include "UniqueHelpers.hh"
#include "CompositeVectorFunctionFactory.hh"
#include "EvaluatorIndependentFunction.hh"

namespace Amanzi {

// ---------------------------------------------------------------------------
// Constructor
// ---------------------------------------------------------------------------
EvaluatorIndependentFunction::EvaluatorIndependentFunction(
    Teuchos::ParameterList& plist) :
    EvaluatorIndependent(plist) {}


// ---------------------------------------------------------------------------
// Virtual Copy constructor
// ---------------------------------------------------------------------------
Teuchos::RCP<Evaluator> EvaluatorIndependentFunction::Clone() const {
  return Teuchos::rcp(new EvaluatorIndependentFunction(*this));
}

// ---------------------------------------------------------------------------
// Operator=
// ---------------------------------------------------------------------------
Evaluator&
EvaluatorIndependentFunction::operator=(const Evaluator& other) {
  if (this != &other) {
    const EvaluatorIndependentFunction* other_p =
        dynamic_cast<const EvaluatorIndependentFunction*>(&other);
    ASSERT(other_p != NULL);
    *this = *other_p;
  }
  return *this;
}

EvaluatorIndependentFunction&
EvaluatorIndependentFunction::operator=(const EvaluatorIndependentFunction& other) {
  if (this != &other) {
    ASSERT(my_key_ == other.my_key_);
    requests_ = other.requests_;
  }
  return *this;
}


// ---------------------------------------------------------------------------
// Update the value in the state.
// ---------------------------------------------------------------------------
void EvaluatorIndependentFunction::Update_(State& S) {
  if (!computed_once_) {
    // Create the function.
    auto& cv = S.Get<CompositeVector>(my_key_, my_tag_);
    ASSERT(plist_.isSublist("function"));
    func_ = Functions::CreateCompositeVectorFunction(plist_.sublist("function"),
            cv.Map());
  }

  // NOTE: EvaluatorIndependentFunctions own their own data.
  CompositeVector& cv =
      S.GetW<CompositeVector>(my_key_, my_tag_, my_key_);
  time_ = S.time();
  func_->Compute(time_, cv);
}


} // namespace
