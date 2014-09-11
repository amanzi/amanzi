/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

A field evaluator with no dependencies specified by a function.

------------------------------------------------------------------------- */

#include "composite_vector_function_factory.hh"
#include "independent_variable_field_evaluator_fromfunction.hh"

namespace Amanzi {

// ---------------------------------------------------------------------------
// Constructor
// ---------------------------------------------------------------------------
IndependentVariableFieldEvaluatorFromFunction::IndependentVariableFieldEvaluatorFromFunction(Teuchos::ParameterList& plist) :
    IndependentVariableFieldEvaluator(plist) {}

// ---------------------------------------------------------------------------
// Copy constructor
// ---------------------------------------------------------------------------
IndependentVariableFieldEvaluatorFromFunction::IndependentVariableFieldEvaluatorFromFunction(const IndependentVariableFieldEvaluatorFromFunction& other) :
    IndependentVariableFieldEvaluator(other),
    func_(other.func_) {}


// ---------------------------------------------------------------------------
// Virtual Copy constructor
// ---------------------------------------------------------------------------
Teuchos::RCP<FieldEvaluator> IndependentVariableFieldEvaluatorFromFunction::Clone() const {
  return Teuchos::rcp(new IndependentVariableFieldEvaluatorFromFunction(*this));
}


// ---------------------------------------------------------------------------
// Update the value in the state.
// ---------------------------------------------------------------------------
void IndependentVariableFieldEvaluatorFromFunction::UpdateField_(const Teuchos::Ptr<State>& S) {
  if (!computed_once_) {
    // Create the function.
    Teuchos::RCP<const CompositeVector> cv = S->GetFieldData(my_key_);
    ASSERT(plist_.isSublist("function"));
    func_ = Functions::CreateCompositeVectorFunction(plist_.sublist("function"), cv->Map());
  }

  // NOTE: IndependentVariableFieldEvaluatorFromFunctions own their own data.
  Teuchos::RCP<CompositeVector> cv = S->GetFieldData(my_key_, my_key_);
  time_ = S->time();
  func_->Compute(time_, cv.ptr());
}


} // namespace
