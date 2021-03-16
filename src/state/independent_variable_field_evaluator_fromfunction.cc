/*
  State

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

  A field evaluator with no dependencies specified by a function.
*/

#include "CompositeVectorFunctionFactory.hh"
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
  if (!func_.get()) {
    // Create the function.
    Teuchos::RCP<const CompositeVector> cv = S->GetFieldData(my_key_);
    AMANZI_ASSERT(plist_.isSublist("function"));
    func_ = Functions::CreateCompositeVectorFunction(plist_.sublist("function"), cv->Map());
  }

  // NOTE: IndependentVariableFieldEvaluatorFromFunctions own their own data.
  Teuchos::RCP<CompositeVector> cv = S->GetFieldData(my_key_, my_key_);
  func_->Compute(time_, cv.ptr());
}


} // namespace
