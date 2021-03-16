/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.
  Authors:
      Ethan Coon
*/

//! An evaluator with no dependencies specified by a constant value.

#include "independent_variable_field_evaluator_constant.hh"

namespace Amanzi {

// ---------------------------------------------------------------------------
// Constructor
// ---------------------------------------------------------------------------
IndependentVariableFieldEvaluatorConstant::IndependentVariableFieldEvaluatorConstant(
  Teuchos::ParameterList& plist)
    : IndependentVariableFieldEvaluator(plist)
{
  temporally_variable_ = false;
}

// ---------------------------------------------------------------------------
// Virtual Copy constructor
// ---------------------------------------------------------------------------
Teuchos::RCP<FieldEvaluator>
IndependentVariableFieldEvaluatorConstant::Clone() const
{
  return Teuchos::rcp(new IndependentVariableFieldEvaluatorConstant(*this));
}


// ---------------------------------------------------------------------------
// Update the value in the state.
// ---------------------------------------------------------------------------
void
IndependentVariableFieldEvaluatorConstant::UpdateField_(
  const Teuchos::Ptr<State>& S)
{
  S->GetField(my_key_, my_key_)->Initialize(plist_);
}


} // namespace Amanzi
