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
    : IndependentVariableFieldEvaluator(plist),
      computed_once_(false)
{}

// ---------------------------------------------------------------------------
// Virtual Copy constructor
// ---------------------------------------------------------------------------
Teuchos::RCP<FieldEvaluator>
IndependentVariableFieldEvaluatorConstant::Clone() const
{
  return Teuchos::rcp(new IndependentVariableFieldEvaluatorConstant(*this));
}

// ---------------------------------------------------------------------------
// Operator=
// ---------------------------------------------------------------------------
void
IndependentVariableFieldEvaluatorConstant::operator=(const FieldEvaluator& other)
{
  if (this != &other) {
    const IndependentVariableFieldEvaluatorConstant* other_p =
      dynamic_cast<const IndependentVariableFieldEvaluatorConstant*>(&other);
    AMANZI_ASSERT(other_p != NULL);
    *this = *other_p;
  }
}

void
IndependentVariableFieldEvaluatorConstant::
operator=(const IndependentVariableFieldEvaluatorConstant& other)
{
  if (this != &other) {
    AMANZI_ASSERT(my_key_ == other.my_key_);
    requests_ = other.requests_;
  }
}


// ---------------------------------------------------------------------------
// Update the value in the state.
// ---------------------------------------------------------------------------
void
IndependentVariableFieldEvaluatorConstant::UpdateField_(
  const Teuchos::Ptr<State>& S)
{
  if (!computed_once_) {
    S->GetField(my_key_, my_key_)->Initialize(plist_);
    computed_once_ = true;
  }
}


} // namespace Amanzi
