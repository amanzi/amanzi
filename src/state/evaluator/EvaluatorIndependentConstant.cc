/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//! An evaluator with no dependencies specified by a constant.

/*!

.. todo:
    This needs a test and documentation! --etc

*/

#include "EvaluatorIndependentConstant.hh"

namespace Amanzi {

// ---------------------------------------------------------------------------
// Constructor
// ---------------------------------------------------------------------------
EvaluatorIndependentConstant::EvaluatorIndependentConstant(
  Teuchos::ParameterList& plist)
    : EvaluatorIndependent<CompositeVector, CompositeVectorSpace>(plist),
      value_(plist.get<double>("value")),
      computed_once_(false)
{}

// ---------------------------------------------------------------------------
// Virtual Copy constructor
// ---------------------------------------------------------------------------
Teuchos::RCP<Evaluator>
EvaluatorIndependentConstant::Clone() const
{
  return Teuchos::rcp(new EvaluatorIndependentConstant(*this));
}

// ---------------------------------------------------------------------------
// Operator=
// ---------------------------------------------------------------------------
Evaluator&
EvaluatorIndependentConstant::operator=(const Evaluator& other)
{
  if (this != &other) {
    const EvaluatorIndependentConstant* other_p =
      dynamic_cast<const EvaluatorIndependentConstant*>(&other);
    AMANZI_ASSERT(other_p != NULL);
    *this = *other_p;
  }
  return *this;
}

EvaluatorIndependentConstant&
EvaluatorIndependentConstant::
operator=(const EvaluatorIndependentConstant& other)
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
EvaluatorIndependentConstant::Update_(State& S)
{
  if (!computed_once_) {
    S.GetW<CompositeVector>(my_key_, my_tag_, my_key_).putScalar(value_);
    computed_once_ = true;
  }
}

} // namespace Amanzi
