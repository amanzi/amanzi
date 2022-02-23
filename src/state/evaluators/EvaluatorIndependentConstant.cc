/*
  State

  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

  A field evaluator with no dependencies specified by a constant value.
*/

#include "EvaluatorIndependentConstant.hh"

namespace Amanzi {

// ---------------------------------------------------------------------------
// Constructor
// ---------------------------------------------------------------------------
EvaluatorIndependentConstant::EvaluatorIndependentConstant(
    Teuchos::ParameterList& plist)
  : EvaluatorIndependent<CompositeVector, CompositeVectorSpace>(plist),
    value_(plist.get<double>("value"))
{
  temporally_variable_ = false;
};


// ---------------------------------------------------------------------------
// Virtual Copy constructor
// ---------------------------------------------------------------------------
Teuchos::RCP<Evaluator> EvaluatorIndependentConstant::Clone() const {
  return Teuchos::rcp(new EvaluatorIndependentConstant(*this));
}


// ---------------------------------------------------------------------------
// Update the value in the state.
// ---------------------------------------------------------------------------
void EvaluatorIndependentConstant::Update_(State& S) {
  CompositeVector& cv = S.GetW<CompositeVector>(my_key_, my_tag_, my_key_);
  cv.PutScalar(value_);
}

} // namespace Amanzi
