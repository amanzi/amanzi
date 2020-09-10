/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.
  Authors:
      Ethan Coon
*/

//! An evaluator with no dependencies specified by a constant value.

/*!

An independent variable evaluator that is constant in both space and time.  Really, just a number.

This is a quality of life addition -- it is doable as a standard independent
variable evaluator from a function using the function-constant, but it makes
life much easier on users who just want a single value.

This is used by providing a field evaluator type == "independent variable constant"

* `"value`" ``[double]`` Provide the value used in PutScalar() to set the vector's value.

*/


#pragma once

#include "independent_variable_field_evaluator.hh"
#include "FieldEvaluator_Factory.hh"

namespace Amanzi {

class IndependentVariableFieldEvaluatorConstant :
    public IndependentVariableFieldEvaluator {

public:
  // ---------------------------------------------------------------------------
  // Constructors
  // ---------------------------------------------------------------------------
  explicit IndependentVariableFieldEvaluatorConstant(Teuchos::ParameterList& plist);
  IndependentVariableFieldEvaluatorConstant(const IndependentVariableFieldEvaluatorConstant& other) =
    default;

  virtual Teuchos::RCP<FieldEvaluator> Clone() const override;
  virtual void operator=(const FieldEvaluator& other) override;

  void operator=(const IndependentVariableFieldEvaluatorConstant& other);

protected:
  // ---------------------------------------------------------------------------
  // Update the value in the state.
  // ---------------------------------------------------------------------------
  virtual void UpdateField_(const Teuchos::Ptr<State>& S) override;

 protected:
  double value_;
  bool computed_once_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,IndependentVariableFieldEvaluatorConstant> fac_;
};

} // namespace Amanzi
