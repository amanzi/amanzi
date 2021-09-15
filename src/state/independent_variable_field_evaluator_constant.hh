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

An independent variable evaluator that is constant in time.  This uses the same
infrastructure as initial conditions for time-independent data.

This is used by providing:

`"field evaluator type`" == `"independent variable constant`"

See the InitialConditions_ description for all parameters.
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

protected:
  // ---------------------------------------------------------------------------
  // Update the value in the state.
  // ---------------------------------------------------------------------------
  virtual void UpdateField_(const Teuchos::Ptr<State>& S) override;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,IndependentVariableFieldEvaluatorConstant> fac_;
};

} // namespace Amanzi
