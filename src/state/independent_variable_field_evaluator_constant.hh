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
.. todo:
    This needs a test and documentation! --etc
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
