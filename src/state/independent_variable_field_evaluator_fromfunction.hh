/*
  State

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

  A field evaluator with no dependencies specified by a function.
*/

#ifndef AMANZI_INDEPENDENT_FIELD_EVALUATOR_FROMFUNCTION_
#define AMANZI_INDEPENDENT_FIELD_EVALUATOR_FROMFUNCTION_

#include "CompositeVectorFunction.hh"
#include "independent_variable_field_evaluator.hh"
#include "FieldEvaluator_Factory.hh"

namespace Amanzi {

class IndependentVariableFieldEvaluatorFromFunction :
      public IndependentVariableFieldEvaluator
{

 public:

  // ---------------------------------------------------------------------------
  // Constructors
  // ---------------------------------------------------------------------------
  explicit
  IndependentVariableFieldEvaluatorFromFunction(Teuchos::ParameterList& plist);
  IndependentVariableFieldEvaluatorFromFunction(const IndependentVariableFieldEvaluatorFromFunction& other);

  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

 protected:
  // ---------------------------------------------------------------------------
  // Update the value in the state.
  // ---------------------------------------------------------------------------
  virtual void UpdateField_(const Teuchos::Ptr<State>& S);

 protected:
  Teuchos::RCP<Functions::CompositeVectorFunction> func_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,IndependentVariableFieldEvaluatorFromFunction> fac_;
};

} // namespace


#endif
