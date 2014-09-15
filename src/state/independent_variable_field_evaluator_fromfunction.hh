/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

A field evaluator with no dependencies specified by a function.

------------------------------------------------------------------------- */

#ifndef AMANZI_INDEPENDENT_FIELD_EVALUATOR_FROMFUNCTION_
#define AMANZI_INDEPENDENT_FIELD_EVALUATOR_FROMFUNCTION_

#include "composite_vector_function.hh"
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
