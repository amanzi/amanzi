/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/*
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)
*/
//!  A field evaluator with no dependencies specified by a function.

/*!

This evaluator is typically used for providing data that are functions of space
and time.  The evaluator consists of a list of region,function pairs, and the
functions are evaluated across that region at each timestep.  If the problem is
time-independent, the `"constant in time`" option results in a performance
boost (as the functions need only be evaluated once).  This leverages the
exaustive functional format capability provided in Amanzi's Functions_ library.

This evaluator is used by providing the option:

`"field evaluator type`" == `"independent variable`"

.. _independent-variable-evaluator-spec:
.. admonition:: independent-variable-evaluator-spec

   * `"constant in time`" ``[bool]`` **false** If true, only evaluate the
      functions once as they are time-independent.
   * `"function`" ``[composite-vector-function-spec-list]``

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
