/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/
/*!

This evaluator is typically used for providing data that are functions of space
and time.  The evaluator consists of a list of (region, function) pairs, and the
functions are evaluated across that region at each timestep.  If the problem is
time-independent, the `"constant in time`" option results in a performance
improvement (as the functions need only be evaluated once).  This leverages the
exaustive functional format capability provided in Amanzi's :ref:`Functions` library.

`"evaluator type`" == `"independent variable`"

.. _evaluator-independent-variable-spec:
.. admonition:: evaluator-independent-variable-spec

   * `"constant in time`" ``[bool]`` **false** If true, only evaluate the
     functions once as they are time-independent.
   * `"function`" ``[composite-vector-function-spec-list]``

*/

#ifndef AMANZI_STATE_INDEPENDENT_EVALUATOR_FUNCTION_
#define AMANZI_STATE_INDEPENDENT_EVALUATOR_FUNCTION_

#include "CompositeVectorFunction.hh"
#include "EvaluatorIndependent.hh"
#include "Evaluator_Factory.hh"

namespace Amanzi {

class EvaluatorIndependentFunction
  : public EvaluatorIndependent<CompositeVector, CompositeVectorSpace> {
 public:
  explicit EvaluatorIndependentFunction(Teuchos::ParameterList& plist);
  EvaluatorIndependentFunction(const EvaluatorIndependentFunction& other) = default;

  virtual Teuchos::RCP<Evaluator> Clone() const override;

  virtual Evaluator& operator=(const Evaluator& other) override;

  EvaluatorIndependentFunction& operator=(const EvaluatorIndependentFunction& other);

 protected:
  virtual void Update_(State& S) override;

 protected:
  Teuchos::RCP<Functions::CompositeVectorFunction> func_;
  bool dot_with_normal_;

 private:
  static Utils::RegisteredFactory<Evaluator, EvaluatorIndependentFunction> fac_;
};

} // namespace Amanzi

#endif
