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

   * `"dot with normal`" ``[bool]`` **false** If true, this expects a 1-dof,
     face-based vector and a N-D function (N is mesh space dim) and will dot
     the normal of the face with the function to compute the field.

   * `"spatial distribution method`" ``[string]`` **none** Describes a strategy
     for mapping the function values into the field.

     - `"none`" is the 1-1 map, just evaluate the function

     - `"volume`" Divides each function spec value by a region total volume.
       This is used when, for instance, the function represents a total {mass,
       energy, etc} flux of a source that is uniformly distributed across the
       region such that the integral is equal to the function value.
     
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
  std::string spatial_dist_method_;

 private:
  static Utils::RegisteredFactory<Evaluator, EvaluatorIndependentFunction> fac_;
};

} // namespace Amanzi

#endif
