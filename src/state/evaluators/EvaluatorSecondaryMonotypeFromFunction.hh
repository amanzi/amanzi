/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*
  State

*/
/*!
Uses functions to evaluate arbitrary algebraic functions of its dependencies.

For example, one might write a dependency:

  VARNAME = 0.2 * DEP1 - DEP2 + 3

`"evaluator type`" = `"secondary variable from function`"

.. _secondary-variable-from-function-evaluator-spec:
.. admonition:: secondary-variable-from-function-evaluator-spec

   ONE OF:

   * `"functions`" ``[composite-vector-function-spec-list]`` Note this is used
     for multiple Degress of Freedom.

   OR:

   * `"function`" ``[composite-vector-function-spec]`` Used for a single degree
     of freedom.

Example:

.. code:: xml

   <ParameterList name="VARNAME">
     <Parameter name="field evaluator type" type="string" value="algebraic variable from function"/>
     <Parameter name="evaluator dependencies" type="Array{string}" value="{DEP1, DEP2}"/>
     <ParameterList name="function">
       <ParameterList name="function-linear">
         <Parameter name="x0" type="Array(double)" value="{0.0,0.0}" />
         <Parameter name="y0" type="double" value="3." />
         <Parameter name="gradient" type="Array(double)" value="{0.2, -1}" />
       </ParameterList>
     </ParameterList>
   </ParameterList>


*/

#ifndef STATE_EVALUATOR_ALGEBRAIC_FROMFUNCTION_HH_
#define STATE_EVALUATOR_ALGEBRAIC_FROMFUNCTION_HH_

#include <string>
#include <vector>

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

#include "CompositeVector.hh"
#include "CompositeVectorSpace.hh"
#include "EvaluatorSecondaryMonotype.hh"
#include "Evaluator_Factory.hh"
#include "State.hh"

namespace Amanzi {

class Function;

class EvaluatorSecondaryMonotypeFromFunction
  : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
 public:
  explicit EvaluatorSecondaryMonotypeFromFunction(Teuchos::ParameterList& plist);
  EvaluatorSecondaryMonotypeFromFunction(const EvaluatorSecondaryMonotypeFromFunction& other);
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  virtual std::string getType() const override {
    return "secondary from function";
  }

 protected:
  // These do the actual work
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override;

  // This should get some careful thought of the right strategy.  Punting for
  // now --etc
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override
  {
    for (auto& r : results) { r->putScalar(0.0); }
  }

 protected:
  std::vector<Teuchos::RCP<const Function>> funcs_;

 private:
  static Utils::RegisteredFactory<Evaluator, EvaluatorSecondaryMonotypeFromFunction> fac_;
};

} // namespace Amanzi

#endif
