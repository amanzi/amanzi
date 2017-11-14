/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.
  See $ATS_DIR/COPYRIGHT

  Author: Ethan Coon
*/

//! A secondary variable evaluator which evaluates functions on its dependenecies.

/*!
Uses functions to evaluate arbitrary algebraic functions of its dependencies.

For example, one might write a dependency:

  VARNAME = 0.2 * DEP1 - DEP2 + 3

as:

Note this is not done by region currently, but could easily be extended to do
so if it was found useful.

Example:
..xml:
    <ParameterList name="VARNAME">
      <Parameter name="field evaluator type" type="string" value="secondary variable from function"/>
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

#ifndef STATE_SECONDARY_VARIABLE_FIELD_EVALUATOR_FROMFILE_HH_
#define STATE_SECONDARY_VARIABLE_FIELD_EVALUATOR_FROMFILE_HH_

#include <string>
#include <vector>

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

#include "State.hh"
#include "FieldEvaluator.hh"
#include "FieldEvaluator_Factory.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {

class Function;

class SecondaryVariableFieldEvaluatorFromFunction : public SecondaryVariableFieldEvaluator {

 public:
  explicit
  SecondaryVariableFieldEvaluatorFromFunction(Teuchos::ParameterList& plist);
  SecondaryVariableFieldEvaluatorFromFunction(const SecondaryVariableFieldEvaluatorFromFunction& other);
  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

 protected:
  // These do the actual work
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);

  // This should get some careful thought of the right strategy.  Punting for now --etc
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {
    result->PutScalar(0.);
  }

 protected:
  Teuchos::RCP<Function> func_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,SecondaryVariableFieldEvaluatorFromFunction> fac_;

}; // class FieldEvaluator

} // namespace Amanzi

#endif
