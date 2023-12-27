/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*
  State

  A field evaluator with no dependencies on a patch, specified by a function.
*/

/*!

This evaluator is typically used for providing data that are functions of space
and time.  The evaluator consists of a list of region,function pairs, and the
functions are evaluated across that region at each timestep.  If the problem is
time-independent, the `"constant in time`" option results in a performance
boost (as the functions need only be evaluated once).  This leverages the
exaustive functional format capability provided in Amanzi's Functions_ library.

It evaluates into a MultiPatch object

This evaluator is used by providing the option:

`"evaluator type`" == `"independent variable, patch`"

.. _independent-variable-patch-function-evaluator-spec:
.. admonition:: independent-variable-patch-function-evaluator-spec

   * `"constant in time`" ``[bool]`` **false** If true, only evaluate the
     functions once as they are time-independent.
   * `"function`" ``[mesh-function-spec-list]``

*/

#pragma once

#include "MeshFunction.hh"
#include "EvaluatorIndependent.hh"
#include "Evaluator_Factory.hh"

namespace Amanzi {

class EvaluatorIndependentPatchFunction
  : public EvaluatorIndependent<MultiPatch<double>, MultiPatchSpace> {
 public:
  explicit EvaluatorIndependentPatchFunction(const Teuchos::RCP<Teuchos::ParameterList>& plist);
  EvaluatorIndependentPatchFunction(const EvaluatorIndependentPatchFunction& other) = default;

  virtual Teuchos::RCP<Evaluator> Clone() const override;

  virtual Evaluator& operator=(const Evaluator& other) override;

  EvaluatorIndependentPatchFunction& operator=(const EvaluatorIndependentPatchFunction& other);

  virtual std::string getType() const override { return "independent variable patch"; }

  virtual void EnsureCompatibility(State& S) override;

 protected:
  virtual void Update_(State& S) override;

 protected:
  Teuchos::RCP<Functions::MeshFunction> func_;
  std::string function_outer_name_, function_inner_name_;
  AmanziMesh::Entity_kind entity_kind_;

 private:
  static Utils::RegisteredFactory<Evaluator, EvaluatorIndependentPatchFunction> fac_;
};

} // namespace Amanzi
