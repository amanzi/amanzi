/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/*
  Copyright 2010-201x held jointly, see COPYRIGHT.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon
*/

//! A secondary variable evaluator which evaluates functions on its
//! dependenecies.

/*!
Uses functions to evaluate arbitrary algebraic functions of its dependencies.

For example, one might write a dependency:

  VARNAME = 0.2 * DEP1 - DEP2 + 3

as:

Example:
..xml:
    <ParameterList name="VARNAME">
      <Parameter name="field evaluator type" type="string" value="secondary
variable from function"/> <Parameter name="evaluator dependencies"
type="Array{string}" value="{DEP1, DEP2}"/> <ParameterList name="function">
        <ParameterList name="function-linear">
          <Parameter name="x0" type="Array(double)" value="{0.0,0.0}" />
          <Parameter name="y0" type="double" value="3." />
          <Parameter name="gradient" type="Array(double)" value="{0.2, -1}" />
        </ParameterList>
      </ParameterList>
    </ParameterList>

Note this is not done by region currently, but could easily be extended to do
so if it was found useful.
*/

#include "EvaluatorSecondaryFromFunction.hh"

#include "Function.hh"
#include "FunctionFactory.hh"

namespace Amanzi {

EvaluatorSecondaryFromFunction::EvaluatorSecondaryFromFunction(
    Teuchos::ParameterList &plist)
    : EvaluatorSecondary<CompositeVector, CompositeVectorSpace>(plist) {
  FunctionFactory fac;
  func_ = Teuchos::rcp(fac.Create(plist.sublist("function")));
}

EvaluatorSecondaryFromFunction::EvaluatorSecondaryFromFunction(
    const EvaluatorSecondaryFromFunction &other)
    : EvaluatorSecondary(other), func_(other.func_->Clone()) {}

Teuchos::RCP<Evaluator> EvaluatorSecondaryFromFunction::Clone() const {
  return Teuchos::rcp(new EvaluatorSecondaryFromFunction(*this));
}

// These do the actual work
void EvaluatorSecondaryFromFunction::Evaluate_(const State &S,
                                               CompositeVector &result) {
  int ndeps = dependencies_.size();
  std::vector<Teuchos::Ptr<const CompositeVector>> deps;
  for (auto &dep : dependencies_) {
    deps.emplace_back(S.GetPtr<CompositeVector>(dep.first, dep.second).ptr());
  }

  for (auto comp : result) {
    Epetra_MultiVector &result_vec = *result.ViewComponent(comp, false);
    std::vector<Teuchos::Ptr<const Epetra_MultiVector>> dep_vecs;
    for (auto &dep : deps) {
      dep_vecs.emplace_back(dep->ViewComponent(comp, false).ptr());
    }
    for (int i = 0; i != result_vec.MyLength(); ++i) {
      std::vector<double> p(ndeps);
      for (int j = 0; j != ndeps; ++j) {
        p[j] = (*dep_vecs[j])[0][i];
      }
      result_vec[0][i] = (*func_)(p);
    }
  }
}

} // namespace Amanzi
