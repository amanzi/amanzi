/*
  State

  Copyright 2010-202x held jointly, see COPYRIGHT.
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

#include "EvaluatorSecondaryMonotypeFromFunction.hh"

#include "Function.hh"
#include "FunctionFactory.hh"

namespace Amanzi {

EvaluatorSecondaryMonotypeFromFunction::EvaluatorSecondaryMonotypeFromFunction(
    Teuchos::ParameterList& plist)
    : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(plist) {
  FunctionFactory fac;
  if (plist.isSublist("functions")) {
    auto& flist = plist.sublist("functions");
    for (auto& lcv : flist) {
      if (flist.isSublist(lcv.first)) {
        funcs_.push_back(Teuchos::rcp(fac.Create(flist.sublist(lcv.first))));
      }
    }
  } else if (plist.isSublist("function")) {
    funcs_.push_back(Teuchos::rcp(fac.Create(plist.sublist("function"))));
  } else {
    Errors::Message m;
    m << "EvaluatorSecondaryMonotypeFromFunction: " << my_keys_[0].first << ","
      << my_keys_[0].second.get() << ": missing list \"function\" or \"functions\"";
    throw(m);
  }
}


EvaluatorSecondaryMonotypeFromFunction::EvaluatorSecondaryMonotypeFromFunction(
    const EvaluatorSecondaryMonotypeFromFunction& other)
    : EvaluatorSecondaryMonotype(other) {
  // for (const auto& fp : other.funcs_) {
  //   funcs_.emplace_back(fp->Clone());
  // }
  AMANZI_ASSERT(false);
}


Teuchos::RCP<Evaluator> EvaluatorSecondaryMonotypeFromFunction::Clone() const
{
  return Teuchos::rcp(new EvaluatorSecondaryMonotypeFromFunction(*this));
}


// These do the actual work
void EvaluatorSecondaryMonotypeFromFunction::Evaluate_(
    const State& S, const std::vector<CompositeVector*>& results)
{
  int ndeps = dependencies_.size();
  std::vector<Teuchos::Ptr<const CompositeVector>> deps;
  for (auto& dep : dependencies_) {
    deps.emplace_back(S.GetPtr<CompositeVector>(dep.first, dep.second).ptr());
  }

  for (auto comp : *results[0]) {
    std::vector<Teuchos::Ptr<const Epetra_MultiVector>> dep_vecs;
    for (const auto& dep : deps) {
      dep_vecs.emplace_back(dep->ViewComponent(comp, false).ptr());
    }

    std::vector<Teuchos::Ptr<Epetra_MultiVector>> result_vecs;
    for (auto& result : results) {
      result_vecs.emplace_back(result->ViewComponent(comp, false).ptr());
    }

    for (int i = 0; i != result_vecs[0]->MyLength(); ++i) {
      std::vector<double> p(ndeps);
      for (int j = 0; j != ndeps; ++j)
        p[j] = (*dep_vecs[j])[0][i];

      for (int k = 0; k != funcs_.size(); ++k)
        (*result_vecs[k])[0][i] = (*funcs_[k])(p);
    }
  }
}

} // namespace Amanzi
