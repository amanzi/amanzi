/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
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
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(plist)
{
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
      << my_keys_[0].second << ": missing list \"function\" or \"functions\"";
    throw(m);
  }

  auto deriv_plist = plist.sublist("derivative functions");
  for (const auto& entry : deriv_plist) {
    if (!deriv_plist.isSublist(entry.first)) {
      Errors::Message msg("EvaluatorSecondaryMonotypeFromFunction: derivative functions entries must be sublists, named by the wrt key.");
    }
    deriv_funcs_[entry.first].push_back(Teuchos::rcp(fac.Create(deriv_plist.sublist(entry.first))));
  }  
}

EvaluatorSecondaryMonotypeFromFunction::EvaluatorSecondaryMonotypeFromFunction(
  const EvaluatorSecondaryMonotypeFromFunction& other)
  : EvaluatorSecondaryMonotype(other)
{
  for (const auto& fp : other.funcs_) { funcs_.emplace_back(fp->Clone()); }
}

Teuchos::RCP<Evaluator>
EvaluatorSecondaryMonotypeFromFunction::Clone() const
{
  return Teuchos::rcp(new EvaluatorSecondaryMonotypeFromFunction(*this));
}

// These do the actual work
void
EvaluatorSecondaryMonotypeFromFunction::Evaluate_(
  const State& S, const std::vector<CompositeVector*>& results)
{
  int ndeps = dependencies_.size();
  std::vector<Teuchos::Ptr<const CompositeVector>> deps;
  for (auto& dep : dependencies_) {
    deps.emplace_back(S.GetPtr<CompositeVector>(dep.first, dep.second).ptr());
  }

  for (auto comp : *results[0]) {
    // Two possible implementations here -- either create the temporary data as
    // a view, copy into that view, and then call apply.  Alternatively, could
    // loop over entities, construct the local "point", then call the
    // operator() instead of apply.  Unclear which is better -- we'll try the
    // first.

    // create temporary
    // space to hold the dependencies
    Kokkos::View<double**> in(
      "tmp_in",
      deps.size(),
      results[0]->GetComponent(comp, false)->getLocalLength());

    std::size_t i = 0;
    for (const auto& dep : deps) {
      if (dep->getNumVectors(comp) != 1) {
        Errors::Message msg(
          "EvaluatorSecondaryMonotypeFromFunction: Currently cannot handle "
          "true MultiVectors, so all dependencies must be MultiVectors with "
          "only 1 vector.");
        throw(msg);
      }
      Kokkos::deep_copy(Kokkos::subview(in, i, Kokkos::ALL),
                        dep->ViewComponent(comp, 0, false));
      i++;
    }

    // loop over results and evaluate the function
    for (std::size_t j = 0; j != results.size(); ++j) {
      if (results[j]->getNumVectors(comp) != 1) {
        Errors::Message msg("EvaluatorSecondaryMonotypeFromFunction: Currently "
                            "cannot handle true MultiVectors, so all result "
                            "MultiVectors must have only 1 vector.");
        throw(msg);
      }
      Kokkos::View<double*> comp_view =
        results[j]->ViewComponent(comp, 0, false);
      funcs_[j]->apply(in, comp_view);
    }
  }
}


// These do the actual work
void
EvaluatorSecondaryMonotypeFromFunction::EvaluatePartialDerivative_(
  const State& S, const Key& wrt_key, const Key& wrt_tag,
  const std::vector<CompositeVector*>& results)
{
  int ndeps = dependencies_.size();
  std::vector<Teuchos::Ptr<const CompositeVector>> deps;
  for (auto& dep : dependencies_) {
    deps.emplace_back(S.GetPtr<CompositeVector>(dep.first, dep.second).ptr());
  }

  // specifically not using tag, because it can't be in the plist.  Not sure if
  // that is a design flaw yet... --etc
  if (deriv_funcs_.count(wrt_key)) {
    auto& func = deriv_funcs_.at(wrt_key);
  
    for (auto comp : *results[0]) {
      // Two possible implementations here -- either create the temporary data as
      // a view, copy into that view, and then call apply.  Alternatively, could
      // loop over entities, construct the local "point", then call the
      // operator() instead of apply.  Unclear which is better -- we'll try the
      // first.

      // create temporary
      // space to hold the dependencies
      Kokkos::View<double**> in(
          "tmp_in",
          deps.size(),
          results[0]->GetComponent(comp, false)->getLocalLength());

      std::size_t i = 0;
      for (const auto& dep : deps) {
        if (dep->getNumVectors(comp) != 1) {
          Errors::Message msg(
              "EvaluatorSecondaryMonotypeFromFunction: Currently cannot handle "
              "true MultiVectors, so all dependencies must be MultiVectors with "
              "only 1 vector.");
          throw(msg);
        }
        Kokkos::deep_copy(Kokkos::subview(in, i, Kokkos::ALL),
                          dep->ViewComponent(comp, 0, false));
        i++;
      }

      // loop over results and evaluate the function
      for (std::size_t j = 0; j != results.size(); ++j) {
        if (results[j]->getNumVectors(comp) != 1) {
          Errors::Message msg("EvaluatorSecondaryMonotypeFromFunction: Currently "
                  "cannot handle true MultiVectors, so all result "
                  "MultiVectors must have only 1 vector.");
          throw(msg);
        }
        Kokkos::View<double*> comp_view =
            results[j]->ViewComponent(comp, 0, false);
        func[j]->apply(in, comp_view);
      }
    }

  } else {
    for (auto& result : results) {
      result->putScalar(0.);
    }
  }
}

} // namespace Amanzi
