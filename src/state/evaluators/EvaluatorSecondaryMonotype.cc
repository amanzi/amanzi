/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*
  State

  A field evaluator base for secondary variables.  A Evaluator is a node in the
  Phalanx-like dependency tree.

  Secondary variable evaluators, such as equations of state, water retention
  evaluators, internal energy evaluators, etc should inherit this class,
  implementing the missing Update_() and UpdateFieldDerivative_() methods.

  Secondary secondary evaluator where all dependencies and this are
  doubles.
*/

#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {

// ---------------------------------------------------------------------------
// Updates the derivative for doubles
// ---------------------------------------------------------------------------
template <>
void
EvaluatorSecondaryMonotype<double>::UpdateDerivative_(State& S,
                                                      const Key& wrt_key,
                                                      const Tag& wrt_tag)
{
  std::vector<double*> results(my_keys_.size());
  int j = 0;
  for (const auto& keytag : my_keys_) {
    results[j] =
      &S.GetDerivativeW<double>(keytag.first, keytag.second, wrt_key, wrt_tag, keytag.first);
    *results[j] = 0.;
    ++j;
  }

  // if provides key, then the result is 1
  if (ProvidesKey(wrt_key, wrt_tag)) {
    auto keytag = std::make_pair(wrt_key, wrt_tag);
    int i = std::find(my_keys_.begin(), my_keys_.end(), keytag) - my_keys_.begin();
    AMANZI_ASSERT(i < my_keys_.size()); // ensured by IsDifferentiableWRT() check previously
    *results[i] = 1.0;
    return;
  }

  // dF/dx = sum_(deps) partial F/ partial dep * ddep/dx + partial F/partial x
  for (auto& dep : dependencies_) {
    if (wrt_key == dep.first && wrt_tag == dep.second) {
      // partial F / partial x
      std::vector<double> tmp_data(my_keys_.size(), 0.);
      std::vector<double*> tmp(my_keys_.size());
      for (int i = 0; i != my_keys_.size(); ++i) { tmp[i] = &tmp_data[i]; }
      EvaluatePartialDerivative_(S, wrt_key, wrt_tag, tmp);
      for (int i = 0; i != my_keys_.size(); ++i) (*results[i]) += tmp_data[i];

    } else if (S.GetEvaluator(dep.first, dep.second).IsDependency(S, wrt_key, wrt_tag)) {
      // partial F / partial dep * ddep/dx
      // note this has already been Updated in the public version of this
      // function
      // -- ddep/dx
      const auto& ddep = S.GetDerivative<double>(dep.first, dep.second, wrt_key, wrt_tag);

      // -- partial F / partial dep
      std::vector<double> tmp_data(my_keys_.size(), 0.);
      std::vector<double*> tmp(my_keys_.size());
      for (int i = 0; i != my_keys_.size(); ++i) { tmp[i] = &tmp_data[i]; }
      EvaluatePartialDerivative_(S, dep.first, dep.second, tmp);

      // sum
      for (int i = 0; i != my_keys_.size(); ++i) (*results[i]) += ddep * tmp_data[i];
    }
  }
}


// ---------------------------------------------------------------------------
// Updates the derivative for CompositeVectors
// ---------------------------------------------------------------------------
template <>
void
EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>::UpdateDerivative_(
  State& S,
  const Key& wrt_key,
  const Tag& wrt_tag)
{
  std::vector<CompositeVector*> results(my_keys_.size());
  int j = 0;
  for (const auto& keytag : my_keys_) {
    results[j] = &S.GetDerivativeW<CompositeVector>(
      keytag.first, keytag.second, wrt_key, wrt_tag, keytag.first);
    results[j]->PutScalarMasterAndGhosted(0.0);
    ++j;
  }

  // if provides key, then the result is 1
  if (ProvidesKey(wrt_key, wrt_tag)) {
    auto keytag = std::make_pair(wrt_key, wrt_tag);
    int i = std::find(my_keys_.begin(), my_keys_.end(), keytag) - my_keys_.begin();
    AMANZI_ASSERT(i < my_keys_.size()); // ensured by IsDifferentiableWRT() check previously
    results[i]->PutScalar(1.);
    return;
  }

  // dF/dx = sum_(deps) partial F/ partial dep * ddep/dx + partial F/partial x
  for (auto& dep : dependencies_) {
    if (wrt_key == dep.first && wrt_tag == dep.second) {
      // partial F / partial x
      std::vector<CompositeVector> tmp_data(my_keys_.size(), *results[0]);
      std::vector<CompositeVector*> tmp(my_keys_.size());
      for (int i = 0; i != my_keys_.size(); ++i) { tmp[i] = &tmp_data[i]; }
      EvaluatePartialDerivative_(S, wrt_key, wrt_tag, tmp);
      for (int i = 0; i != my_keys_.size(); ++i) results[i]->Update(1., tmp_data[i], 1.);

    } else if (S.GetEvaluator(dep.first, dep.second).IsDifferentiableWRT(S, wrt_key, wrt_tag)) {
      // partial F / partial dep * ddep/dx
      // note this has already been Updated in the public version of this
      // function

      // -- ddep/dx
      const auto& ddep = S.GetDerivative<CompositeVector>(dep.first, dep.second, wrt_key, wrt_tag);

      // -- partial F / partial dep
      std::vector<CompositeVector> tmp_data(my_keys_.size(), *results[0]);
      std::vector<CompositeVector*> tmp(my_keys_.size());
      for (int i = 0; i != my_keys_.size(); ++i) { tmp[i] = &tmp_data[i]; }
      EvaluatePartialDerivative_(S, dep.first, dep.second, tmp);

      // sum
      for (int i = 0; i != my_keys_.size(); ++i) results[i]->Multiply(1., ddep, tmp_data[i], 1.);
    }
  }
}


// ---------------------------------------------------------------------------
// If multiple of my_keys are evaluated, it can be useful to make sure that all
// of my keys share the same structure.
// ---------------------------------------------------------------------------
template <>
void
EvaluatorSecondaryMonotype<CompositeVector,
                           CompositeVectorSpace>::EnsureCompatibility_StructureSame_(State& S)
{
  for (const auto& keytag : my_keys_) {
    const auto& my_fac =
      S.Require<CompositeVector, CompositeVectorSpace>(keytag.first, keytag.second);
    for (const auto& keytag_other : my_keys_) {
      if (keytag != keytag_other) {
        auto& other_fac =
          S.Require<CompositeVector, CompositeVectorSpace>(keytag_other.first, keytag_other.second);
        other_fac.Update(my_fac);
      }
    }
  }
}


// ---------------------------------------------------------------------------
// Push metadata from my_keys to my_derivs
//
// Simply Updates deriv facs with fac of the variable.  Assumes the structure
// of my_keys has already been set, either by user code, or by calls to
// dependencies EnsureCompatibility_FromDeps_()
// ---------------------------------------------------------------------------
template <>
void
EvaluatorSecondaryMonotype<CompositeVector,
                           CompositeVectorSpace>::EnsureCompatibility_DerivStructure_(State& S)
{
  for (const auto& keytag : my_keys_) {
    if (S.HasDerivativeSet(keytag.first, keytag.second)) {
      auto my_fac = S.Require<CompositeVector, CompositeVectorSpace>(keytag.first, keytag.second);
      for (const auto& deriv : S.GetDerivativeSet(keytag.first, keytag.second)) {
        auto wrt = Keys::splitKeyTag(deriv.first.get());
        S.RequireDerivative<CompositeVector, CompositeVectorSpace>(
           keytag.first, keytag.second, wrt.first, wrt.second, keytag.first)
          .Update(my_fac);
      }
    }
  }
}


// ---------------------------------------------------------------------------
// Get vector structure for my_keys from dependencies.
// ---------------------------------------------------------------------------
template <>
void
EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>::EnsureCompatibility_FromDeps_(
  State& S,
  const std::string& consistency_policy)
{
  // take my requirements as the intersection or union of my children
  CompositeVectorSpace combined_dep_fac;
  for (const auto& dep : dependencies_) {
    const auto& dep_fac = S.Require<CompositeVector, CompositeVectorSpace>(dep.first, dep.second);
    if (combined_dep_fac.size() == 0) {
      combined_dep_fac = dep_fac;
    } else {
      if (!combined_dep_fac.SameAs(dep_fac)) {
        if ((consistency_policy == "take from child: intersection" &&
             dep_fac.SubsetOf(combined_dep_fac)) ||
            (consistency_policy == "take from child: union" &&
             combined_dep_fac.SubsetOf(dep_fac))) {
          combined_dep_fac = dep_fac;
        }
      }
    }
  }

  // update my facs with this info
  for (auto keytag : my_keys_) {
    S.Require<CompositeVector, CompositeVectorSpace>(keytag.first, keytag.second, keytag.first)
      .Update(combined_dep_fac);
  }
}


// ---------------------------------------------------------------------------
// Push my structure to my dependencies.  Default method, using my factory.
//
// Note, this is the most frequently overridden entry point for other
// evaluators.  Most overrides simply choose a different factory to send to
// their dependencies using the other overload of ToDeps_.  For instance, an
// upwinding evaluator, which looks to calculate coefficents on faces based on
// those from cells and boundary faces, might create a CVS on the same mesh but
// with CELL and BOUNDARY_FACE entries, then call
//
//   EnsureCompatibility_ToDeps_(S, cell_bf_fac);
// ---------------------------------------------------------------------------
template <>
void
EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>::EnsureCompatibility_ToDeps_(
  State& S)
{
  auto akeytag = my_keys_.front();
  const auto& my_fac =
    S.Require<CompositeVector, CompositeVectorSpace>(akeytag.first, akeytag.second);
  EnsureCompatibility_ToDeps_(S, my_fac);
}

// ---------------------------------------------------------------------------
// Push fac structure to dependencies, potentially change the Mesh based on the
// domain name.
// ---------------------------------------------------------------------------
template <>
void
EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>::EnsureCompatibility_ToDeps_(
  State& S,
  const CompositeVectorSpace& fac)
{
  for (const auto& dep : dependencies_) {
    auto& dep_fac = S.Require<CompositeVector, CompositeVectorSpace>(dep.first, dep.second);
    dep_fac.Update(fac);
  }
}


template <>
void
EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>::EnsureCompatibility_ToDeps_(
  State& S,
  const std::vector<std::string>& names,
  const std::vector<AmanziMesh::Entity_kind>& locations,
  const std::vector<int>& num_dofs)
{
  for (const auto& dep : dependencies_) {
    auto& dep_fac = S.Require<CompositeVector, CompositeVectorSpace>(dep.first, dep.second);
    auto domain = Keys::getDomain(dep.first);
    dep_fac.SetMesh(S.GetMesh(domain))->AddComponents(names, locations, num_dofs);
  }
}

template <>
Teuchos::Ptr<const Comm_type>
EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>::get_comm_(const State& S) const
{
  return S.Get<CompositeVector>(my_keys_.front().first, my_keys_.front().second).Comm().ptr();
}


} // namespace Amanzi
