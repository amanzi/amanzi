/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//! A secondary evaluator based on a user-provided Model, computed on patches.

/*!

  This implements a very generic EvaluatorSecondary which uses a
  user-provided Model.  That model provides the list of keys calculated,
  dependencies, and operator() methods that implement the model (but not yet
  derivatives) to act as a Kokkos functor, on patches.

  This is typically going to be used as an evaluator for nonlinear sources and
  boundary conditions which need a model to compute a flux or primary value.
  Unlike an EvaluatorModelCV, whose computed things are CVs and whose
  dependencies are CVs, here we assume that the dependencies are CVs, but that
  the computed thing is on a MultiPatch.  This is because this is typically
  used for very small subsets of the domain -- e.g. boundary faces for BCs, or
  small regions for sources.
  
  Currently, this class assumes that tags of the keys calculated and
  dependencies are all the same.

  Models are expected to implement a fairly precise set of variables and
  methods.  Examples of these are shown (for a simple synthetic case) in
  `src/state/test/dag_domain_function_models.hh`.

  Use of this evaluator is tested in
  src/state/test/state_domain_function_evaluators.cc.

*/


#pragma once

#include "AmanziTypes.hh"
#include "CompositeVectorSpace.hh"
#include "CompositeVector.hh"
#include "Factory.hh"
#include "StateDefs.hh"
#include "EvaluatorModelLauncher.hh"
#include "EvaluatorSecondary.hh"

namespace Amanzi {


//
// The actual generic Evaluator class.
//
template <template <class, class> class Model, class Device_type = DefaultDevice>
class EvaluatorModelPatch : public EvaluatorSecondary {
 public:
  using View_type = Patch<double>::View_type;
  using cView_type = CompositeVector::cView_type;
  using Model_type = Model<cView_type, View_type>;

  EvaluatorModelPatch(const Teuchos::RCP<Teuchos::ParameterList>& plist);

  virtual Teuchos::RCP<Evaluator> Clone() const override;
  virtual std::string getType() const override { return model_->eval_type + " on patch"; }

  virtual bool
  IsDifferentiableWRT(const State& S, const Key& wrt_key, const Tag& wrt_tag) const override
  {
    return false; // derivatives are not implemented
  }

 protected:
  virtual void Update_(State& S) override;

  virtual void UpdateDerivative_(State& S, const Key& wrt_key, const Tag& wrt_tag) override
  {
    AMANZI_ASSERT(false); // never called
  }

  // note, this follows that of EvaluatorSecondaryMonotype, but since this is
  // not monotype, we have to reimplement it here.
  virtual void EnsureCompatibility(State& S) override;

 protected:
  Teuchos::RCP<Model_type> model_;
  std::string name_;

 private:
  // registration in the evaluator factory
  static Utils::RegisteredFactory<Evaluator, EvaluatorModelPatch<Model, Device_type>> reg_;
};


//
// Implementation
// -----------------------------------------------------------------------------
template <template <class, class> class Model, class Device_type>
EvaluatorModelPatch<Model, Device_type>::EvaluatorModelPatch(
  const Teuchos::RCP<Teuchos::ParameterList>& plist)
  : EvaluatorSecondary(plist),
    model_(Teuchos::rcp(new Model_type(plist))),
    name_(Keys::cleanPListName(plist->name()))
{
  for (const KeyTag& dep : model_->getDependencies()) dependencies_.insert(dep);
}

template <template <class, class> class Model, class Device_type>
Teuchos::RCP<Evaluator>
EvaluatorModelPatch<Model, Device_type>::Clone() const
{
  return Teuchos::rcp(new EvaluatorModelPatch<Model, Device_type>(*this));
}

template <template <class, class> class Model, class Device_type>
void
EvaluatorModelPatch<Model, Device_type>::Update_(State& S)
{
  KeyTag my_key = my_keys_.front();
  MultiPatch<double>& result =
    S.GetW<MultiPatch<double>>(my_key.first, my_key.second, my_key.first);

  for (auto& res_patch : result) {
    std::vector<View_type> result_views{ res_patch.data };
    AmanziMesh::Entity_kind entity_kind = res_patch.space->entity_kind;
    std::string comp = to_string(entity_kind);

    std::vector<cView_type> dependency_views;
    for (const auto& dep : dependencies_) {
      const auto& vec = S.Get<CompositeVector>(dep.first, dep.second);
      dependency_views.emplace_back(vec.viewComponent(comp, false));
    }

    // set up the model and range and then dispatch
    model_->setViews(dependency_views, result_views, S);

    auto mat_ids = result.space->mesh->getSetEntities(
      res_patch.space->region, entity_kind, AmanziMesh::Parallel_kind::OWNED);
    model_->setAccessor(mat_ids);

    Kokkos::RangePolicy<typename Device_type::execution_space> range(0, result_views[0].extent(0));

    Kokkos::parallel_for(name_, range, *model_);
    Kokkos::fence();

    // Reset views
    model_ = Teuchos::rcp<Model_type>(new Model_type(plist_));
  }
  //Debug_(S);
}


// note, this follows that of EvaluatorSecondaryMonotype, but since this is
// not monotype, we have to reimplement it here.
template <template <class, class> class Model, class Device_type>
void
EvaluatorModelPatch<Model, Device_type>::EnsureCompatibility(State& S)
{
  // EnsureCompatibility_ClaimOwnership_
  AMANZI_ASSERT(my_keys_.size() == 1);
  auto my_keytag = my_keys_.front();
  auto& my_fac = S.Require<MultiPatch<double>, MultiPatchSpace>(
    my_keytag.first, my_keytag.second, my_keytag.first);

  EnsureCompatibility_Flags_(S);

  // EnsureCompatibility_Deps_, Structure_
  if (my_fac.mesh != Teuchos::null) {
    for (const auto& dep : dependencies_) {
      for (auto& ps : my_fac) {
        S.Require<CompositeVector, CompositeVectorSpace>(dep.first, dep.second)
          .SetMesh(my_fac.mesh)
          ->AddComponent(to_string(ps->entity_kind), ps->entity_kind, ps->num_vectors);
      }
    }
  }
}


} // namespace Amanzi
