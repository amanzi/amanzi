/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//! A secondary evaluator based on a user-provided Model.

/*!

  This implements a very generic EvaluatorSecondaryMonotype which uses a
  user-provided Model.  The implementation of this is a bit tricky, mainly
  because of the need for a map from the variable that we wish to differentiate
  with to a functor for that partial derivative.  Obvious implementations are
  not able to take a run-time partial derivative dependency to a compile-time
  function in a generic way (with an arbitrary number of potential
  dependencies).

  This solution uses an auxilary class to do a compile-time "loop" (through
  recursion and partial class specialization) to generate, at compile-time, the
  correct Tags needed to execute the correct derivative operation.  Each
  Launcher object's launch() call checks if wrt is dependency I, (starting at
  n_dependencies - 1), and calls Kokkos for_each if so, or recurses to I-1 if
  not.  When I-1 == -1, then wrt is not in dependencies, and it throws an error
  (ending the recursion).  Deriv<I> is the tag corresponding to dependency I,
  as indexed in the Model's list of dependencies.

  Models are expected to implement a fairly precise set of variables and
  methods.  Examples of these are shown (for a simple synthetic case) in
  `src/state/test/dag_models.hh`.

  Use of this evaluator is tested in
  src/state/test/state_evaluators_dag_models.cc.

*/

#ifndef STATE_EVALUATOR_MODEL_BY_MATERIAL_HH_
#define STATE_EVALUATOR_MODEL_BY_MATERIAL_HH_

#include "AmanziTypes.hh"
#include "StateDefs.hh"
#include "EvaluatorModelLauncher.hh"
#include "EvaluatorSecondaryMonotype.hh"
#include "Factory.hh"

namespace Amanzi {

//
// A generic evaluator for acting on regions
//
template <template <class, class> class Model, class Device_type = DefaultDevice>
class EvaluatorModelCVByMaterial
  : public EvaluatorSecondaryMonotypeCV {
 public:
  using cView_type = cMultiVectorView_type<Device_type>;
  using View_type = MultiVectorView_type<Device_type>;
  using Model_type = Model<cView_type, View_type>;
  static const std::string eval_type;

  EvaluatorModelCVByMaterial(const Teuchos::RCP<Teuchos::ParameterList>& plist);

  virtual Teuchos::RCP<Evaluator> Clone() const override;
  virtual std::string getType() const override
  {
    return eval_type;
  }

  std::vector<std::pair<std::string, Teuchos::RCP<Model_type>>>& getModels() { return models_; }

  // This function needs to be public for Kokkos::CUDA backend
  // Function calling kernel cannot be protected/private
 public:
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override;

 protected:
  virtual void EnsureCompatibility_Structure_(State& S) override
  {
    // make sure all my_keys have the same CVS
    EnsureCompatibility_StructureSame_(S);
  }

 protected:
  std::vector<std::pair<std::string, Teuchos::RCP<Model_type>>> models_;
  std::string name_;

 private:
  // registration in the evaluator factory
  static Utils::RegisteredFactory<Evaluator, EvaluatorModelCVByMaterial<Model, Device_type>> reg_;
};


template <template <class, class> class Model, class Device_type>
EvaluatorModelCVByMaterial<Model, Device_type>::EvaluatorModelCVByMaterial(
  const Teuchos::RCP<Teuchos::ParameterList>& plist)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(plist),
    name_(Keys::cleanPListName(plist->name()))
{
  const Teuchos::ParameterList& region_plist = plist->sublist("model parameters");

  for (const auto& region : region_plist) {
    std::string region_name = region.first;
    if (region_plist.isSublist(region_name)) {
      // make a copy for use in the Model
      auto model_list = Teuchos::rcp(new Teuchos::ParameterList(*plist));
      model_list->set("model parameters", region_plist.sublist(region_name));
      model_list->sublist("model parameters").set<std::string>("region", region_name);
      models_.emplace_back(std::make_pair(region_name, Teuchos::rcp(new Model_type(model_list))));
    } else {
      Errors::Message msg(
        "EvaluatorModelCVByMaterial: \"model parameters\" sublist should only contain "
        "sublists, one per material region.");
      throw(msg);
    }
  }

  if (models_.size() == 0) {
    Errors::Message msg("EvaluatorModelCVByMaterial: \"model parameters\" sublist was "
                        "empty, must have at least one model.");
    throw(msg);
  }

  // Take my_keys and dependencies from the first model.  Could check that they
  // are all the same?  Shouldn't be necessary since they all got constructed
  // from the same list.
  my_keys_.clear();
  my_keys_ = models_[0].second->getMyKeys();
  for (const KeyTag& dep : models_[0].second->getDependencies()) dependencies_.insert(dep);
}


template <template <class, class> class Model, class Device_type>
Teuchos::RCP<Evaluator>
EvaluatorModelCVByMaterial<Model, Device_type>::Clone() const
{
  return Teuchos::rcp(new EvaluatorModelCVByMaterial<Model, Device_type>(*this));
}

template <template <class, class> class Model, class Device_type>
void
EvaluatorModelCVByMaterial<Model, Device_type>::Evaluate_(
  const State& S,
  const std::vector<CompositeVector*>& results)
{
  AMANZI_ASSERT(results.size() > 0);

  for (const auto& comp : *results[0]) {
    // get the list of dependency views
    std::vector<cView_type> dependency_views;
    for (const auto& dep : dependencies_) {
      const auto& vec = S.Get<CompositeVector>(dep.first, dep.second);
      dependency_views.emplace_back(vec.viewComponent(comp, false));
    }

    // get the list of result views
    std::vector<View_type> result_views;
    for (auto result : results) {
      AMANZI_ASSERT(result->hasComponent(comp));
      result_views.emplace_back(result->viewComponent(comp, false));
    }

    // get the location on which we are executing
    auto mesh_location = results[0]->getMap()->getLocation(comp);

    // loop over region/model pairs and execute the model
    for (const auto& region_model : models_) {
      // get a view of the entities to execute over
      auto mat_ids = results[0]->getMap()->getMesh()->getSetEntities(
        region_model.first, mesh_location, AmanziMesh::Parallel_kind::OWNED);

      // set up the model and range and then dispatch
      region_model.second->setViews(dependency_views, result_views, S);
      auto&& f = *(region_model.second);

      Kokkos::RangePolicy<typename Device_type::execution_space> range(0, mat_ids.extent(0));

      Kokkos::parallel_for(
        name_, range, KOKKOS_LAMBDA(const int& i) { f(mat_ids(i)); });

      // must clear the views -- this tells tpetra's dual view that it can sync
      // if needed.
      region_model.second->freeViews();
    }
  }
  //  Debug_(S);
}

template <template <class, class> class Model, class Device_type>
void
EvaluatorModelCVByMaterial<Model, Device_type>::EvaluatePartialDerivative_(
  const State& S,
  const Key& wrt_key,
  const Tag& wrt_tag,
  const std::vector<CompositeVector*>& results)
{
  AMANZI_ASSERT(results.size() > 0);

  for (const auto& comp : *results[0]) {
    // get the list of dependency views
    std::vector<cView_type> dependency_views;
    for (const auto& dep : dependencies_) {
      const auto& vec = S.Get<CompositeVector>(dep.first, dep.second);
      if (vec.getNumVectors(comp) != 1) {
        Errors::Message msg("EvaluatorModel: expects only dependencies with one DoF.");
        throw(msg);
      }
      dependency_views.emplace_back(vec.viewComponent(comp, false));
    }

    // get the list of results views
    std::vector<View_type> result_views;
    for (auto result : results) {
      if (result->getNumVectors(comp) != 1) {
        Errors::Message msg("EvaluatorModel: expects only results with one DoF.");
        throw(msg);
      }
      result_views.emplace_back(result->viewComponent(comp, false));
    }

    // get the location on which we are executing
    auto mesh_location = results[0]->getMap()->getLocation(comp);

    // loop over region/model pairs and execute the model
    for (const auto& region_model : models_) {
      // get a view of the entities to execute over
      auto mat_ids = results[0]->getMap()->getMesh()->getSetEntities(
        region_model.first, mesh_location, AmanziMesh::Parallel_kind::OWNED);

      // set up the model and range and then dispatch
      region_model.second->setViews(dependency_views, result_views, S);

      KeyTag wrt{ wrt_key, wrt_tag };
      Impl::EvaluatorModelLauncher<Model_type::n_dependencies - 1, Model_type, Device_type>
        launcher(name_, wrt, dependencies_, *region_model.second);
      launcher.launch(mat_ids);

      region_model.second->freeViews();
    }
  }
}


template <template <class, class> class Model, class Device_type>
const std::string EvaluatorModelCVByMaterial<Model, Device_type>::eval_type =
  EvaluatorModelCVByMaterial<Model, Device_type>::Model_type::eval_type + " by material";
} // namespace Amanzi


#endif
