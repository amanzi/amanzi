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

namespace Amanzi {

//
// A generic evaluator for acting on regions
//
template <template <class, class> class Model,
          class Device_type = AmanziDefaultDevice>
class EvaluatorModelByMaterial
  : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
 public:
  using Model_type =
    Model<cVectorView_type<Device_type>, VectorView_type<Device_type>>;
  using View_type = VectorView_type<Device_type>;
  using cView_type = cVectorView_type<Device_type>;

  EvaluatorModelByMaterial(Teuchos::ParameterList& plist);

  virtual Teuchos::RCP<Evaluator> Clone() const override;


// This function needs to be public for Kokkos::CUDA backend
// Function calling kernel cannot be protected/private 
public:
  virtual void Evaluate_(const State& S,
                         const std::vector<CompositeVector*>& results) override;

protected:
  virtual void EvaluatePartialDerivative_(
    const State& S, const Key& wrt_key, const Key& wrt_tag,
    const std::vector<CompositeVector*>& results) override;

  std::vector<std::pair<std::string, Teuchos::RCP<Model_type>>> models_;
  std::string name_;
  Key tag_;

 private:
  // registration in the evaluator factory
  static Utils::RegisteredFactory<Evaluator, EvaluatorModelByMaterial<Model,Device_type>> fac_;
};


template <template <class, class> class Model, class Device_type>
EvaluatorModelByMaterial<Model, Device_type>::EvaluatorModelByMaterial(
  Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(plist),
    name_(plist.name()),
    tag_(plist.get<std::string>("tag"))
{
  Teuchos::ParameterList& region_plist = plist.sublist("model parameters");

  // make a copy for use in the Model
  Teuchos::ParameterList model_list(plist);

  for (const auto& region : region_plist) {
    std::string region_name = region.first;
    if (region_plist.isSublist(region_name)) {
      model_list.set("model parameters", region_plist.sublist(region_name));
      models_.push_back(
        std::make_pair(region_name, Teuchos::rcp(new Model_type(model_list))));
    } else {
      Errors::Message msg(
        "EvaluatorModelByMaterial: \"model parameters\" sublist should only contain "
        "sublists, one per material region.");
      throw(msg);
    }
  }

  if (models_.size() == 0) {
    Errors::Message msg("EvaluatorModelByMaterial: \"regions\" sublist was "
                        "empty, must have at least one model.");
    throw(msg);
  }

  // Take the dependencies from the first model.  Could check that they
  // are all the same?  Shouldn't be necessary since they all got constructed
  // from the same list.
  auto dep_list = models_[0].second->dependencies();
  for (const auto& dep : dep_list) {
    dependencies_.emplace_back(KeyPair(dep, tag_));
  }
}


template <template <class, class> class Model, class Device_type>
Teuchos::RCP<Evaluator>
EvaluatorModelByMaterial<Model, Device_type>::Clone() const
{
  return Teuchos::rcp(new EvaluatorModelByMaterial<Model, Device_type>(*this));
}

template <template <class, class> class Model, class Device_type>
void
EvaluatorModelByMaterial<Model, Device_type>::Evaluate_(
  const State& S, const std::vector<CompositeVector*>& results)
{
  AMANZI_ASSERT(results.size() > 0);
  for (const auto& comp : *results[0]) {
    // get the list of dependency views
    std::vector<cView_type> dependency_views;
    for (const auto& dep : dependencies_) {
      const auto& vec = S.Get<CompositeVector>(dep.first, dep.second);
      if (vec.getNumVectors(comp) != 1) {
        Errors::Message msg(
          "EvaluatorModel: expects only dependencies with one DoF.");
        throw(msg);
      }
      dependency_views.emplace_back(vec.ViewComponent(comp, 0, false));
    }

    // get the list of result views
    std::vector<View_type> result_views;
    for (auto result : results) {
      if (result->getNumVectors(comp) != 1) {
        Errors::Message msg(
          "EvaluatorModel: expects only results with one DoF.");
        throw(msg);
      }
      result_views.emplace_back(result->ViewComponent(comp, 0, false));
    }

    // get the location on which we are executing
    auto mesh_location = results[0]->getMap()->Location(comp);
    AMANZI_ASSERT(
      mesh_location !=
      AmanziMesh::Entity_kind::BOUNDARY_FACE); // not yet implemented

    // loop over region/model pairs and execute the model
    for (const auto& region_model : models_) {
      // get a view of the entities to execute over
      AmanziMesh::Entity_ID_View material_ids;
      results[0]->getMap()->Mesh()->get_set_entities(
        region_model.first,
        mesh_location,
        AmanziMesh::Parallel_type::OWNED,
        material_ids);

      // set up the model and range and then dispatch
      region_model.second->SetViews(dependency_views, result_views);

      Kokkos::RangePolicy<typename Device_type::execution_space> range(
        0, material_ids.extent(0));
      Kokkos::parallel_for(name_, range, KOKKOS_LAMBDA(const int i) {
        (*region_model.second)(material_ids(i));
      });
    }
  }
}

template <template <class, class> class Model, class Device_type>
void
EvaluatorModelByMaterial<Model, Device_type>::EvaluatePartialDerivative_(
  const State& S, const Key& wrt_key, const Key& wrt_tag,
  const std::vector<CompositeVector*>& results)
{
  AMANZI_ASSERT(results.size() > 0);
  auto wrt = std::make_pair(wrt_key, wrt_tag);

  for (const auto& comp : *results[0]) {
    // get the list of dependency views
    std::vector<cView_type> dependency_views;
    for (const auto& dep : dependencies_) {
      const auto& vec = S.Get<CompositeVector>(dep.first, dep.second);
      if (vec.getNumVectors(comp) != 1) {
        Errors::Message msg(
          "EvaluatorModel: expects only dependencies with one DoF.");
        throw(msg);
      }
      dependency_views.emplace_back(vec.ViewComponent(comp, 0, false));
    }

    // get the list of results views
    std::vector<View_type> result_views;
    for (auto result : results) {
      if (result->getNumVectors(comp) != 1) {
        Errors::Message msg(
          "EvaluatorModel: expects only results with one DoF.");
        throw(msg);
      }
      result_views.emplace_back(result->ViewComponent(comp, 0, false));
    }

    // get the location on which we are executing
    auto mesh_location = results[0]->getMap()->Location(comp);
    AMANZI_ASSERT(
      mesh_location !=
      AmanziMesh::Entity_kind::BOUNDARY_FACE); // not yet implemented

    // loop over region/model pairs and execute the model
    for (const auto& region_model : models_) {
      // get a view of the entities to execute over
      AmanziMesh::Entity_ID_View material_ids;
      results[0]->getMap()->Mesh()->get_set_entities(
        region_model.first,
        mesh_location,
        AmanziMesh::Parallel_type::OWNED,
        material_ids);

      // set up the model and range and then dispatch
      region_model.second->SetViews(dependency_views, result_views);

      Impl::EvaluatorModelLauncher<Model_type::n_dependencies - 1,
                                   Model_type,
                                   Device_type>
        launcher(name_, wrt, dependencies_, *region_model.second);
      launcher.launch(material_ids);
    }
  }
}


} // namespace Amanzi


#endif
