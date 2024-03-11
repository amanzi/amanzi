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
  user-provided Model.  That model provides the list of keys calculated,
  dependencies, and operator() methods that implement the model (and
  derivatives) to act as a Kokkos functor.

  Currently, this class assumes that tags of the keys calculated and
  dependencies are all the same.

  Models are expected to implement a fairly precise set of variables and
  methods.  Examples of these are shown (for a simple synthetic case) in
  `src/state/test/dag_models.hh`.

  Use of this evaluator is tested in
  src/state/test/state_evaluators_dag_models.cc.

  Developer note:

  The implementation of this is a bit tricky, mainly
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

*/


#ifndef STATE_EVALUATOR_MODEL_HH_
#define STATE_EVALUATOR_MODEL_HH_

#include "AmanziTypes.hh"
#include "CompositeVectorSpace.hh"
#include "CompositeVector.hh"
#include "Factory.hh"
#include "StateDefs.hh"
#include "EvaluatorModelLauncher.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {


//
// The actual generic Evaluator class.
//
template <template <class, class> class Model, class Device_type = DefaultDevice>
class EvaluatorModelCV : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
 public:
  using View_type = MultiVectorView_type<Device_type>;
  using cView_type = cMultiVectorView_type<Device_type>;
  using Model_type = Model<cView_type, View_type>;

  EvaluatorModelCV(const Teuchos::RCP<Teuchos::ParameterList>& plist);

  virtual Teuchos::RCP<Evaluator> Clone() const override;
  virtual std::string getType() const override { return Model_type::eval_type; }

 protected:
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override;

  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override;

  virtual void EnsureCompatibility_Structure_(State& S) override
  {
    // make sure all my_keys have the same CVS
    EnsureCompatibility_StructureSame_(S);
  }

 protected:
  Teuchos::RCP<Model_type> model_;
  std::string name_;

 private:
  // registration in the evaluator factory
  static Utils::RegisteredFactory<Evaluator, EvaluatorModelCV<Model, Device_type>> reg_;
};


//
// Implementation
// -----------------------------------------------------------------------------
template <template <class, class> class Model, class Device_type>
EvaluatorModelCV<Model, Device_type>::EvaluatorModelCV(
  const Teuchos::RCP<Teuchos::ParameterList>& plist)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(plist),
    model_(Teuchos::rcp(new Model_type(plist))),
    name_(Keys::cleanPListName(plist->name()))
{
  my_keys_.clear();
  my_keys_ = model_->getMyKeys();
  for (const KeyTag& dep : model_->getDependencies()) dependencies_.insert(dep);
}

template <template <class, class> class Model, class Device_type>
Teuchos::RCP<Evaluator>
EvaluatorModelCV<Model, Device_type>::Clone() const
{
  return Teuchos::rcp(new EvaluatorModelCV<Model, Device_type>(*this));
}

template <template <class, class> class Model, class Device_type>
void
EvaluatorModelCV<Model, Device_type>::Evaluate_(const State& S,
                                                const std::vector<CompositeVector*>& results)
{
  AMANZI_ASSERT(results.size() > 0);

  for (const auto& comp : *results[0]) {
    std::vector<cView_type> dependency_views;
    for (const auto& dep : dependencies_) {
      const auto& vec = S.Get<CompositeVector>(dep.first, dep.second);
      dependency_views.emplace_back(vec.viewComponent(comp, false));
    }

    std::vector<View_type> result_views;
    for (auto result : results) { result_views.emplace_back(result->viewComponent(comp, false)); }

    // set up the model and range and then dispatch
    model_->setViews(dependency_views, result_views, S);
    Kokkos::RangePolicy<typename Device_type::execution_space> range(0, result_views[0].extent(0));
    Kokkos::parallel_for(name_, range, *model_);
    Kokkos::fence();

    // Reset views
    // model_ = Teuchos::rcp<Model_type>(new Model_type(plist_));
  }
  //Debug_(S);
}

template <template <class, class> class Model, class Device_type>
void
EvaluatorModelCV<Model, Device_type>::EvaluatePartialDerivative_(
  const State& S,
  const Key& wrt_key,
  const Tag& wrt_tag,
  const std::vector<CompositeVector*>& results)
{
  AMANZI_ASSERT(results.size() > 0);
  for (const auto& comp : *results[0]) {
    std::vector<cView_type> dependency_views;
    for (const auto& dep : dependencies_) {
      const auto& vec = S.Get<CompositeVector>(dep.first, dep.second);
      dependency_views.emplace_back(vec.viewComponent(comp, false));
    }

    std::vector<View_type> result_views;
    for (auto result : results) { result_views.emplace_back(result->viewComponent(comp, false)); }

    // set up the model and range and then dispatch
    model_->setViews(dependency_views, result_views, S);

    KeyTag wrt{wrt_key, wrt_tag};
    Impl::EvaluatorModelLauncher<Model_type::n_dependencies - 1, Model_type, Device_type> launcher(
      name_, wrt, dependencies_, *model_);
    launcher.launch(result_views[0].extent(0));
    Kokkos::fence();

    // Reset views
    //model_ = Teuchos::rcp<Model_type>(new Model_type(plist_));
  }
}

} // namespace Amanzi


#endif
