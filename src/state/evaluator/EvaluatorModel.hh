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


#ifndef STATE_EVALUATOR_MODEL_HH_
#define STATE_EVALUATOR_MODEL_HH_

#include "AmanziTypes.hh"
#include "StateDefs.hh"

namespace Amanzi {

//
// A Launcher class that launches a parallel_for to calculate the partial
// derivative kernel if I is the index of wrt in the model's dependencies, or
// recurses to I-1 if not.
template<int I, class Model_type, class Device_type>
class Launcher {
 public:
  Launcher(const std::string& name,
           const KeyPair& wrt,
           const KeyPairVector& dependencies,
           const Model_type& model,
           const int extent) :
      name_(name),
      wrt_(wrt),
      dependencies_(dependencies),
      model_(model),
      extent_(extent) {}

  void launch() {
    // std::cout << "wrt: " << wrt_.first << " matches " << dependencies_[I].first << " in slot " << I << "?";
    if (wrt_ == dependencies_[I]) {
      // std::cout << "  YUP, go!" << std::endl;
      Kokkos::RangePolicy<Deriv<I>, typename Device_type::execution_space> range(0, extent_);
      Kokkos::parallel_for(name_, range, model_);
    } else {
      // std::cout << "  nope, recurse!" << std::endl;
      Launcher<I-1,Model_type,Device_type> launcher(name_, wrt_, dependencies_, model_, extent_);
      launcher.launch();
    }              
  }

 private:
  const std::string& name_;
  const KeyPair& wrt_;
  const KeyPairVector& dependencies_;
  const Model_type& model_;
  const int extent_;  
};

//
// Partial specialization of this class to terminate the recursion.
template<class Model_type, class Device_type>
class Launcher<-1, Model_type, Device_type> {
 public:
  Launcher(const std::string& name,
           const KeyPair& wrt,
           const KeyPairVector& dependencies,
           const Model_type& model,
           const int extent)
      : name_(name),
        wrt_(wrt) {}

  void launch() {
    Errors::Message msg;
    msg << "EvaluatorModel (" << name_ << "): requested derivative with respect to ( \""
        << wrt_.first << "\",\"" << wrt_.second << "\" ) is not provided.";
    throw(msg);
  }

 private:
  const std::string& name_;
  const KeyPair& wrt_;
  
};


//
// The actual generic Evaluator class.
//
template<template<class,class> class Model, class Device_type=AmanziDefaultDevice>
class EvaluatorModel_CompositeVector : public EvaluatorSecondaryMonotype<CompositeVector,CompositeVectorSpace> {
 public:
  using Model_type = Model<cVectorView_type<Device_type>, VectorView_type<Device_type>>;
  using View_type = VectorView_type<Device_type>;
  using cView_type = cVectorView_type<Device_type>;
  
  EvaluatorModel_CompositeVector(Teuchos::ParameterList& plist)
      : EvaluatorSecondaryMonotype<CompositeVector,CompositeVectorSpace>(plist),
        model_(Teuchos::rcp(new Model_type(plist))),
        name_(plist.name())
  {
    dependencies_ = model_->dependencies();
  }

  virtual Teuchos::RCP<Evaluator> Clone() const override
  {
    return Teuchos::rcp(new EvaluatorModel_CompositeVector<Model,Device_type>(*this));
  }

  virtual void
  Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override
  {
    AMANZI_ASSERT(results.size() > 0);
    for (const auto& comp : *results[0]) {
      std::vector<cView_type> dependency_views;
      for (const auto& dep : dependencies_) {
        const auto& vec = S.Get<CompositeVector>(dep.first, dep.second);
        if (vec.getNumVectors(comp) != 1) {
          Errors::Message msg("EvaluatorModel: expects only dependencies with one DoF.");
          throw(msg);
        }        
        dependency_views.emplace_back(vec.ViewComponent(comp, 0, false));
      }

      std::vector<View_type> result_views;
      for (auto result : results) {
        if (result->getNumVectors(comp) != 1) {
          Errors::Message msg("EvaluatorModel: expects only results with one DoF.");
          throw(msg);
        }        
        result_views.emplace_back(result->ViewComponent(comp, 0, false));
      }

      // set up the model and range and then dispatch
      model_->SetViews(dependency_views, result_views);
      Kokkos::RangePolicy<typename Device_type::execution_space> range(0, result_views[0].extent(0));
      Kokkos::parallel_for(name_, range, *model_);
    }
  }

  virtual void
  EvaluatePartialDerivative_(const State& S, const Key& wrt_key,
                             const Key& wrt_tag,
                             const std::vector<CompositeVector*>& results) override {
    AMANZI_ASSERT(results.size() > 0);
    for (const auto& comp : *results[0]) {
      std::vector<cView_type> dependency_views;
      for (const auto& dep : dependencies_) {
        const auto& vec = S.Get<CompositeVector>(dep.first, dep.second);
        if (vec.getNumVectors(comp) != 1) {
          Errors::Message msg("EvaluatorModel: expects only dependencies with one DoF.");
          throw(msg);
        }        
        dependency_views.emplace_back(vec.ViewComponent(comp, 0, false));
      }

      std::vector<View_type> result_views;
      for (auto result : results) {
        if (result->getNumVectors(comp) != 1) {
          Errors::Message msg("EvaluatorModel: expects only results with one DoF.");
          throw(msg);
        }        
        result_views.emplace_back(result->ViewComponent(comp, 0, false));
      }

      // set up the model and range and then dispatch
      model_->SetViews(dependency_views, result_views);
      auto wrt = std::make_pair(wrt_key, wrt_tag);

      Launcher<Model_type::n_dependencies-1, Model_type, Device_type> launcher(name_, wrt, dependencies_, *model_, result_views[0].extent(0));
      launcher.launch();
    }
  }

  Teuchos::RCP<Model_type> model_;
  std::string name_;
};

} // namespace Amanzi


#endif
