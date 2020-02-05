/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//! A launcher for constructing Derivative tags based upon a WRT key/value tag.

/*!

  This helps implement a very generic EvaluatorSecondaryMonotype which uses a
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

  Use of this launcher is tested in
  src/state/test/state_evaluators_dag_models.cc.

*/

#ifndef STATE_EVALUATOR_MODEL_LAUNCHER_HH_
#define STATE_EVALUATOR_MODEL_LAUNCHER_HH_

#include "AmanziTypes.hh"
#include "StateDefs.hh"

namespace Amanzi {
namespace Impl {


//
// A Launcher class that launches a parallel_for to calculate the partial
// derivative kernel if I is the index of wrt in the model's dependencies, or
// recurses to I-1 if not.
template <int I, class Model_type, class Device_type>
class EvaluatorModelLauncher {
 public:
  EvaluatorModelLauncher(const std::string& name, const KeyPair& wrt,
                         const KeyPairVector& dependencies,
                         const Model_type& model)
    : name_(name), wrt_(wrt), dependencies_(dependencies), model_(model)
  {}

  void launch(const int extent)
  {
    if (wrt_ == dependencies_[I]) {
      Kokkos::RangePolicy<Deriv<I>, typename Device_type::execution_space>
        range(0, extent);
      Kokkos::parallel_for(
        "EvaluatorModelLauncher::launch",
        name_, range, model_);
    } else {
      EvaluatorModelLauncher<I - 1, Model_type, Device_type> launcher(
        name_, wrt_, dependencies_, model_);
      launcher.launch(extent);
    }
  }

  void launch(const AmanziMesh::Entity_ID_View& material_ids)
  {
    if (wrt_ == dependencies_[I]) {
      Kokkos::RangePolicy<typename Device_type::execution_space> range(
        0, material_ids.extent(0));
      Kokkos::parallel_for(
        "EvaluatorModelLauncher::launch",
        name_, range, KOKKOS_LAMBDA(const int i) {
        model_(Deriv<I>(), material_ids(i));
      });
    } else {
      EvaluatorModelLauncher<I - 1, Model_type, Device_type> launcher(
        name_, wrt_, dependencies_, model_);
      launcher.launch(material_ids);
    }
  }

 private:
  const std::string& name_;
  const KeyPair& wrt_;
  const KeyPairVector& dependencies_;
  const Model_type& model_;
};

//
// Partial specialization of this class to terminate the recursion.
template <class Model_type, class Device_type>
class EvaluatorModelLauncher<-1, Model_type, Device_type> {
 public:
  EvaluatorModelLauncher(const std::string& name, const KeyPair& wrt,
                         const KeyPairVector& dependencies,
                         const Model_type& model)
    : name_(name), wrt_(wrt)
  {}

  void launch(const int extent)
  {
    Errors::Message msg;
    msg << "EvaluatorModel (" << name_
        << "): requested derivative with respect to ( \"" << wrt_.first
        << "\",\"" << wrt_.second << "\" ) is not provided.";
    throw(msg);
  }

  void launch(const AmanziMesh::Entity_ID_View& material_ids)
  {
    Errors::Message msg;
    msg << "EvaluatorModel (" << name_
        << "): requested derivative with respect to ( \"" << wrt_.first
        << "\",\"" << wrt_.second << "\" ) is not provided.";
    throw(msg);
  }

 private:
  const std::string& name_;
  const KeyPair& wrt_;
};

} // namespace Impl
} // namespace Amanzi


#endif
