/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//!

#pragma once

#include "EvaluatorIndependent.hh"
#include "Evaluator_Factory.hh"
#include "TensorVector.hh"

namespace Amanzi {

namespace Functions {
class CompositeVectorFunction;
}

namespace Impl {

//
// Hides need parallel_for in a function to avoid issues with public/protected access.
//
template <class View_type>
void
assignViewToTensorVectorDiag(const View_type& v, int j, TensorVector& tv)
{
  Kokkos::parallel_for(
    "TensorFunctionUpdate", v.extent(0), KOKKOS_LAMBDA(const int& i) {
      auto Ti = tv.at(j + i);
      Ti(0, 0) = v(i, 0);
    });
}
} // namespace Impl

class EvaluatorIndependentTensorFunction
  : public EvaluatorIndependent<TensorVector, TensorVector_Factory> {
 public:
  // ---------------------------------------------------------------------------
  // Constructors
  // ---------------------------------------------------------------------------
  explicit EvaluatorIndependentTensorFunction(const Teuchos::RCP<Teuchos::ParameterList>& plist);
  EvaluatorIndependentTensorFunction(const EvaluatorIndependentTensorFunction& other) = default;

  virtual Teuchos::RCP<Evaluator> Clone() const override;

  virtual void EnsureCompatibility(State& S) override;

  static const std::string eval_type;
  virtual std::string getType() const override { return eval_type; }

 protected:
  // ---------------------------------------------------------------------------
  // Update the value in the state.
  // ---------------------------------------------------------------------------
  virtual void Update_(State& S) override;

 protected:
  Teuchos::RCP<Functions::CompositeVectorFunction> func_;
  int num_funcs_;
  int dimension_;
  int rank_;
  double rescaling_;

 private:
  static Utils::RegisteredFactory<Evaluator, EvaluatorIndependentTensorFunction> reg_;
};

} // namespace Amanzi