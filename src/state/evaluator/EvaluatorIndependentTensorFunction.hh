/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//!

#ifndef AMANZI_INDEPENDENT_TENSOR_FUNCTION_HH_
#define AMANZI_INDEPENDENT_TENSOR_FUNCTION_HH_

#include "CompositeVectorFunction.hh"
#include "EvaluatorIndependent.hh"
#include "Evaluator_Factory.hh"
#include "TensorVector.hh"

namespace Amanzi {

class EvaluatorIndependentTensorFunction
  : public EvaluatorIndependent<TensorVector, TensorVector_Factory> {
 public:
  // ---------------------------------------------------------------------------
  // Constructors
  // ---------------------------------------------------------------------------
  explicit EvaluatorIndependentTensorFunction(Teuchos::ParameterList& plist);
  EvaluatorIndependentTensorFunction(
    const EvaluatorIndependentTensorFunction& other) = default;

  virtual Teuchos::RCP<Evaluator> Clone() const override;

  virtual Evaluator& operator=(const Evaluator& other) override;

  EvaluatorIndependentTensorFunction&
  operator=(const EvaluatorIndependentTensorFunction& other);

  virtual void EnsureCompatibility(State& S) override;

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

 private:
  static Utils::RegisteredFactory<Evaluator, EvaluatorIndependentTensorFunction>
    fac_;
};

} // namespace Amanzi

#endif
