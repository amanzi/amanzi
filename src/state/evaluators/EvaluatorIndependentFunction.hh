/*
  State

  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

  A field evaluator with no dependencies specified by a function.
*/

#ifndef AMANZI_STATE_INDEPENDENT_EVALUATOR_FUNCTION_
#define AMANZI_STATE_INDEPENDENT_EVALUATOR_FUNCTION_

#include "CompositeVectorFunction.hh"
#include "EvaluatorIndependent.hh"
#include "Evaluator_Factory.hh"

namespace Amanzi {

class EvaluatorIndependentFunction
    : public EvaluatorIndependent<CompositeVector, CompositeVectorSpace> {

 public:
  explicit EvaluatorIndependentFunction(Teuchos::ParameterList& plist);
  EvaluatorIndependentFunction(const EvaluatorIndependentFunction& other) = default;

  virtual Teuchos::RCP<Evaluator> Clone() const override;

  virtual Evaluator& operator=(const Evaluator& other) override;

  EvaluatorIndependentFunction& operator=(const EvaluatorIndependentFunction& other);

 protected:
  virtual void Update_(State& S) override;

 protected:
  Teuchos::RCP<Functions::CompositeVectorFunction> func_;

 private:
  static Utils::RegisteredFactory<Evaluator, EvaluatorIndependentFunction> fac_;
};

} // namespace Amanzi

#endif
