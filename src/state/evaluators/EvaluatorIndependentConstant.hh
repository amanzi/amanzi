/*
  State

  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

  A field evaluator with no dependencies, a constant value.
*/

#pragma once

#include "CompositeVectorFunction.hh"
#include "EvaluatorIndependent.hh"
#include "Evaluator_Factory.hh"

namespace Amanzi {

class EvaluatorIndependentConstant
    : public EvaluatorIndependent<CompositeVector, CompositeVectorSpace> {

 public:
  explicit EvaluatorIndependentConstant(Teuchos::ParameterList &plist);
  EvaluatorIndependentConstant(const EvaluatorIndependentConstant &other) = default;
  virtual Teuchos::RCP<Evaluator> Clone() const override;

 protected:
  virtual void Update_(State &S) override;

 protected:
  double value_;

 private:
  static Utils::RegisteredFactory<Evaluator, EvaluatorIndependentConstant> fac_;
};

} // namespace Amanzi


