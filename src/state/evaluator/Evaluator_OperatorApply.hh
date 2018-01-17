/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/*
  Copyright 2010-201x held jointly, see COPYRIGHT.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon
*/

//! Evaluator_OperatorApply is a secondary evaluator that calculates r = b - Ax

/*!


*/

#ifndef STATE_EVALUATOR_OPERATOR_APPLY_HH_
#define STATE_EVALUATOR_OPERATOR_APPLY_HH_

#include "EvaluatorSecondary.hh"
#include "Evaluator_Factory.hh"

namespace Amanzi {

class Evaluator_OperatorApply
    : public EvaluatorSecondary<CompositeVector, CompositeVectorSpace> {
public:
  explicit Evaluator_OperatorApply(Teuchos::ParameterList &plist);
  Evaluator_OperatorApply(const Evaluator_OperatorApply &other) = default;

  virtual Teuchos::RCP<Evaluator> Clone() const override {
    return Teuchos::rcp(new Evaluator_OperatorApply(*this));
  }

  void EnsureCompatibility(State &S) override;

protected:
  // These do the actual work
  virtual void Evaluate_(const State &S, CompositeVector &result) override;

  // This should get some careful thought of the right strategy.  Punting for
  // now --etc
  virtual void EvaluatePartialDerivative_(const State &S, const Key &wrt_key,
                                          const Key &wrt_tag,
                                          CompositeVector &result) override;

protected:
  Key rhs_key_;
  Key operator_key_;
  Key x_key_;
  std::vector<Key> local_op_keys_;

  bool is_square_;

private:
  static Utils::RegisteredFactory<Evaluator, Evaluator_OperatorApply> fac_;
};

} // namespace Amanzi

#endif
