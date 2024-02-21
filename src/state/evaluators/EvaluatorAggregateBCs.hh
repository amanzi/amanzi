/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//! Evaluates a function to provide Dirichlet data on faces.

#pragma once

#include "Evaluator_Factory.hh"
#include "EvaluatorSecondary.hh"

namespace Amanzi {

class EvaluatorAggregateBCs : public EvaluatorSecondary {
 public:
  explicit EvaluatorAggregateBCs(const Teuchos::RCP<Teuchos::ParameterList>& plist);

  EvaluatorAggregateBCs(const EvaluatorAggregateBCs& other) = default;
  virtual Teuchos::RCP<Evaluator> Clone() const override
  {
    return Teuchos::rcp(new EvaluatorAggregateBCs(*this));
  }

  virtual bool
  IsDifferentiableWRT(const State& S, const Key& wrt_key, const Tag& wrt_tag) const override
  {
    return false;
  }

  virtual void EnsureCompatibility(State& S) override;

  virtual void Update_(State& S) override;

  virtual void UpdateDerivative_(State& S, const Key& wrt_key, const Tag& wrt_tag) override final
  {
    AMANZI_ASSERT(false); // never called
  }

  static const std::string eval_type;
  virtual std::string getType() const override { return eval_type; }

 protected:
  bool inited_;

 private:
  static Utils::RegisteredFactory<Evaluator, EvaluatorAggregateBCs> fac_;
};


} // namespace Amanzi
