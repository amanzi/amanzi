/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/*
  Copyright 2010-201x held jointly, see COPYRIGHT.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon
*/

//! Evaluates a function to provide Dirichlet data on faces.


#ifndef AMANZI_EVALUATOR_BCS_DIRICHLET_HH_
#define AMANZI_EVALUATOR_BCS_DIRICHLET_HH_

#include "BCs_Factory.hh"
#include "Evaluator_Factory.hh"
#include "EvaluatorSecondary.hh"
#include "BoundaryFunctionFactory.hh"

namespace Amanzi {

class Evaluator_BCs : public EvaluatorSecondary<Operators::BCs, Operators::BCs_Factory> {
 public:
  explicit
  Evaluator_BCs(Teuchos::ParameterList &plist);

  Evaluator_BCs(const Evaluator_BCs &other) = default;
  virtual Teuchos::RCP<Evaluator> Clone() const override {
    return Teuchos::rcp(new Evaluator_BCs(*this));
  }

  virtual bool
  IsDifferentiableWRT(const State &S, const Key &wrt_key,
                      const Key &wrt_tag) const override final;

  virtual void EnsureCompatibility(State &S) override;
  // virtual void EnsureCompatibleDerivative(State &S,
  //         const Key& wrt_key, const Key& wrt_tag) override;

  virtual void Evaluate_(const State &S, Operators::BCs &result) override;

 protected:
  std::vector<int> bc_types_;

 private:
  static Utils::RegisteredFactory<Evaluator, Evaluator_BCs> fac_;

};


} // namespace

#endif
