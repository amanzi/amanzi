/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon (coonet@ornl.gov)
*/

//! Implements accumulation terms.

/*!

*/

#pragma once
#include "Teuchos_ParameterList.hpp"

#include "Evaluator_Factory.hh"
#include "State.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {

class Evaluator_PDE_Accumulation
  : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
 public:
  Evaluator_PDE_Accumulation(Teuchos::ParameterList& plist);

  Evaluator_PDE_Accumulation(const Evaluator_PDE_Accumulation& other) = default;
  virtual Teuchos::RCP<Evaluator> Clone() const override
  {
    return Teuchos::rcp(new Evaluator_PDE_Accumulation(*this));
  };

  virtual void EnsureCompatibility(State& S) override;

  virtual std::string name() const override { return "accumulation operator"; }
  
 protected:
  virtual void EvaluatePartialDerivative_(
    const State& S, const Key& wrt_key, const Key& wrt_tag,
    const std::vector<CompositeVector*>& results) override;

  virtual void Evaluate_(const State& S,
                         const std::vector<CompositeVector*>& results) override;

 protected:
  Key conserved_key_;
  Key cv_key_;
  Key tag_old_, tag_new_;

 private:
  static Utils::RegisteredFactory<Evaluator, Evaluator_PDE_Accumulation> fac_;
};

} // namespace Amanzi
