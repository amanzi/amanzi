/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*!

Computes a reconstructed velocity field given a face-based flux field and
density.

`"evaluator type`" = `"velocity reconstruction`"

.. _evaluator-velocity-reconstruction-spec:
.. admonition:: evaluator-velocity-reconstruction-spec

   KEYS:
   - `"water flux`"
   - `"molar density liquid`"

*/

#pragma once

#include "EvaluatorSecondaryMonotype.hh"
#include "Evaluator_Factory.hh"

namespace Amanzi {
class State;

class EvaluatorVelocityReconstruction
  : public EvaluatorSecondaryMonotype<CompositeVector,CompositeVectorSpace> {

public:
  explicit EvaluatorVelocityReconstruction(Teuchos::ParameterList& plist);
  EvaluatorVelocityReconstruction(const EvaluatorVelocityReconstruction& other) = default;
  virtual Teuchos::RCP<Evaluator> Clone() const override;

 protected:
  // These do the actual work
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override;

  // This should get some careful thought of the right strategy.  Punting for
  // now --etc
  virtual bool IsDifferentiableWRT(const State& S, const Key& wrt_key, const Tag& wrt_tag) const override {
    return false;
  }

  virtual void UpdateDerivative_(State& S, const Key& wrt_key, const Tag& wrt_tag) override {}

  // This should get some careful thought of the right strategy.  Punting for
  // now --etc
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override
  {}

protected:
  // helper function -- calls Require on all dependencies
  // Called in Basics
  virtual void EnsureCompatibility_ToDeps_(State& S) override;

protected:
  Key flux_key_;
  Key molar_dens_key_;
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;

 private:
  static Utils::RegisteredFactory<Evaluator, EvaluatorVelocityReconstruction> fac_;
};

} // namespace Amanzi


