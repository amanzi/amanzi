/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Flow PK

  Evaluator for computing the Darcy velocity.
*/

#ifndef AMANZI_FLOW_DARCY_VELOCITY_EVALUATOR_
#define AMANZI_FLOW_DARCY_VELOCITY_EVALUATOR_

#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace Flow {

class DarcyVelocityEvaluator
  : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
 public:
  explicit DarcyVelocityEvaluator(Teuchos::ParameterList& plist);
  DarcyVelocityEvaluator(const DarcyVelocityEvaluator& other);

  // required inteface functions
  virtual Teuchos::RCP<Evaluator> Clone() const override;

 protected:
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override;

  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override{};

  // since cell requires face, the default behavior is not applicable
  virtual void EnsureCompatibility_ToDeps_(State& S) override{
    // if this were implemented, it should call something like:
    // S.Require<CompositeVector,CompositeVectorSpace>(vol_flowrate_key, tag)
    //  .SetMesh(... the mesh ...)
    //  ->AddComponent("face", FACE, 1);
    //
    // Since it is not implemented, you're relying on someone else to make sure
    // that flux has the expected structure.
  };

 protected:
  Key vol_flowrate_key_;
};

} // namespace Flow
} // namespace Amanzi

#endif
