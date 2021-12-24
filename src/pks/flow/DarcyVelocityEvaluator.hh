/*
  Flow PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov) 
           Konstantin Lipnikov (lipnikov@lanl.gov)

  Evaluator for determining darcy_velocity(darcy_flux).
*/

#ifndef AMANZI_FLOW_DARCY_VELOCITY_EVALUATOR_
#define AMANZI_FLOW_DARCY_VELOCITY_EVALUATOR_

// #include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace Flow {

class DarcyVelocityEvaluator : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
 public:
  explicit DarcyVelocityEvaluator(Teuchos::ParameterList& plist);
  DarcyVelocityEvaluator(const DarcyVelocityEvaluator& other);

  // required inteface functions
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override;

  virtual void EvaluatePartialDerivative_(const State& S, const Key& wrt_key, const Key& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override {};

  // since cell reequires face, the default behavior is not applicable
  virtual void EnsureCompatibility(State& S) override {};

 protected:
  Key darcy_flux_key_;
};

}  // namespace Flow
}  // namespace Amanzi

#endif
