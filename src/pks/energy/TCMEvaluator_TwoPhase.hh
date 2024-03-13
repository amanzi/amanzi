/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  Energy

  Interface for a thermal conductivity model with two phases.
*/

#ifndef AMANZI_ENERGY_TCM_EVALUATOR_TWOPHASE_HH_
#define AMANZI_ENERGY_TCM_EVALUATOR_TWOPHASE_HH_

#include "EvaluatorSecondaryMonotype.hh"
#include "TCM_TwoPhase.hh"

namespace Amanzi {
namespace Energy {

// Equation of State model
class TCMEvaluator_TwoPhase
  : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
 public:
  // constructor format for all derived classes
  TCMEvaluator_TwoPhase(Teuchos::ParameterList& plist);
  TCMEvaluator_TwoPhase(const TCMEvaluator_TwoPhase& other);

  // required inteface functions
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override;

  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override;

 protected:
  Teuchos::RCP<TCM_TwoPhase> tc_;

  // Keys for fields dependencies
  Key porosity_key_, saturation_key_;
  Tag tag_;
};

} // namespace Energy
} // namespace Amanzi

#endif
