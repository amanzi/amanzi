/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Energy

  Thermal conductivity evaluator with one liquid (water) phase.
*/

#ifndef AMANZI_ENERGY_TCM_EVALUATOR_ONEPHASE_HH_
#define AMANZI_ENERGY_TCM_EVALUATOR_ONEPHASE_HH_

#include "EvaluatorSecondaryMonotype.hh"
#include "H2O_ThermalConductivity.hh"

namespace Amanzi {
namespace Energy {

// Equation of State model
class TCMEvaluator_OnePhase
  : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
 public:
  // constructor format for all derived classes
  TCMEvaluator_OnePhase(Teuchos::ParameterList& plist);
  TCMEvaluator_OnePhase(const TCMEvaluator_OnePhase& other);

  // required inteface functions
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override;

  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override;

 protected:
  // We have one model so far; hence, no factory is needed.
  Teuchos::RCP<AmanziEOS::H2O_ThermalConductivity> tc_;
  double k_rock_;

  // Keys for fields dependencies
  Key temperature_key_, porosity_key_;
};

} // namespace Energy
} // namespace Amanzi

#endif
