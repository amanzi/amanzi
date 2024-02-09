/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Flow PK

  Evaluator for volumetric flow rate.
*/

#ifndef AMANZI_FLOW_VOLUMETRIC_FLOW_RATE_EVALUATOR_
#define AMANZI_FLOW_VOLUMETRIC_FLOW_RATE_EVALUATOR_

#include "BCs.hh"
#include "Operator.hh"
#include "Upwind.hh"

#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace Flow {

class VolumetricFlowRateEvaluator
  : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
 public:
  VolumetricFlowRateEvaluator(Teuchos::ParameterList& plist, double molar_rho = 0.0);
  VolumetricFlowRateEvaluator(const VolumetricFlowRateEvaluator& other);

  // special initialization function
  void set_upwind(Teuchos::RCP<Operators::Upwind>& upwind) { upwind_ = upwind; }
  void set_bc(Teuchos::RCP<Operators::BCs>& bc) { bc_ = bc; }

  // required inteface functions
  virtual Teuchos::RCP<Evaluator> Clone() const override;

 protected:
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override;

  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override{};

  // since cell-based stress requires nodes, the default behavior is not applicable
  virtual void EnsureCompatibility_ToDeps_(State& S) override{};

 protected:
  Key domain_;
  Key mol_flowrate_key_, mol_density_key_;
  Teuchos::RCP<Operators::Upwind> upwind_;
  Teuchos::RCP<Operators::BCs> bc_;

 private:
  double molar_rho_;
};

} // namespace Flow
} // namespace Amanzi

#endif
