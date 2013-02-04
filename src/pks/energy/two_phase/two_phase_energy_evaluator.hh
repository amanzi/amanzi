/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -----------------------------------------------------------------------------
ATS

Authors: Ethan Coon (ecoon@lanl.gov)

FieldEvaluator for internal energy.

Wrapping this conserved quantity as a field evaluator makes it easier to take
derivatives, keep updated, and the like.  The equation for this is simply:

IE = phi * (s_liquid * n_liquid * u_liquid + s_gas * n_gas * u_gas )
  + (1 - phi) * rho_rock * u_rock

This is simply the conserved quantity in the energy equation.
----------------------------------------------------------------------------- */


#ifndef AMANZI_TWO_PHASE_ENERGY_EVALUATOR_HH_
#define AMANZI_TWO_PHASE_ENERGY_EVALUATOR_HH_

#include "Teuchos_ParameterList.hpp"

#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace Energy {

class TwoPhaseEnergyEvaluator : public SecondaryVariableFieldEvaluator {

public:
  explicit
  TwoPhaseEnergyEvaluator(Teuchos::ParameterList& energy_plist);
  TwoPhaseEnergyEvaluator(const TwoPhaseEnergyEvaluator& other);

  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

};

} // namespace
} // namespace

#endif
