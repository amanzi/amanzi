/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  Interface for a thermal conductivity model with three phases.

  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_ENERGY_RELATIONS_TC_THREEPHASE_EVALUATOR_HH_
#define AMANZI_ENERGY_RELATIONS_TC_THREEPHASE_EVALUATOR_HH_

#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace Energy {
namespace EnergyRelations {

// Equation of State model
class ThermalConductivityThreePhaseEvaluator :
    public SecondaryVariableFieldEvaluator {

 public:
  // constructor format for all derived classes
  ThermalConductivityThreePhaseEvaluator(Teuchos::ParameterList& tc_plist);
  ThermalConductivityThreePhaseEvaluator(const ThermalConductivityThreePhaseEvaluator& other);

  Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldModel
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

 protected:

  // PList
  Teuchos::ParameterList tc_plist_;

  Teuchos::RCP<ThermalConductivityThreePhase> tc_;

  // Keys for fields
  // dependencies
  Key poro_key_;
  Key sat_key_;
  Key sat2_key_;
};

} // namespace
} // namespace
} // namespace

#endif
