/*
  This is the energy component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)

  Interface for a thermal conductivity model with two phases.
*/

#ifndef AMANZI_ENERGY_RELATIONS_TC_TWOPHASE_EVALUATOR_HH_
#define AMANZI_ENERGY_RELATIONS_TC_TWOPHASE_EVALUATOR_HH_

#include "secondary_variable_field_evaluator.hh"
#include "twophase_thermal_conductivity.hh"

namespace Amanzi {
namespace Energy {

// Equation of State model
class ThermalConductivityTwoPhaseEvaluator : public SecondaryVariableFieldEvaluator {
 public:
  // constructor format for all derived classes
  ThermalConductivityTwoPhaseEvaluator(Teuchos::ParameterList& plist);
  ThermalConductivityTwoPhaseEvaluator(const ThermalConductivityTwoPhaseEvaluator& other);

  Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldModel
  virtual void EvaluateField_(
      const Teuchos::Ptr<State>& S, const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(
      const Teuchos::Ptr<State>& S, Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

 protected:
  Teuchos::RCP<ThermalConductivityTwoPhase> tc_;

  // Keys for fields dependencies
  Key poro_key_;
  Key sat_key_;
};

}  // namespace Energy
}  // namespace Amanzi

#endif
