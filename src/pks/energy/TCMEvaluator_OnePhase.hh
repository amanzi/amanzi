/*
  Energy

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Thermal conductivity evaluator with one liquid (water) phase.
*/

#ifndef AMANZI_ENERGY_TCM_EVALUATOR_ONEPHASE_HH_
#define AMANZI_ENERGY_TCM_EVALUATOR_ONEPHASE_HH_

#include "secondary_variable_field_evaluator.hh"
#include "ThermalConductivity_Water.hh"

namespace Amanzi {
namespace Energy {

// Equation of State model
class TCMEvaluator_OnePhase : public SecondaryVariableFieldEvaluator {
 public:
  // constructor format for all derived classes
  TCMEvaluator_OnePhase(Teuchos::ParameterList& plist);
  TCMEvaluator_OnePhase(const TCMEvaluator_OnePhase& other);

  Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldModel
  virtual void EvaluateField_(
      const Teuchos::Ptr<State>& S, const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(
      const Teuchos::Ptr<State>& S, Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

 protected:
  // We have one model so far; hence, no factory is needed.
  Teuchos::RCP<AmanziEOS::ThermalConductivity_Water> tc_;
  double k_rock_;

  // Keys for fields dependencies
  Key temperature_key_, porosity_key_;
};

}  // namespace Energy
}  // namespace Amanzi

#endif
