/*
  Energy

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)

  Interface for a thermal conductivity model with two phases.
*/

#ifndef AMANZI_ENERGY_TCM_EVALUATOR_TWOPHASE_HH_
#define AMANZI_ENERGY_TCM_EVALUATOR_TWOPHASE_HH_

#include "secondary_variable_field_evaluator.hh"
#include "TCM_TwoPhase.hh"

namespace Amanzi {
namespace Energy {

// Equation of State model
class TCMEvaluator_TwoPhase : public SecondaryVariableFieldEvaluator {
 public:
  // constructor format for all derived classes
  TCMEvaluator_TwoPhase(Teuchos::ParameterList& plist);
  TCMEvaluator_TwoPhase(const TCMEvaluator_TwoPhase& other);

  Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldModel
  virtual void EvaluateField_(
      const Teuchos::Ptr<State>& S, const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(
      const Teuchos::Ptr<State>& S, Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

 protected:
  Teuchos::RCP<TCM_TwoPhase> tc_;

  // Keys for fields dependencies
  Key porosity_key_;
  Key saturation_key_;
};

}  // namespace Energy
}  // namespace Amanzi

#endif
