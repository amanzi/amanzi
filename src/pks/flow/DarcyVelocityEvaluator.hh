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

// #include "factory.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace Flow {

class DarcyVelocityEvaluator : public SecondaryVariableFieldEvaluator {
 public:
  explicit
  DarcyVelocityEvaluator(Teuchos::ParameterList& plist);
  DarcyVelocityEvaluator(const DarcyVelocityEvaluator& other);

  virtual void EnsureCompatibility(const Teuchos::Ptr<State>& S) {};
  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

 protected:
  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {};  // should not be called

 protected:
  Key darcy_flux_key_;
};

}  // namespace Flow
}  // namespace Amanzi

#endif
