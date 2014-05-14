/*
  Evaluator for determining darcy_velocity(darcy_flux).

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_FLOW_DARCY_VELOCITY_EVALUATOR_
#define AMANZI_FLOW_DARCY_VELOCITY_EVALUATOR_

// #include "factory.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace AmanziFlow {

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

}  // namespace AmanziFlow
}  // namespace Amanzi

#endif
