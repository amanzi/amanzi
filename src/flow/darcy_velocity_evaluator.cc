/*
  Evaluator for darcy_velocity(darcy_flux)

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include "darcy_velocity_evaluator.hh"

namespace Amanzi {
namespace AmanziFlow {

/* ******************************************************************
* A simple constructor: create dependencies.
****************************************************************** */
DarcyVelocityEvaluator::DarcyVelocityEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist) {
  // hard-coded keys
  my_key_ = "darcy_velocity";
  darcy_flux_key_ = "darcy_flux";
  dependencies_.insert(darcy_flux_key_);
}


/* ******************************************************************
* Clone with unclear yet purpose.
****************************************************************** */
Teuchos::RCP<FieldEvaluator>
DarcyVelocityEvaluator::Clone() const {
  return Teuchos::rcp(new DarcyVelocityEvaluator(*this));
}


/* ******************************************************************
* Required member function: basic algorithm.
****************************************************************** */
void DarcyVelocityEvaluator::EvaluateField_(
    const Teuchos::Ptr<State>& S, const Teuchos::Ptr<CompositeVector>& result) {

  const CompositeVector& flux = *S->GetFieldData(darcy_flux_key_);
  Epetra_MultiVector& result_c = *(result->ViewComponent("cell", false));
}

}  // namespace AmanziFlow
}  // namespace Amanzi
