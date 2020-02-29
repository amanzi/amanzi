/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

#include "subgrid_manning_coefficient_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

SubgridManningCoefficientEvaluator::SubgridManningCoefficientEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist) {

  Key domain = Keys::getDomain(my_key_);
  mann_key_ = Keys::readKey(plist_, domain, "manning coefficient", "manning_coefficient");
  dependencies_.insert(mann_key_); 

  frac_cond_key_ = Keys::readKey(plist_, domain, "fractional conductance", "fractional_conductance");
  dependencies_.insert(frac_cond_key_);
  
  beta_key_ = Keys::readKey(plist_, domain, "drag exponent", "drag_exponent");
  dependencies_.insert(beta_key_);
}

void SubgridManningCoefficientEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result) {

  Epetra_MultiVector& res = *result->ViewComponent("cell",false);
  
  const Epetra_MultiVector& mann = *S->GetFieldData(mann_key_)->ViewComponent("cell",false);
  const Epetra_MultiVector& fc = *S->GetFieldData(frac_cond_key_)->ViewComponent("cell",false);
  const Epetra_MultiVector& beta = *S->GetFieldData(beta_key_)->ViewComponent("cell",false);
  
  int ncells = res.MyLength();
  for (int c=0; c!=ncells; ++c) {
    if (fc[0][c] <= 0.) {
      res[0][c] = mann[0][c]; // note this will result in a zero anyway
    } else {
      res[0][c] = mann[0][c] * std::pow(fc[0][c], -(beta[0][c]+1));
    }
  }
}


void SubgridManningCoefficientEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {

  Epetra_MultiVector& res = *result->ViewComponent("cell",false);
  const Epetra_MultiVector& mann = *S->GetFieldData(mann_key_)->ViewComponent("cell",false);
  const Epetra_MultiVector& fc = *S->GetFieldData(frac_cond_key_)->ViewComponent("cell",false);
  const Epetra_MultiVector& beta = *S->GetFieldData(beta_key_)->ViewComponent("cell",false);

  int ncells = res.MyLength();

  if (wrt_key == frac_cond_key_) {
    for (int c=0; c!=ncells; ++c) {
      if (fc[0][c] <= 0.0) {
        res[0][c] = 0.;
      } else {
        res[0][c] = -(beta[0][c] + 1) * mann[0][c] * std::pow(fc[0][c], -(beta[0][c]+2));
      }
    }
  } else {
    Errors::Message msg("SubgridManningCoefficientEvaluator: Not Implemented: no derivatives implemented other than fractional conductance.");
    Exceptions::amanzi_throw(msg);
  }
  
}

} //namespace
} //namespace
} //namespace
