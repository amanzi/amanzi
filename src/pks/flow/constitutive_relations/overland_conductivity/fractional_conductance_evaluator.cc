/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

#include "fractional_conductance_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

FractionalConductanceEvaluator::FractionalConductanceEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist) {

  Key domain = Keys::getDomain(my_key_);

  vpd_key_ = Keys::readKey(plist_, domain, "volumetric ponded depth", "volumetric_ponded_depth");
  dependencies_.insert(vpd_key_); 

  pdd_key_ = Keys::readKey(plist_, domain, "ponded depth minus depression depth", "ponded_depth_minus_depression_depth");
  dependencies_.insert(pdd_key_);
  
  delta_max_key_ = Keys::readKey(plist_, domain, "microtopographic relief", "microtopographic_relief");
  dependencies_.insert(delta_max_key_);

  delta_ex_key_ = plist_.get<std::string>("excluded volume key", Keys::getKey(domain,"excluded_volume"));
  dependencies_.insert(delta_ex_key_);

  depr_depth_key_ = plist_.get<std::string>("depression depth key", Keys::getKey(domain,"depression_depth"));
  dependencies_.insert(depr_depth_key_);
}


FractionalConductanceEvaluator::FractionalConductanceEvaluator(const FractionalConductanceEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    pdd_key_(other.pdd_key_),
    vpd_key_(other.vpd_key_),
    delta_ex_key_(other.delta_ex_key_),
    delta_max_key_(other.delta_max_key_),
    depr_depth_key_(other.depr_depth_key_)
{};

Teuchos::RCP<FieldEvaluator>
FractionalConductanceEvaluator::Clone() const {
  return Teuchos::rcp(new FractionalConductanceEvaluator(*this));
}


void FractionalConductanceEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result) {

  Epetra_MultiVector& res = *result->ViewComponent("cell",false);
  
  const Epetra_MultiVector& vpd = *S->GetFieldData(vpd_key_)->ViewComponent("cell",false);
  const Epetra_MultiVector& delta_max_v = *S->GetFieldData(delta_max_key_)->ViewComponent("cell", false);
  const Epetra_MultiVector& delta_ex_v = *S->GetFieldData(delta_ex_key_)->ViewComponent("cell", false);
  const Epetra_MultiVector& depr_depth_v = *S->GetFieldData(depr_depth_key_)->ViewComponent("cell", false);
  const Epetra_MultiVector& pdd_v = *S->GetFieldData(pdd_key_)->ViewComponent("cell",false);
  
  int ncells = res.MyLength();

  for (int c=0; c!=ncells; ++c) {
    double depr_depth = depr_depth_v[0][c];
    double delta_max = delta_max_v[0][c];
    double delta_ex = delta_ex_v[0][c];
    double fixed_depth = std::pow(depr_depth,2)*(2*delta_max - 3*delta_ex)/std::pow(delta_max,2) + std::pow(depr_depth,3)*(2*delta_ex - delta_max)/std::pow(delta_max,3);
    
    if (pdd_v[0][c] <= 0.0)
      res[0][c] = 0;
    else{
      res[0][c] = (vpd[0][c] - fixed_depth) / (pdd_v[0][c]);
    }
  }
}


void FractionalConductanceEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {

  Epetra_MultiVector& res = *result->ViewComponent("cell",false);
  
  const Epetra_MultiVector& vpd = *S->GetFieldData(vpd_key_)->ViewComponent("cell",false);
  const Epetra_MultiVector& delta_max_v = *S->GetFieldData(delta_max_key_)->ViewComponent("cell", false);
  const Epetra_MultiVector& delta_ex_v = *S->GetFieldData(delta_ex_key_)->ViewComponent("cell", false);
  const Epetra_MultiVector& depr_depth_v = *S->GetFieldData(depr_depth_key_)->ViewComponent("cell", false);
  const Epetra_MultiVector& pdd_v = *S->GetFieldData(pdd_key_)->ViewComponent("cell",false);
  
  int ncells = res.MyLength();
  //  assert(depr_depth_v[0][0] > 0.);

  if (wrt_key == pdd_key_) {
    for (int c=0; c!=ncells; ++c) {
      double depr_depth = depr_depth_v[0][c];
      double delta_max = delta_max_v[0][c];
      double delta_ex = delta_ex_v[0][c];
      double fixed_depth = std::pow(depr_depth,2)*(2*delta_max - 3*delta_ex)/std::pow(delta_max,2) + std::pow(depr_depth,3)*(2*delta_ex - delta_max)/std::pow(delta_max,3);
      
      if (pdd_v[0][c] <= 0.0)
        res[0][c] = 0;
      else{
        res[0][c] = - (vpd[0][c] - fixed_depth) / std::pow(pdd_v[0][c],2.);
      }
    }
    
  }
  else if (wrt_key == vpd_key_) {
    for (int c=0; c!=ncells; ++c) {
      double depr_depth = depr_depth_v[0][c];
      double delta_max = delta_max_v[0][c];
      double delta_ex = delta_ex_v[0][c];
      double fixed_depth = std::pow(depr_depth,2)*(2*delta_max - 3*delta_ex)/std::pow(delta_max,2) + std::pow(depr_depth,3)*(2*delta_ex - delta_max)/std::pow(delta_max,3);
      
      if (pdd_v[0][c] <= 0.0)
        res[0][c] = 0;
      else{
        res[0][c] =  1.0 / pdd_v[0][c];
      }
    }
    
  }
  else {
    Errors::Message msg("VolumetricHeightSubgridEvaluator: Not Implemented: no derivatives implemented other than ponded depth.");
    Exceptions::amanzi_throw(msg);
  }
  
}

} //namespace
} //namespace
} //namespace
