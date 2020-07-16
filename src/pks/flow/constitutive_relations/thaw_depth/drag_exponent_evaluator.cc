/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  The dynamic subgrid model evaluator gets the subgrid parameters and evolve polygons.

  Authors: Ahmad Jan (jana@ornl.gov)
*/

#include "microtopographic_relief_evaluator.hh"

namespace Amanzi {
namespace Flow {
  

DepressionDepthEvaluator::DepressionDepthEvaluator(Teuchos::ParameterList& plist) :
   SecondaryVariableFieldEvaluator(plist){

  std::string domain_name=Keys::getDomain(my_key_);
  if (my_key_.empty())
    my_key_ = Keys::getKey(domain_name, "drag_exponent");
  delta_init_key_ = plist_.get<std::string>("drag exponent initial key");
  dependencies_.insert(delta_init_key_);
  delta_evolve_key_ = plist_.get<std::string>("drag exponent evolution key");
  dependencies_.insert(delta_evolve_key_);
  sg_entity_key_ = plist_.get<std::string>("subgrid entity key");
  dependencies_.insert(sg_entity_key_);
}

  
DepressionDepthEvaluator::DepressionDepthEvaluator(const DepressionDepthEvaluator& other) :
  SecondaryVariableFieldEvaluator(other),
    delta_init_key_(other.delta_init_key_),
  delta_evolve_key_(other.delta_evolve_key_),
  sg_entity_key_(other.sg_entity_key_) {}


Teuchos::RCP<FieldEvaluator>
DepressionDepthEvaluator::Clone() const {
  return Teuchos::rcp(new DepressionDepthEvaluator(*this));
}
void DepressionDepthEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result) {
 
  Epetra_MultiVector& res_c = *result->ViewComponent("cell",false);
  const Epetra_MultiVector& delta_init_c = *S->GetFieldData(delta_init_key_)->ViewComponent("cell", false);
  const Epetra_MultiVector& delta_evolve_c = *S->GetFieldData(delta_evolve_key_)->ViewComponent("cell",false);
  const Epetra_MultiVector& sg_entity_c = *S->GetFieldData(sg_entity_key_)->ViewComponent("cell",false);

  Key domain = Keys::getDomain(delta_init_key_);
  assert(!domain.empty());

  int ncells = res_c.MyLength();
  for (int c=0; c!=ncells; c++){
    if (sg_entity_c[0][c] == 1) // 1 is for HCP, 0 is for LCP
      res_c[0][c] = delta_init_c[0][c];
    else {
      res_c[0][c] = delta_init_c[0][c] + delta_evolve_c[0][c];
    }
  }

}
    
void DynamicDepressionDepthEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
               Key wrt_key, const Teuchos::Ptr<CompositeVector>& result){}

} //namespace
} //namespace



  
