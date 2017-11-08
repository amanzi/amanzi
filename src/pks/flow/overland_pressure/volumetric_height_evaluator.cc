/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  Evaluator for determining height( rho, head )

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "volumetric_height_model.hh"
#include "volumetric_height_evaluator.hh"


namespace Amanzi {
namespace Flow {


VolumetricHeightEvaluator::VolumetricHeightEvaluator(Teuchos::ParameterList& plist) :
     SecondaryVariableFieldEvaluator(plist) {
  Key domain;

  if(!my_key_.empty())
    domain = Keys::getDomain(my_key_);
  else if (my_key_.empty())
    my_key_ = plist_.get<std::string>("volumetric ponded depth key", "surface_star-volumetric_ponded_depth"); 
 
  // Key domain = Keys::getDomain(my_key_);
  // my extra dependencies
  pd_key_ = plist_.get<std::string>("height key", Keys::getKey(domain,"ponded_depth"));
  dependencies_.insert(pd_key_);

  delta_max_key_ = plist_.get<std::string>("maximum ponded depth key", Keys::getKey(domain,"maximum_ponded_depth"));
  dependencies_.insert(delta_max_key_);
  delta_ex_key_ = plist_.get<std::string>("excluded volume key", Keys::getKey(domain,"excluded_volume"));
  dependencies_.insert(delta_ex_key_);
  // model
  Teuchos::ParameterList model_plist = plist_.sublist("height model parameters");
  vol_model_ = Teuchos::rcp(new VolumetricHeightModel(model_plist));

}


VolumetricHeightEvaluator::VolumetricHeightEvaluator(const VolumetricHeightEvaluator& other) :
  SecondaryVariableFieldEvaluator(other),
  pd_key_(other.pd_key_),
  delta_max_key_(other.delta_max_key_),
  delta_ex_key_(other.delta_ex_key_),
  vol_model_(other.vol_model_) {}


Teuchos::RCP<FieldEvaluator>
VolumetricHeightEvaluator::Clone() const {
  return Teuchos::rcp(new VolumetricHeightEvaluator(*this));
}


void VolumetricHeightEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result) {

  Epetra_MultiVector& res_c = *result->ViewComponent("cell",false);
  const Epetra_MultiVector& pd_c = *S->GetFieldData(pd_key_)->ViewComponent("cell",false);

  int ncells = res_c.MyLength();

  Teuchos::RCP<const CompositeVector> max_pd = S->GetFieldData(delta_max_key_);
  Teuchos::RCP<const CompositeVector> ex_vol = S->GetFieldData(delta_ex_key_);
  // cell values
  const Epetra_MultiVector& max_pd_v = *max_pd->ViewComponent("cell", false);
  const Epetra_MultiVector& ex_vol_v = *ex_vol->ViewComponent("cell", false);

  for (int c=0; c!=ncells; ++c){
    double delta_max = max_pd_v[0][c];
    double delta_ex = ex_vol_v[0][c];
    res_c[0][c] = std::pow(pd_c[0][c],2)*(2*delta_max - 3*delta_ex)/std::pow(delta_max,2) + std::pow(pd_c[0][c],3)*(2*delta_ex - delta_max)/std::pow(delta_max,3);
  }

}


void VolumetricHeightEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {
 
  // this is rather hacky.  surface_pressure is a mixed field vector -- it has
  // pressure on cells and ponded depth on faces.
  // -- NO FACE DERIVATIVES
  //  result->ViewComponent("face",false)->PutScalar(1.0);

  const Epetra_MultiVector& res_c = *result->ViewComponent("cell",false);
  const Epetra_MultiVector& pd_c = *S->GetFieldData(pd_key_)->ViewComponent("cell",false);

  Teuchos::RCP<const CompositeVector> max_pd = S->GetFieldData(delta_max_key_);
  Teuchos::RCP<const CompositeVector> ex_vol = S->GetFieldData(delta_ex_key_);
  // cell values
  const Epetra_MultiVector& max_pd_v = *max_pd->ViewComponent("cell", false);
  const Epetra_MultiVector& ex_vol_v = *ex_vol->ViewComponent("cell", false);

  if(wrt_key == pd_key_){
      int ncells = res_c.MyLength();
      for (int c=0; c!=ncells; ++c) {
        double delta_max = max_pd_v[0][c];
        double delta_ex = ex_vol_v[0][c];
        if (pd_c[0][c] <= delta_max)
          res_c[0][c] = 2*pd_c[0][c]*(2*delta_max - 3*delta_ex )/ std::pow(delta_max,2) 
            + 3*pd_c[0][c]*pd_c[0][c]*(2*delta_ex - delta_max)/std::pow(delta_max,3);
        else
          res_c[0][c] = 1.0;
        //        res_c[0][c] = icy_model_->DHeightDPressure(pres_c[0][c], eta[0][c],
        //        rho_l[0][c], rho_i[0][c], p_atm, gz);
      }
    }
  else{
    std::cout<<"VOLUMETRIC HEIGHT EVALUATOR: NO DERIVATIVE EXISTS: "<<wrt_key<<"\n";
    ASSERT(0);
  }

}


} //namespace
} //namespace
