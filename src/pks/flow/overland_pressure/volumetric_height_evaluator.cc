/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  Evaluator for determining height( rho, head )

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "volumetric_height_model.hh"
#include "volumetric_height_evaluator.hh"


namespace Amanzi {
namespace Flow {
namespace FlowRelations {


VolumetricHeightEvaluator::VolumetricHeightEvaluator(Teuchos::ParameterList& plist) :
     SecondaryVariableFieldEvaluator(plist) {
  Key domain;

  if(!my_key_.empty())
    domain = getDomain(my_key_);
  else if (my_key_.empty())
    my_key_ = plist_.get<std::string>("volumetric ponded depth key", "surface_star-volumetric_ponded_depth"); 
 
  // Key domain = getDomain(my_key_);
  // my extra dependencies
  pd_key_ = plist_.get<std::string>("height key", getKey(domain,"ponded_depth"));
  dependencies_.insert(pd_key_);

  delta_max_ = plist_.get<double>("maximum ponded depth");//,0.483);
  delta_ex_ = plist_.get<double>("excluded volume");//,0.23);
  // model
  Teuchos::ParameterList model_plist = plist_.sublist("height model parameters");
  vol_model_ = Teuchos::rcp(new VolumetricHeightModel(model_plist));

}


VolumetricHeightEvaluator::VolumetricHeightEvaluator(const VolumetricHeightEvaluator& other) :
  SecondaryVariableFieldEvaluator(other),
  pd_key_(other.pd_key_),
  delta_max_(other.delta_max_),
  delta_ex_(other.delta_ex_),
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
  bool bar_ = false;
  for (int c=0; c!=ncells; ++c){
    res_c[0][c] = std::pow(pd_c[0][c],2)*(2*delta_max_ - 3*delta_ex_)/std::pow(delta_max_,2) + std::pow(pd_c[0][c],3)*(2*delta_ex_ - delta_max_)/std::pow(delta_max_,3);
  }
  /*
  if (bar_) {
    for (int c=0; c!=ncells; ++c) {
      res_c[0][c] = vol_model_->Height(pd_c[0][c], delta_max_,delta_ex_);
    }
  } else {
    for (int c=0; c!=ncells; ++c) {
      res_c[0][c] = pd_c[0][c] <= 0.0 ? 0. :
          vol_model_->Height(pd_c[0][c], delta_max_,delta_ex_);
    }
  }
  */
}


void VolumetricHeightEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {
 
  // this is rather hacky.  surface_pressure is a mixed field vector -- it has
  // pressure on cells and ponded depth on faces.
  // -- NO FACE DERIVATIVES
  //  result->ViewComponent("face",false)->PutScalar(1.0);

  const Epetra_MultiVector& res_c = *result->ViewComponent("cell",false);
  const Epetra_MultiVector& pd_c = *S->GetFieldData(pd_key_)->ViewComponent("cell",false);

  if(wrt_key == pd_key_){
      int ncells = res_c.MyLength();
      for (int c=0; c!=ncells; ++c) {
        if (pd_c[0][c] <= delta_max_)
          res_c[0][c] = 2*pd_c[0][c]*(2*delta_max_ - 3*delta_ex_ )/ std::pow(delta_max_,2) 
            + 3*pd_c[0][c]*pd_c[0][c]*(2*delta_ex_ - delta_max_)/std::pow(delta_max_,3);
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
  /*
  // -- cells need the function eval
  const Epetra_MultiVector& res_c = *result->ViewComponent("cell",false);
  const Epetra_MultiVector& pres_c = *S->GetFieldData(pres_key_)
      ->ViewComponent("cell",false);
   const Epetra_MultiVector& rho_l = *S->GetFieldData(dens_key_)
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& rho_i = *S->GetFieldData(dens_ice_key_)
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& eta = *S->GetFieldData(unfrozen_frac_key_)
      ->ViewComponent("cell",false);

  const double& p_atm = *S->GetScalarData(patm_key_);
  const Epetra_Vector& gravity = *S->GetConstantVectorData(gravity_key_);
  double gz = -gravity[2];  // check this

  // For derivatives, the height is always assumed to be non-negative.  If it
  // is negative, the term gets zeroed later.
  if (bar_) {
    if (wrt_key == pres_key_) {
      int ncells = res_c.MyLength();
      for (int c=0; c!=ncells; ++c) {
        res_c[0][c] = icy_model_->DHeightDPressure(pres_c[0][c], eta[0][c],
                rho_l[0][c], rho_i[0][c], p_atm, gz);
      }
    } else if (wrt_key == dens_key_) {
      int ncells = res_c.MyLength();
      for (int c=0; c!=ncells; ++c) {
        res_c[0][c] = icy_model_->DHeightDRho_l(pres_c[0][c], eta[0][c],
                rho_l[0][c], rho_i[0][c], p_atm, gz);
      }
    } else if (wrt_key == dens_ice_key_) {
      int ncells = res_c.MyLength();
      for (int c=0; c!=ncells; ++c) {
        res_c[0][c] = icy_model_->DHeightDRho_i(pres_c[0][c], eta[0][c],
                rho_l[0][c], rho_i[0][c], p_atm, gz);
      }
    } else if (wrt_key == unfrozen_frac_key_) {
      int ncells = res_c.MyLength();
      for (int c=0; c!=ncells; ++c) {
        res_c[0][c] = icy_model_->DHeightDEta(pres_c[0][c], eta[0][c],
                rho_l[0][c], rho_i[0][c], p_atm, gz);
      }
    } else {
      ASSERT(0);
    }
  } else {
    if (wrt_key == pres_key_) {
      int ncells = res_c.MyLength();
      for (int c=0; c!=ncells; ++c) {
        res_c[0][c] =  pres_c[0][c] < p_atm ? 0. :
            icy_model_->DHeightDPressure(pres_c[0][c], eta[0][c],
                rho_l[0][c], rho_i[0][c], p_atm, gz);
      }
    } else if (wrt_key == dens_key_) {
      int ncells = res_c.MyLength();
      for (int c=0; c!=ncells; ++c) {
        res_c[0][c] =  pres_c[0][c] < p_atm ? 0. :
            icy_model_->DHeightDRho_l(pres_c[0][c], eta[0][c],
                rho_l[0][c], rho_i[0][c], p_atm, gz);
      }
    } else if (wrt_key == dens_ice_key_) {
      int ncells = res_c.MyLength();
      for (int c=0; c!=ncells; ++c) {
        res_c[0][c] =  pres_c[0][c] < p_atm ? 0. :
            icy_model_->DHeightDRho_i(pres_c[0][c], eta[0][c],
                rho_l[0][c], rho_i[0][c], p_atm, gz);
      }
    } else if (wrt_key == unfrozen_frac_key_) {
      int ncells = res_c.MyLength();
      for (int c=0; c!=ncells; ++c) {
        res_c[0][c] =  pres_c[0][c] < p_atm ? 0. :
            icy_model_->DHeightDEta(pres_c[0][c], eta[0][c],
                rho_l[0][c], rho_i[0][c], p_atm, gz);
      }
    } else {
      ASSERT(0);
    }
    }*/
}



} //namespace
} //namespace
} //namespace
