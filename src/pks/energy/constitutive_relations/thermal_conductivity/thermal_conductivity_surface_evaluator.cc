/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Interface for a thermal conductivity model with two phases.

  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "thermal_conductivity_surface_evaluator.hh"

namespace Amanzi {
namespace Energy {

ThermalConductivitySurfaceEvaluator::ThermalConductivitySurfaceEvaluator(
      Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist) {
  if (my_key_ == std::string("")) {
    my_key_ = plist_.get<std::string>("thermal conductivity key",
            "surface-thermal_conductivity");
  }


  Key domain = getDomain(my_key_);
  uf_key_ = plist_.get<std::string>("unfrozen fraction key", getKey(domain,"unfrozen_fraction"));
  dependencies_.insert(uf_key_);

  height_key_ = plist_.get<std::string>("height key", getKey(domain, "ponded_depth"));

  dependencies_.insert(height_key_);


  sg_model_ = plist.get<bool>("subgrid model",false);
  
  if(sg_model_){
    vpd_key_ = plist.get<std::string>("volumetric height key", getKey(domain,"volumetric_ponded_depth"));    
    dependencies_.insert(vpd_key_);
  }

  ASSERT(plist_.isSublist("thermal conductivity parameters"));
  Teuchos::ParameterList sublist = plist_.sublist("thermal conductivity parameters");
  K_liq_ = sublist.get<double>("thermal conductivity of water [W/(m-K)]", 0.58);
  K_ice_ = sublist.get<double>("thermal conductivity of ice [W/(m-K)]", 2.18);

  min_K_ = sublist.get<double>("minimum thermal conductivity", 1.e-14);
}


ThermalConductivitySurfaceEvaluator::ThermalConductivitySurfaceEvaluator(
      const ThermalConductivitySurfaceEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    uf_key_(other.uf_key_),
    height_key_(other.height_key_),
    K_liq_(other.K_liq_),
    K_ice_(other.K_ice_),
    vpd_key_(other.vpd_key_),
    min_K_(other.min_K_),
    sg_model_(other.sg_model_){}


Teuchos::RCP<FieldEvaluator>
ThermalConductivitySurfaceEvaluator::Clone() const {
  return Teuchos::rcp(new ThermalConductivitySurfaceEvaluator(*this));
}


void ThermalConductivitySurfaceEvaluator::EvaluateField_(
      const Teuchos::Ptr<State>& S,
      const Teuchos::Ptr<CompositeVector>& result) {
  // pull out the dependencies
  Teuchos::RCP<const CompositeVector> eta = S->GetFieldData(uf_key_);
  Teuchos::RCP<const CompositeVector> height = S->GetFieldData(height_key_);

  for (CompositeVector::name_iterator comp=result->begin();
       comp!=result->end(); ++comp) {
    // much more efficient to pull out vectors first
    const Epetra_MultiVector& eta_v = *eta->ViewComponent(*comp,false);
    const Epetra_MultiVector& height_v = *height->ViewComponent(*comp,false);
    Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

    int ncomp = result->size(*comp, false);
    if(!sg_model_){
      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = std::max(min_K_,
                                  height_v[0][i] * (K_liq_ * eta_v[0][i] + K_ice_ * (1. - eta_v[0][i])));
      }
    }
    else{
      const Epetra_MultiVector& vpd = *S->GetFieldData(vpd_key_)->ViewComponent("cell",false);
      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = std::max(min_K_,
                                  vpd[0][i] * (K_liq_ * eta_v[0][i] + K_ice_ * (1. - eta_v[0][i])));
      }
    }
  }

  result->Scale(1.e-6); // convert to MJ
}


void ThermalConductivitySurfaceEvaluator::EvaluateFieldPartialDerivative_(
      const Teuchos::Ptr<State>& S, Key wrt_key,
      const Teuchos::Ptr<CompositeVector>& result) {
  std::cout<<"THERMAL CONDUCITIVITY: Derivative not implemented yet!"<<wrt_key<<"\n";
  ASSERT(0); // not implemented, not yet needed
  result->Scale(1.e-6); // convert to MJ
}

} //namespace
} //namespace
