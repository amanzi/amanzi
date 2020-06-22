/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Evaluates the conductivity of surface flow.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "independent_variable_field_evaluator.hh"
#include "overland_conductivity_evaluator.hh"
#include "manning_conductivity_model.hh"
#include "split_denominator_conductivity_model.hh"
#include "ponded_depth_passthrough_conductivity_model.hh"

namespace Amanzi {
namespace Flow {

OverlandConductivityEvaluator::OverlandConductivityEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist) {
  Key domain = Keys::getDomain(my_key_);

  if (plist_.isParameter("height key") || plist_.isParameter("ponded depth key")
      || plist_.isParameter("height key suffix") || plist_.isParameter("ponded depth key suffix")) {
    Errors::Message message("OverlandConductivity: only use \"depth key\" or \"depth key suffix\", not \"height key\" or \"ponded depth key\".");
    Exceptions::amanzi_throw(message);
  }
  depth_key_ = Keys::readKey(plist_, domain, "depth", "ponded_depth");
  dependencies_.insert(depth_key_);

  slope_key_ = Keys::readKey(plist_, domain, "slope", "slope_magnitude");
  dependencies_.insert(slope_key_);

  coef_key_ = Keys::readKey(plist_, domain, "coefficient", "manning_coefficient");
  dependencies_.insert(coef_key_);

  dt_ = plist_.get<bool>("include dt factor", false);
  if (dt_) {
    factor_ = plist_.get<double>("dt factor");
  }
  dens_ = plist_.get<bool>("include density factor", true);

  if (dens_) {
    dens_key_ = Keys::readKey(plist_, domain, "molar density liquid", "molar_density_liquid");
    dependencies_.insert(dens_key_);
  }

  sg_model_ =  plist_.get<bool>("subgrid model", false);
  AMANZI_ASSERT(!sg_model_);
  if(sg_model_){
    pdd_key_ = Keys::readKey(plist_, domain, "ponded depth minus depression depth", "ponded_depth_minus_depression_depth");
    dependencies_.insert(pdd_key_);
    
    frac_cond_key_ = Keys::readKey(plist_, domain, "fractional conductance", "fractional_conductance");
    dependencies_.insert(frac_cond_key_); 

    vpd_key_ = Keys::readKey(plist_, domain, "volumetric ponded depth", "volumetric_ponded_depth");
    dependencies_.insert(vpd_key_); 

    drag_exp_key_ = Keys::readKey(plist_, domain, "drag exponent", "drag_exponent");
    dependencies_.insert(drag_exp_key_); 
  }

  // create the model
  AMANZI_ASSERT(plist_.isSublist("overland conductivity model"));
  Teuchos::ParameterList sublist = plist_.sublist("overland conductivity model");
  std::string model_type = sublist.get<std::string>("overland conductivity type", "manning");
  if ((model_type == "manning") || (model_type == "manning harmonic mean") ||
      (model_type == "manning cell centered")){
    model_ = Teuchos::rcp(new ManningConductivityModel(sublist));
  } else if (model_type == "manning split denominator") {
    model_ = Teuchos::rcp(new SplitDenominatorConductivityModel(sublist));
  } else if (model_type == "manning ponded depth passthrough") {
    model_ = Teuchos::rcp(new PondedDepthPassthroughConductivityModel(sublist));
  } else {
    AMANZI_ASSERT(0);
  }

  if (sg_model_ && model_type != "manning") {
    Errors::Message message("Subgrid Overland Conductivity: only 'manning' model type is implemented with subgrid models!");
    Exceptions::amanzi_throw(message);
  }
  
}


OverlandConductivityEvaluator::OverlandConductivityEvaluator(const OverlandConductivityEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    depth_key_(other.depth_key_),
    slope_key_(other.slope_key_),
    coef_key_(other.coef_key_),
    dens_key_(other.dens_key_),
    dt_(other.dt_),
    factor_(other.factor_),
    dens_(other.dens_),
    model_(other.model_),
    pdd_key_(other.pdd_key_),
    vpd_key_(other.vpd_key_),
    frac_cond_key_(other.frac_cond_key_),
    sg_model_(other.sg_model_),
    drag_exp_key_(other.drag_exp_key_) {}


Teuchos::RCP<FieldEvaluator>
OverlandConductivityEvaluator::Clone() const {
  return Teuchos::rcp(new OverlandConductivityEvaluator(*this));
}


// Required methods from SecondaryVariableFieldEvaluator
void OverlandConductivityEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result) {

  Teuchos::RCP<const CompositeVector> depth = S->GetFieldData(depth_key_);
  Teuchos::RCP<const CompositeVector> slope = S->GetFieldData(slope_key_);
  Teuchos::RCP<const CompositeVector> coef = S->GetFieldData(coef_key_);

#ifdef ENABLE_DBC  
  double min_coef = 1.;
  coef->MinValue(&min_coef);
  if (min_coef <= 0.) {
    Errors::Message message("Overland Conductivity Evaluator: Manning coeficient has at least one value that is non-positive.  Perhaps you forgot to set the \"boundary_face\" component?"); 
    Exceptions::amanzi_throw(message);
  }
#endif

  if (!sg_model_){
    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      const Epetra_MultiVector& depth_v = *depth->ViewComponent(*comp,false);
      const Epetra_MultiVector& slope_v = *slope->ViewComponent(*comp,false);
      const Epetra_MultiVector& coef_v = *coef->ViewComponent(*comp,false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

      int ncomp = result->size(*comp, false);
      if (dt_ && factor_ < 0.) {
        double dt = S->time() - S->last_time();
        for (int i=0; i!=ncomp; ++i) {
          result_v[0][i] = model_->Conductivity(dt*depth_v[0][i], slope_v[0][i], coef_v[0][i]);
        }
      } else if (dt_) {
        double dt = S->time() - S->last_time();
        for (int i=0; i!=ncomp; ++i) {
          result_v[0][i] = model_->Conductivity(10*factor_*depth_v[0][i], slope_v[0][i], coef_v[0][i]);
        }
      } else {
        for (int i=0; i!=ncomp; ++i) {
          result_v[0][i] = model_->Conductivity(depth_v[0][i], slope_v[0][i], coef_v[0][i]);
        }
      }
      if (dens_) {
        const Epetra_MultiVector& dens_v = *S->GetFieldData(dens_key_)->ViewComponent(*comp,false);
        for (int i=0; i!=ncomp; ++i)
          result_v[0][i] *= dens_v[0][i];
      }
    }

    
  }
  else{
    Teuchos::RCP<const CompositeVector> pd_depth = S->GetFieldData(pdd_key_);
    Teuchos::RCP<const CompositeVector> frac_cond = S->GetFieldData(frac_cond_key_);
    Teuchos::RCP<const CompositeVector> vpd = S->GetFieldData(vpd_key_);

    Teuchos::RCP<const CompositeVector> drag = S->GetFieldData(drag_exp_key_);

    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      const Epetra_MultiVector& depth_v = *depth->ViewComponent(*comp,false);
      
      const Epetra_MultiVector& pd_depth_v = *pd_depth->ViewComponent(*comp,false);
      const Epetra_MultiVector& frac_cond_v = *frac_cond->ViewComponent(*comp,false);
      const Epetra_MultiVector& vpd_v = *vpd->ViewComponent(*comp,false);
      const Epetra_MultiVector& drag_v = *drag->ViewComponent(*comp,false);
      
      const Epetra_MultiVector& slope_v = *slope->ViewComponent(*comp,false);
      const Epetra_MultiVector& coef_v = *coef->ViewComponent(*comp,false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);
      
      int ncomp = result->size(*comp, false);
      if (dt_ && factor_ < 0.) {
        double dt = S->time() - S->last_time();
        for (int i=0; i!=ncomp; ++i) {
          result_v[0][i] = model_->Conductivity(dt*depth_v[0][i], slope_v[0][i],coef_v[0][i], pd_depth_v[0][i],frac_cond_v[0][i],drag_v[0][i]);
        }
      } else if (dt_) {
        double dt = S->time() - S->last_time();
        for (int i=0; i!=ncomp; ++i) {
          result_v[0][i] = model_->Conductivity(10*factor_*depth_v[0][i], slope_v[0][i],coef_v[0][i], pd_depth_v[0][i],frac_cond_v[0][i],drag_v[0][i]);
        }
      } else {
        for (int i=0; i!=ncomp; ++i) {
          result_v[0][i] = model_->Conductivity(depth_v[0][i], slope_v[0][i], coef_v[0][i], pd_depth_v[0][i], frac_cond_v[0][i],drag_v[0][i]);
        }
      }
      if (dens_) {
        const Epetra_MultiVector& dens_v = *S->GetFieldData(dens_key_)->ViewComponent(*comp,false);
        for (int i=0; i!=ncomp; ++i)
          result_v[0][i] *= dens_v[0][i];
      }
    }
  }
}


void OverlandConductivityEvaluator::EvaluateFieldPartialDerivative_(
    const Teuchos::Ptr<State>& S,
    Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {
   //never called otherwsise it needs to changed for subgrid model
  if (sg_model_){
    Errors::Message message("Overland Conductivity Evaluator: Evaluate partial derivaritve not implemented for the Subgrid Model."); 
    Exceptions::amanzi_throw(message);
  }
  Teuchos::RCP<const CompositeVector> depth = S->GetFieldData(depth_key_);
  Teuchos::RCP<const CompositeVector> slope = S->GetFieldData(slope_key_);
  Teuchos::RCP<const CompositeVector> coef = S->GetFieldData(coef_key_);

  if (wrt_key == depth_key_) {
    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      const Epetra_MultiVector& depth_v = *depth->ViewComponent(*comp,false);
      const Epetra_MultiVector& slope_v = *slope->ViewComponent(*comp,false);
      const Epetra_MultiVector& coef_v = *coef->ViewComponent(*comp,false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

      int ncomp = result->size(*comp, false);
      if (dt_ && factor_ < 0.) {
        double dt = S->time() - S->last_time();
        for (int i=0; i!=ncomp; ++i) {
          result_v[0][i] = model_->DConductivityDDepth(dt*depth_v[0][i], slope_v[0][i], coef_v[0][i]) * dt;
        }
      } else if (dt_) {
        double dt = S->time() - S->last_time();
        for (int i=0; i!=ncomp; ++i) {
          result_v[0][i] = model_->DConductivityDDepth(10*factor_*depth_v[0][i], slope_v[0][i], coef_v[0][i]) * 10*factor_;
        }
      } else {
        for (int i=0; i!=ncomp; ++i) {
          result_v[0][i] = model_->DConductivityDDepth(depth_v[0][i], slope_v[0][i], coef_v[0][i]);
        }
      }

      if (dens_) {
        const Epetra_MultiVector& dens_v = *S->GetFieldData(dens_key_)->ViewComponent(*comp,false);
        for (int i=0; i!=ncomp; ++i) {
          result_v[0][i] *= dens_v[0][i];
        }
      }
    }

  } else if (wrt_key == dens_key_) {
    AMANZI_ASSERT(dens_);
    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      const Epetra_MultiVector& depth_v = *depth->ViewComponent(*comp,false);
      const Epetra_MultiVector& slope_v = *slope->ViewComponent(*comp,false);
      const Epetra_MultiVector& coef_v = *coef->ViewComponent(*comp,false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

      int ncomp = result->size(*comp, false);
      if (dt_ && factor_ < 0.) {
        double dt = S->time() - S->last_time();
        for (int i=0; i!=ncomp; ++i) {
          result_v[0][i] = model_->Conductivity(dt*depth_v[0][i], slope_v[0][i], coef_v[0][i]) * dt;
        }
      } else if (dt_) {
        double dt = S->time() - S->last_time();
        for (int i=0; i!=ncomp; ++i) {
          result_v[0][i] = model_->Conductivity(10*factor_*depth_v[0][i], slope_v[0][i], coef_v[0][i]) * 10*factor_;
        }
      } else {
        for (int i=0; i!=ncomp; ++i) {
          result_v[0][i] = model_->Conductivity(depth_v[0][i], slope_v[0][i], coef_v[0][i]);
        }
      }
    }

  } else {
    // FIX ME -- need to add derivatives of conductivity model wrt slope, coef --etc
    result->PutScalar(0.);
  }
}


} //namespace
} //namespace

