/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  Evaluator for determining height( rho, head )

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "overland_head_water_content_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {


OverlandHeadWaterContentEvaluator::OverlandHeadWaterContentEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist) {
  // my keys are for saturation and rel perm.
  my_key_ = plist_.get<string>("water content key", "surface_water_content");
  setLinePrefix(my_key_+std::string(" evaluator"));

  // my dependencies
  dens_key_ = plist_.get<string>("molar density key", "surface_molar_density_liquid");
  dependencies_.insert(dens_key_);

  height_key_ = plist_.get<string>("height key", "ponded_depth");
  dependencies_.insert(height_key_);

}


OverlandHeadWaterContentEvaluator::OverlandHeadWaterContentEvaluator(const OverlandHeadWaterContentEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    dens_key_(other.dens_key_),
    height_key_(other.height_key_) {}


Teuchos::RCP<FieldEvaluator>
OverlandHeadWaterContentEvaluator::Clone() const {
  return Teuchos::rcp(new OverlandHeadWaterContentEvaluator(*this));
}


void OverlandHeadWaterContentEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result) {

  Teuchos::RCP<const CompositeVector> height = S->GetFieldData(height_key_);
  Teuchos::RCP<const CompositeVector> dens = S->GetFieldData(dens_key_);

  for (CompositeVector::name_iterator comp=result->begin();
       comp!=result->end(); ++comp) {
    Epetra_MultiVector& res_v = *result->ViewComponent(*comp,false);
    const Epetra_MultiVector& height_v = *height->ViewComponent(*comp,false);
    const Epetra_MultiVector& dens_v = *dens->ViewComponent(*comp,false);
    res_v.Multiply(1.0, height_v, dens_v, 0.);
  }
}


void OverlandHeadWaterContentEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {
  Teuchos::RCP<const CompositeVector> height = S->GetFieldData(height_key_);
  Teuchos::RCP<const CompositeVector> dens = S->GetFieldData(dens_key_);

  if (wrt_key == height_key_) {
    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      *result->ViewComponent(*comp,false) = *dens->ViewComponent(*comp,false);
    }
  } else if (wrt_key == dens_key_) {
    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      *result->ViewComponent(*comp,false) = *height->ViewComponent(*comp,false);
    }
  }
}



} //namespace
} //namespace
} //namespace
