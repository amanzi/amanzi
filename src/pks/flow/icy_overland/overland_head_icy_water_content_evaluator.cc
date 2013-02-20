/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  Evaluator for determining height( rho, head )

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "overland_head_icy_water_content_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {


OverlandHeadIcyWaterContentEvaluator::OverlandHeadIcyWaterContentEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist) {
  // my keys are for saturation and rel perm.
  my_key_ = plist_.get<string>("water content key", "surface_water_content");
  setLinePrefix(my_key_+std::string(" evaluator"));

  // my dependencies
  dens_key_ = plist_.get<string>("liquid molar density key", "surface_molar_density_liquid");
  dependencies_.insert(dens_key_);

  dens_ice_key_ = plist_.get<string>("ice molar density key", "surface_molar_density_ice");
  dependencies_.insert(dens_ice_key_);

  unfrozen_frac_key_ = plist_.get<string>("unfrozen fraction key", "unfrozen_fraction");
  dependencies_.insert(unfrozen_frac_key_);

  height_key_ = plist_.get<string>("height key", "ponded_depth");
  dependencies_.insert(height_key_);

  //  dependencies_.insert(std::string("surface_cell_volume"));
}


OverlandHeadIcyWaterContentEvaluator::OverlandHeadIcyWaterContentEvaluator(const OverlandHeadIcyWaterContentEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    dens_key_(other.dens_key_),
    dens_ice_key_(other.dens_ice_key_),
    height_key_(other.height_key_),
    unfrozen_frac_key_(other.unfrozen_frac_key_) {}


Teuchos::RCP<FieldEvaluator>
OverlandHeadIcyWaterContentEvaluator::Clone() const {
  return Teuchos::rcp(new OverlandHeadIcyWaterContentEvaluator(*this));
}


void OverlandHeadIcyWaterContentEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result) {

  Epetra_MultiVector& res = *result->ViewComponent("cell",false);
  const Epetra_MultiVector& height = *S->GetFieldData(height_key_)
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& dens_l = *S->GetFieldData(dens_key_)
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& dens_i = *S->GetFieldData(dens_ice_key_)
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& eta = *S->GetFieldData(unfrozen_frac_key_)
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& cv = *S->GetFieldData("surface_cell_volume")
      ->ViewComponent("cell",false);

  int ncells = res.MyLength();
  for (int c=0; c!=ncells; ++c) {
    res[0][c] = cv[0][c] * height[0][c]
        * (eta[0][c] * dens_l[0][c] + (1.-eta[0][c]) * dens_i[0][c]);
  }
}


void OverlandHeadIcyWaterContentEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {

  Epetra_MultiVector& res = *result->ViewComponent("cell",false);
  const Epetra_MultiVector& height = *S->GetFieldData(height_key_)
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& dens_l = *S->GetFieldData(dens_key_)
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& dens_i = *S->GetFieldData(dens_ice_key_)
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& eta = *S->GetFieldData(unfrozen_frac_key_)
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& cv = *S->GetFieldData("surface_cell_volume")
      ->ViewComponent("cell",false);

  if (wrt_key == height_key_) {
    int ncells = res.MyLength();
    for (int c=0; c!=ncells; ++c) {
      res[0][c] = cv[0][c] * (eta[0][c] * dens_l[0][c] + (1.-eta[0][c]) * dens_i[0][c]);
    }
  } else if (wrt_key == dens_key_) {
    int ncells = res.MyLength();
    for (int c=0; c!=ncells; ++c) {
      res[0][c] = cv[0][c] * height[0][c] * eta[0][c];
    }
  } else if (wrt_key == dens_ice_key_) {
    int ncells = res.MyLength();
    for (int c=0; c!=ncells; ++c) {
      res[0][c] = cv[0][c] * height[0][c] * (1. - eta[0][c]);
    }
  } else if (wrt_key == unfrozen_frac_key_) {
    int ncells = res.MyLength();
    for (int c=0; c!=ncells; ++c) {
      res[0][c] = cv[0][c] * height[0][c] * (dens_l[0][c] - dens_i[0][c]);
    }
  }
}

} //namespace
} //namespace
} //namespace
