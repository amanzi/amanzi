/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  Evaluator for determining height( rho, head )

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "overland_head_icy_water_content_evaluator2.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {


OverlandHeadIcyWaterContentEvaluator2::OverlandHeadIcyWaterContentEvaluator2(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist) {
  // my keys are for saturation and rel perm.
  my_key_ = plist_.get<string>("water content key", "surface_water_content");

  // my dependencies
  pres_key_ = plist_.get<string>("pressure key", "surface_pressure");
  dependencies_.insert(pres_key_);
  //  dependencies_.insert(std::string("surface_cell_volume"));

  M_ = plist_.get<double>("molar mass", 0.0180153);
}


OverlandHeadIcyWaterContentEvaluator2::OverlandHeadIcyWaterContentEvaluator2(const OverlandHeadIcyWaterContentEvaluator2& other) :
    SecondaryVariableFieldEvaluator(other),
    pres_key_(other.pres_key_),
    M_(other.M_) {}

Teuchos::RCP<FieldEvaluator>
OverlandHeadIcyWaterContentEvaluator2::Clone() const {
  return Teuchos::rcp(new OverlandHeadIcyWaterContentEvaluator2(*this));
}

void OverlandHeadIcyWaterContentEvaluator2::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result) {

  Epetra_MultiVector& res = *result->ViewComponent("cell",false);
  const Epetra_MultiVector& pres = *S->GetFieldData(pres_key_)
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& cv = *S->GetFieldData("surface_cell_volume")
      ->ViewComponent("cell",false);

  const double& p_atm = *S->GetScalarData("atmospheric_pressure");
  const Epetra_Vector& gravity = *S->GetConstantVectorData("gravity");
  double gz = -gravity[2];  // check this

  int ncells = res.MyLength();
  for (int c=0; c!=ncells; ++c) {
    res[0][c] = pres[0][c] < p_atm ? 0. : cv[0][c] * (pres[0][c] - p_atm) / (gz * M_);
  }
}


void OverlandHeadIcyWaterContentEvaluator2::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {

  Epetra_MultiVector& res = *result->ViewComponent("cell",false);
  const Epetra_MultiVector& pres = *S->GetFieldData(pres_key_)
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& cv = *S->GetFieldData("surface_cell_volume")
      ->ViewComponent("cell",false);

  const double& p_atm = *S->GetScalarData("atmospheric_pressure");
  const Epetra_Vector& gravity = *S->GetConstantVectorData("gravity");
  double gz = -gravity[2];  // check this

  if (wrt_key == pres_key_) {
    int ncells = res.MyLength();
    for (int c=0; c!=ncells; ++c) {
      res[0][c] = pres[0][c]  < p_atm ? 0. : cv[0][c] / (gz * M_);
    }
  } else {
    ASSERT(0);
  }
    
}

} //namespace
} //namespace
} //namespace
