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
  M_ = plist_.get<double>("molar mass", 0.0180153);
  bar_ = plist_.get<bool>("water content bar", false);

  my_key_ = "surface_water_content";
  if (bar_) my_key_ += std::string("_bar");
  my_key_ = plist_.get<string>("water content key", my_key_);

  // my dependencies
  pres_key_ = plist_.get<string>("pressure key", "surface_pressure");
  dependencies_.insert(pres_key_);

  //  dependencies_.insert(std::string("surface_cell_volume"));
}


OverlandHeadWaterContentEvaluator::OverlandHeadWaterContentEvaluator(const OverlandHeadWaterContentEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    pres_key_(other.pres_key_),
    M_(other.M_),
    bar_(other.bar_) {}


Teuchos::RCP<FieldEvaluator>
OverlandHeadWaterContentEvaluator::Clone() const {
  return Teuchos::rcp(new OverlandHeadWaterContentEvaluator(*this));
}


void OverlandHeadWaterContentEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
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
  if (bar_) {
    for (int c=0; c!=ncells; ++c) {
      res[0][c] = cv[0][c] * (pres[0][c] - p_atm) / (gz * M_);
    }
  } else {
    for (int c=0; c!=ncells; ++c) {
      res[0][c] = pres[0][c] < p_atm ? 0. :
          cv[0][c] * (pres[0][c] - p_atm) / (gz * M_);
    }
  }
}


void OverlandHeadWaterContentEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {
  ASSERT(wrt_key == pres_key_);

  Epetra_MultiVector& res = *result->ViewComponent("cell",false);
  const Epetra_MultiVector& pres = *S->GetFieldData(pres_key_)
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& cv = *S->GetFieldData("surface_cell_volume")
      ->ViewComponent("cell",false);

  const double& p_atm = *S->GetScalarData("atmospheric_pressure");
  const Epetra_Vector& gravity = *S->GetConstantVectorData("gravity");
  double gz = -gravity[2];  // check this

  int ncells = res.MyLength();
  if (bar_) {
    for (int c=0; c!=ncells; ++c) {
      res[0][c] = cv[0][c] / (gz * M_);
    }
  } else {
    for (int c=0; c!=ncells; ++c) {
      res[0][c] = pres[0][c] < p_atm ? 0. :
          cv[0][c] / (gz * M_);
    }
  }
}


} //namespace
} //namespace
} //namespace
