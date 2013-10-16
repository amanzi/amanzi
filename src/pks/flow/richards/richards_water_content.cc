/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -----------------------------------------------------------------------------
ATS

Authors: Ethan Coon (ecoon@lanl.gov)

FieldEvaluator for water content.

Wrapping this conserved quantity as a field evaluator makes it easier to take
derivatives, keep updated, and the like.  The equation for this is simply:

WC = phi * (s_liquid * n_liquid + omega_gas * s_gas * n_gas)

This is simply the conserved quantity in Richards equation.
----------------------------------------------------------------------------- */


#include "richards_water_content.hh"

namespace Amanzi {
namespace Flow {

RichardsWaterContent::RichardsWaterContent(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist) {
  my_key_ = std::string("water_content");
  setLinePrefix(my_key_+std::string(" evaluator"));

  dependencies_.insert(std::string("porosity"));

  dependencies_.insert(std::string("saturation_liquid"));
  dependencies_.insert(std::string("molar_density_liquid"));

  dependencies_.insert(std::string("saturation_gas"));
  dependencies_.insert(std::string("molar_density_gas"));
  dependencies_.insert(std::string("mol_frac_gas"));
  dependencies_.insert(std::string("cell_volume"));
};

RichardsWaterContent::RichardsWaterContent(const RichardsWaterContent& other) :
    SecondaryVariableFieldEvaluator(other) {};

Teuchos::RCP<FieldEvaluator>
RichardsWaterContent::Clone() const {
  return Teuchos::rcp(new RichardsWaterContent(*this));
};


void RichardsWaterContent::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result) {
  const Epetra_MultiVector& s_l = *S->GetFieldData("saturation_liquid")->ViewComponent("cell",false);
  const Epetra_MultiVector& n_l = *S->GetFieldData("molar_density_liquid")->ViewComponent("cell",false);

  const Epetra_MultiVector& s_g = *S->GetFieldData("saturation_gas")->ViewComponent("cell",false);
  const Epetra_MultiVector& n_g = *S->GetFieldData("molar_density_gas")->ViewComponent("cell",false);
  const Epetra_MultiVector& omega_g = *S->GetFieldData("mol_frac_gas")->ViewComponent("cell",false);

  const Epetra_MultiVector& phi = *S->GetFieldData("porosity")->ViewComponent("cell",false);
  const Epetra_MultiVector& cell_volume = *S->GetFieldData("cell_volume")->ViewComponent("cell",false);
  Epetra_MultiVector& result_v = *result->ViewComponent("cell",false);

  int ncells = result->size("cell",false);
  for (int c=0; c!=ncells; ++c) {
    result_v[0][c] = phi[0][c] * ( s_l[0][c]*n_l[0][c]
            + s_g[0][c]*n_g[0][c]*omega_g[0][c] );
    result_v[0][c] *= cell_volume[0][c];
  }
};


void RichardsWaterContent::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {
  const Epetra_MultiVector& s_l = *S->GetFieldData("saturation_liquid")->ViewComponent("cell",false);
  const Epetra_MultiVector& n_l = *S->GetFieldData("molar_density_liquid")->ViewComponent("cell",false);

  const Epetra_MultiVector& s_g = *S->GetFieldData("saturation_gas")->ViewComponent("cell",false);
  const Epetra_MultiVector& n_g = *S->GetFieldData("molar_density_gas")->ViewComponent("cell",false);
  const Epetra_MultiVector& omega_g = *S->GetFieldData("mol_frac_gas")->ViewComponent("cell",false);

  const Epetra_MultiVector& phi = *S->GetFieldData("porosity")->ViewComponent("cell",false);
  const Epetra_MultiVector& cell_volume = *S->GetFieldData("cell_volume")->ViewComponent("cell",false);
  Epetra_MultiVector& result_v = *result->ViewComponent("cell",false);

  int ncells = result->size("cell",false);
  if (wrt_key == "porosity") {
    for (int c=0; c!=ncells; ++c) {
      result_v[0][c] = s_l[0][c]*n_l[0][c]
          + s_g[0][c]*n_g[0][c]*omega_g[0][c];
    }
  } else if (wrt_key == "saturation_liquid") {
    for (int c=0; c!=ncells; ++c) {
      result_v[0][c] = phi[0][c] * n_l[0][c];
    }
  } else if (wrt_key == "molar_density_liquid") {
    for (int c=0; c!=ncells; ++c) {
      result_v[0][c] = phi[0][c] * s_l[0][c];
    }
  } else if (wrt_key == "saturation_gas") {
    for (int c=0; c!=ncells; ++c) {
      result_v[0][c] = phi[0][c] * n_g[0][c]*omega_g[0][c];
    }
  } else if (wrt_key == "molar_density_gas") {
    for (int c=0; c!=ncells; ++c) {
      result_v[0][c] = phi[0][c] * s_g[0][c]*omega_g[0][c];
    }
  } else if (wrt_key == "mol_frac_gas") {
    for (int c=0; c!=ncells; ++c) {
      result_v[0][c] = phi[0][c] * s_g[0][c]*n_g[0][c];
    }
  }

  for (int c=0; c!=ncells; ++c) {
    result_v[0][c] *= cell_volume[0][c];
  }
};


} //namespace
} //namespace
