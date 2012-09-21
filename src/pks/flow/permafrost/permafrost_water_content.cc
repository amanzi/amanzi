/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -----------------------------------------------------------------------------
ATS

Authors: Ethan Coon (ecoon@lanl.gov)

FieldEvaluator for water content.

Wrapping this conserved quantity as a field evaluator makes it easier to take
derivatives, keep updated, and the like.  The equation for this is simply:

WC = phi * (s_ice * n_ice + s_liquid * n_liquid + omega_gas * s_gas * n_gas)

This is simply the conserved quantity in permafrost-Richards equation.
----------------------------------------------------------------------------- */


#include "permafrost_water_content.hh"

namespace Amanzi {
namespace Flow {

PermafrostWaterContent::PermafrostWaterContent(Teuchos::ParameterList& wc_plist) {
  my_key_ = std::string("water_content");
  dependencies_.insert(std::string("porosity"));

  dependencies_.insert(std::string("saturation_liquid"));
  dependencies_.insert(std::string("molar_density_liquid"));

  dependencies_.insert(std::string("saturation_ice"));
  dependencies_.insert(std::string("molar_density_ice"));

  dependencies_.insert(std::string("saturation_gas"));
  dependencies_.insert(std::string("molar_density_gas"));
  dependencies_.insert(std::string("mol_frac_gas"));
  //  dependencies_.insert(std::string("cell_volume"));
};

PermafrostWaterContent::PermafrostWaterContent(const PermafrostWaterContent& other) :
    SecondaryVariableFieldEvaluator(other) {};

Teuchos::RCP<FieldEvaluator>
PermafrostWaterContent::Clone() const {
  return Teuchos::rcp(new PermafrostWaterContent(*this));
};


void PermafrostWaterContent::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result) {
  Teuchos::RCP<const CompositeVector> phi = S->GetFieldData("porosity");
  Teuchos::RCP<const CompositeVector> s_l = S->GetFieldData("saturation_liquid");
  Teuchos::RCP<const CompositeVector> n_l = S->GetFieldData("molar_density_liquid");
  Teuchos::RCP<const CompositeVector> s_i = S->GetFieldData("saturation_ice");
  Teuchos::RCP<const CompositeVector> n_i = S->GetFieldData("molar_density_ice");
  Teuchos::RCP<const CompositeVector> s_g = S->GetFieldData("saturation_gas");
  Teuchos::RCP<const CompositeVector> n_g = S->GetFieldData("molar_density_gas");
  Teuchos::RCP<const CompositeVector> omega_g = S->GetFieldData("mol_frac_gas");

  Teuchos::RCP<const CompositeVector> cell_volume = S->GetFieldData("cell_volume");

  for (int c=0; c!=result->size("cell"); ++c) {
    (*result)("cell",c) = (*phi)("cell",c) * ( (*s_l)("cell",c)*(*n_l)("cell",c)
            + (*s_i)("cell",c)*(*n_i)("cell",c)
            + (*s_g)("cell",c)*(*n_g)("cell",c)*(*omega_g)("cell",c) );
    (*result)("cell",c) *= (*cell_volume)("cell",c);
  }
};


void PermafrostWaterContent::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {
  Teuchos::RCP<const CompositeVector> phi = S->GetFieldData("porosity");
  Teuchos::RCP<const CompositeVector> s_l = S->GetFieldData("saturation_liquid");
  Teuchos::RCP<const CompositeVector> n_l = S->GetFieldData("molar_density_liquid");
  Teuchos::RCP<const CompositeVector> s_i = S->GetFieldData("saturation_ice");
  Teuchos::RCP<const CompositeVector> n_i = S->GetFieldData("molar_density_ice");
  Teuchos::RCP<const CompositeVector> s_g = S->GetFieldData("saturation_gas");
  Teuchos::RCP<const CompositeVector> n_g = S->GetFieldData("molar_density_gas");
  Teuchos::RCP<const CompositeVector> omega_g = S->GetFieldData("mol_frac_gas");

  Teuchos::RCP<const CompositeVector> cell_volume = S->GetFieldData("cell_volume");

  if (wrt_key == "porosity") {
    for (int c=0; c!=result->size("cell"); ++c) {
      (*result)("cell",c) = (*s_l)("cell",c)*(*n_l)("cell",c)
          + (*s_g)("cell",c)*(*n_g)("cell",c)*(*omega_g)("cell",c);
    }
  } else if (wrt_key == "saturation_liquid") {
    for (int c=0; c!=result->size("cell"); ++c) {
      (*result)("cell",c) = (*phi)("cell",c) * (*n_l)("cell",c);
    }
  } else if (wrt_key == "molar_density_liquid") {
    for (int c=0; c!=result->size("cell"); ++c) {
      (*result)("cell",c) = (*phi)("cell",c) * (*s_l)("cell",c);
    }
  } else if (wrt_key == "saturation_ice") {
    for (int c=0; c!=result->size("cell"); ++c) {
      (*result)("cell",c) = (*phi)("cell",c) * (*n_i)("cell",c);
    }
  } else if (wrt_key == "molar_density_ice") {
    for (int c=0; c!=result->size("cell"); ++c) {
      (*result)("cell",c) = (*phi)("cell",c) * (*s_i)("cell",c);
    }
  } else if (wrt_key == "saturation_gas") {
    for (int c=0; c!=result->size("cell"); ++c) {
      (*result)("cell",c) = (*phi)("cell",c) * (*n_g)("cell",c)*(*omega_g)("cell",c);
    }
  } else if (wrt_key == "molar_density_gas") {
    for (int c=0; c!=result->size("cell"); ++c) {
      (*result)("cell",c) = (*phi)("cell",c) * (*s_g)("cell",c)*(*omega_g)("cell",c);
    }
  } else if (wrt_key == "mol_frac_gas") {
    for (int c=0; c!=result->size("cell"); ++c) {
      (*result)("cell",c) = (*phi)("cell",c) * (*s_g)("cell",c)*(*n_g)("cell",c);
    }
  }

  for (int c=0; c!=result->size("cell"); ++c) {
    (*result)("cell",c) *= (*cell_volume)("cell",c);
  }
};


} //namespace
} //namespace
