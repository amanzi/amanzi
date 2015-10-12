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
namespace Relations {

PermafrostWaterContent::PermafrostWaterContent(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist) {
  //my_key_ = std::string("water_content");
  
  //Key domain_name = getDomainPrefix(my_key_);
  Key domain_name = getDomain(my_key_);
  phi_key = std::string(getKey(domain_name, "porosity"));
  sl_key = std::string(getKey(domain_name,"saturation_liquid"));
  mdl_key = std::string(getKey(domain_name,"molar_density_liquid"));
  si_key = std::string(getKey(domain_name,"saturation_ice"));
  mdi_key = std::string(getKey(domain_name,"molar_density_ice"));
  sg_key = std::string(getKey(domain_name,"saturation_gas"));
  mdg_key = std::string(getKey(domain_name,"molar_density_gas"));
  mfg_key = std::string(getKey(domain_name,"mol_frac_gas"));
  cv_key = std::string(getKey(domain_name,"cell_volume"));

  dependencies_.insert(std::string(phi_key));

  dependencies_.insert(std::string(sl_key));
  dependencies_.insert(std::string(mdl_key));

  dependencies_.insert(std::string(si_key));
  dependencies_.insert(std::string(mdi_key));

  dependencies_.insert(std::string(sg_key));
  dependencies_.insert(std::string(mdg_key));
  dependencies_.insert(std::string(mfg_key));
  dependencies_.insert(std::string(cv_key));

  //  check_derivative_ = true;
  
};

PermafrostWaterContent::PermafrostWaterContent(const PermafrostWaterContent& other) :
    SecondaryVariableFieldEvaluator(other) {};

Teuchos::RCP<FieldEvaluator>
PermafrostWaterContent::Clone() const {
  return Teuchos::rcp(new PermafrostWaterContent(*this));
};


void PermafrostWaterContent::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result) {
  Key domain_name = getDomain(my_key_);
  phi_key = std::string(getKey(domain_name, "porosity"));
  sl_key = std::string(getKey(domain_name,"saturation_liquid"));
  mdl_key = std::string(getKey(domain_name,"molar_density_liquid"));
  si_key = std::string(getKey(domain_name,"saturation_ice"));
  mdi_key = std::string(getKey(domain_name,"molar_density_ice"));
  sg_key = std::string(getKey(domain_name,"saturation_gas"));
  mdg_key = std::string(getKey(domain_name,"molar_density_gas"));
  mfg_key = std::string(getKey(domain_name,"mol_frac_gas"));
  cv_key = std::string(getKey(domain_name,"cell_volume"));

  const Epetra_MultiVector& s_l = *S->GetFieldData(sl_key)->ViewComponent("cell",false);
  const Epetra_MultiVector& n_l = *S->GetFieldData(mdl_key)->ViewComponent("cell",false);

  const Epetra_MultiVector& s_i = *S->GetFieldData(si_key)->ViewComponent("cell",false);
  const Epetra_MultiVector& n_i = *S->GetFieldData(mdi_key)->ViewComponent("cell",false);

  const Epetra_MultiVector& s_g = *S->GetFieldData(sg_key)->ViewComponent("cell",false);
  const Epetra_MultiVector& n_g = *S->GetFieldData(mdg_key)->ViewComponent("cell",false);
  const Epetra_MultiVector& omega_g = *S->GetFieldData(mfg_key)->ViewComponent("cell",false);

  const Epetra_MultiVector& phi = *S->GetFieldData(phi_key)->ViewComponent("cell",false);
  const Epetra_MultiVector& cell_volume = *S->GetFieldData(cv_key)->ViewComponent("cell",false);
  Epetra_MultiVector& result_v = *result->ViewComponent("cell",false);

  /* const Epetra_MultiVector& s_l = *S->GetFieldData(sl_key)->ViewComponent("cell",false);
  const Epetra_MultiVector& n_l = *S->GetFieldData(mdl_key)->ViewComponent("cell",false);

  const Epetra_MultiVector& s_i = *S->GetFieldData(si_key)->ViewComponent("cell",false);
  const Epetra_MultiVector& n_i = *S->GetFieldData(mdi_key)->ViewComponent("cell",false);

  const Epetra_MultiVector& s_g = *S->GetFieldData(sg_key)->ViewComponent("cell",false);
  const Epetra_MultiVector& n_g = *S->GetFieldData(mdg_key)->ViewComponent("cell",false);
  const Epetra_MultiVector& omega_g = *S->GetFieldData(mfg_key)->ViewComponent("cell",false);

  const Epetra_MultiVector& phi = *S->GetFieldData(phi_key)->ViewComponent("cell",false);
  const Epetra_MultiVector& cell_volume = *S->GetFieldData("cell_volume")->ViewComponent("cell",false);
  Epetra_MultiVector& result_v = *result->ViewComponent("cell",false);
  */
  int ncells = result->size("cell",false);
  for (int c=0; c!=ncells; ++c) {
    result_v[0][c] = phi[0][c] * ( s_l[0][c]*n_l[0][c]
            + s_i[0][c]*n_i[0][c]
            + s_g[0][c]*n_g[0][c]*omega_g[0][c] );
    result_v[0][c] *= cell_volume[0][c];
  }
};


void PermafrostWaterContent::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {
  /*
const Epetra_MultiVector& s_l = *S->GetFieldData("saturation_liquid")->ViewComponent("cell",false);
  const Epetra_MultiVector& n_l = *S->GetFieldData("molar_density_liquid")->ViewComponent("cell",false);

  const Epetra_MultiVector& s_i = *S->GetFieldData("saturation_ice")->ViewComponent("cell",false);
  const Epetra_MultiVector& n_i = *S->GetFieldData("molar_density_ice")->ViewComponent("cell",false);

  const Epetra_MultiVector& s_g = *S->GetFieldData("saturation_gas")->ViewComponent("cell",false);
  const Epetra_MultiVector& n_g = *S->GetFieldData("molar_density_gas")->ViewComponent("cell",false);
  const Epetra_MultiVector& omega_g = *S->GetFieldData("mol_frac_gas")->ViewComponent("cell",false);

  const Epetra_MultiVector& phi = *S->GetFieldData("porosity")->ViewComponent("cell",false);
  const Epetra_MultiVector& cell_volume = *S->GetFieldData("cell_volume")->ViewComponent("cell",false);
  Epetra_MultiVector& result_v = *result->ViewComponent("cell",false);
  */
  Key domain_name = getDomain(my_key_);
  phi_key = std::string(getKey(domain_name, "porosity"));
  sl_key = std::string(getKey(domain_name,"saturation_liquid"));
  mdl_key = std::string(getKey(domain_name,"molar_density_liquid"));
  si_key = std::string(getKey(domain_name,"saturation_ice"));
  mdi_key = std::string(getKey(domain_name,"molar_density_ice"));
  sg_key = std::string(getKey(domain_name,"saturation_gas"));
  mdg_key = std::string(getKey(domain_name,"molar_density_gas"));
  mfg_key = std::string(getKey(domain_name,"mol_frac_gas"));
  cv_key = std::string(getKey(domain_name,"cell_volume"));
  const Epetra_MultiVector& s_l = *S->GetFieldData(sl_key)->ViewComponent("cell",false);
  const Epetra_MultiVector& n_l = *S->GetFieldData(mdl_key)->ViewComponent("cell",false);
  
  const Epetra_MultiVector& s_i = *S->GetFieldData(si_key)->ViewComponent("cell",false);
  const Epetra_MultiVector& n_i = *S->GetFieldData(mdi_key)->ViewComponent("cell",false);
  
  const Epetra_MultiVector& s_g = *S->GetFieldData(sg_key)->ViewComponent("cell",false);
  const Epetra_MultiVector& n_g = *S->GetFieldData(mdg_key)->ViewComponent("cell",false);
  const Epetra_MultiVector& omega_g = *S->GetFieldData(mfg_key)->ViewComponent("cell",false);

  const Epetra_MultiVector& phi = *S->GetFieldData(phi_key)->ViewComponent("cell",false);
  const Epetra_MultiVector& cell_volume = *S->GetFieldData(cv_key)->ViewComponent("cell",false);
  Epetra_MultiVector& result_v = *result->ViewComponent("cell",false);

  int ncells = result->size("cell",false);
  if (wrt_key == phi_key) {
    for (int c=0; c!=ncells; ++c) {
      result_v[0][c] = s_l[0][c]*n_l[0][c]
          + s_i[0][c]*n_i[0][c]
          + s_g[0][c]*n_g[0][c]*omega_g[0][c];
    }
  } else if (wrt_key == sl_key) {
    for (int c=0; c!=ncells; ++c) {
      result_v[0][c] = phi[0][c] * n_l[0][c];
    }
  } else if (wrt_key == mdl_key) {
    for (int c=0; c!=ncells; ++c) {
      result_v[0][c] = phi[0][c] * s_l[0][c];
    }
  } else if (wrt_key == si_key) {
    for (int c=0; c!=ncells; ++c) {
      result_v[0][c] = phi[0][c] * n_i[0][c];
    }
  } else if (wrt_key == mdi_key) {
    for (int c=0; c!=ncells; ++c) {
      result_v[0][c] = phi[0][c] * s_i[0][c];
    }
  } else if (wrt_key == sg_key) {
    for (int c=0; c!=ncells; ++c) {
      result_v[0][c] = phi[0][c] * n_g[0][c]*omega_g[0][c];
    }
  } else if (wrt_key == mdg_key) {
    for (int c=0; c!=ncells; ++c) {
      result_v[0][c] = phi[0][c] * s_g[0][c]*omega_g[0][c];
    }
  } else if (wrt_key == mfg_key) {
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
} //namespace
