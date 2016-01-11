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

  Key domain_name = getDomain(my_key_);
  phi_key_ = getKey(domain_name, "porosity");
  sl_key_ = getKey(domain_name,"saturation_liquid");
  mdl_key_ = getKey(domain_name,"molar_density_liquid");
  si_key_ = getKey(domain_name,"saturation_ice");
  mdi_key_ = getKey(domain_name,"molar_density_ice");
  sg_key_ = getKey(domain_name,"saturation_gas");
  mdg_key_ = getKey(domain_name,"molar_density_gas");
  mfg_key_ = getKey(domain_name,"mol_frac_gas");
  cv_key_ = getKey(domain_name,"cell_volume");

  dependencies_.insert(phi_key_);

  dependencies_.insert(sl_key_);
  dependencies_.insert(mdl_key_);

  dependencies_.insert(si_key_);
  dependencies_.insert(mdi_key_);

  dependencies_.insert(sg_key_);
  dependencies_.insert(mdg_key_);
  dependencies_.insert(mfg_key_);
  dependencies_.insert(cv_key_);

  //  check_derivative_ = true;
};


PermafrostWaterContent::PermafrostWaterContent(const PermafrostWaterContent& other) :
    SecondaryVariableFieldEvaluator(other),
    phi_key_(other.phi_key_),
    sl_key_(other.sl_key_),
    mdl_key_(other.mdl_key_),
    sg_key_(other.sg_key_),
    mdg_key_(other.mdg_key_),
    mfg_key_(other.mfg_key_),
    si_key_(other.si_key_),
    mdi_key_(other.mdi_key_),
    cv_key_(other.cv_key_)
{};


Teuchos::RCP<FieldEvaluator>
PermafrostWaterContent::Clone() const {
  return Teuchos::rcp(new PermafrostWaterContent(*this));
};


void PermafrostWaterContent::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result) {

  const Epetra_MultiVector& s_l = *S->GetFieldData(sl_key_)
    ->ViewComponent("cell",false);
  const Epetra_MultiVector& n_l = *S->GetFieldData(mdl_key_)
    ->ViewComponent("cell",false);

  const Epetra_MultiVector& s_i = *S->GetFieldData(si_key_)
    ->ViewComponent("cell",false);
  const Epetra_MultiVector& n_i = *S->GetFieldData(mdi_key_)
    ->ViewComponent("cell",false);

  const Epetra_MultiVector& s_g = *S->GetFieldData(sg_key_)
    ->ViewComponent("cell",false);
  const Epetra_MultiVector& n_g = *S->GetFieldData(mdg_key_)
    ->ViewComponent("cell",false);
  const Epetra_MultiVector& omega_g = *S->GetFieldData(mfg_key_)
    ->ViewComponent("cell",false);

  const Epetra_MultiVector& phi = *S->GetFieldData(phi_key_)
    ->ViewComponent("cell",false);
  const Epetra_MultiVector& cell_volume = *S->GetFieldData(cv_key_)
    ->ViewComponent("cell",false);
  Epetra_MultiVector& result_v = *result->ViewComponent("cell",false);

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

  const Epetra_MultiVector& s_l = *S->GetFieldData(sl_key_)
    ->ViewComponent("cell",false);
  const Epetra_MultiVector& n_l = *S->GetFieldData(mdl_key_)
    ->ViewComponent("cell",false);
  
  const Epetra_MultiVector& s_i = *S->GetFieldData(si_key_)
    ->ViewComponent("cell",false);
  const Epetra_MultiVector& n_i = *S->GetFieldData(mdi_key_)
    ->ViewComponent("cell",false);
  
  const Epetra_MultiVector& s_g = *S->GetFieldData(sg_key_)
    ->ViewComponent("cell",false);
  const Epetra_MultiVector& n_g = *S->GetFieldData(mdg_key_)
    ->ViewComponent("cell",false);
  const Epetra_MultiVector& omega_g = *S->GetFieldData(mfg_key_)
    ->ViewComponent("cell",false);

  const Epetra_MultiVector& phi = *S->GetFieldData(phi_key_)
    ->ViewComponent("cell",false);
  const Epetra_MultiVector& cell_volume = *S->GetFieldData(cv_key_)
    ->ViewComponent("cell",false);
  Epetra_MultiVector& result_v = *result->ViewComponent("cell",false);

  int ncells = result->size("cell",false);
  if (wrt_key == phi_key_) {
    for (int c=0; c!=ncells; ++c) {
      result_v[0][c] = s_l[0][c]*n_l[0][c]
          + s_i[0][c]*n_i[0][c]
          + s_g[0][c]*n_g[0][c]*omega_g[0][c];
    }
  } else if (wrt_key == sl_key_) {
    for (int c=0; c!=ncells; ++c) {
      result_v[0][c] = phi[0][c] * n_l[0][c];
    }
  } else if (wrt_key == mdl_key_) {
    for (int c=0; c!=ncells; ++c) {
      result_v[0][c] = phi[0][c] * s_l[0][c];
    }
  } else if (wrt_key == si_key_) {
    for (int c=0; c!=ncells; ++c) {
      result_v[0][c] = phi[0][c] * n_i[0][c];
    }
  } else if (wrt_key == mdi_key_) {
    for (int c=0; c!=ncells; ++c) {
      result_v[0][c] = phi[0][c] * s_i[0][c];
    }
  } else if (wrt_key == sg_key_) {
    for (int c=0; c!=ncells; ++c) {
      result_v[0][c] = phi[0][c] * n_g[0][c]*omega_g[0][c];
    }
  } else if (wrt_key == mdg_key_) {
    for (int c=0; c!=ncells; ++c) {
      result_v[0][c] = phi[0][c] * s_g[0][c]*omega_g[0][c];
    }
  } else if (wrt_key == mfg_key_) {
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
