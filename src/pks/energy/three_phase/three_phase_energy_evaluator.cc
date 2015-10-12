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


#include "three_phase_energy_evaluator.hh"

namespace Amanzi {
namespace Energy {

ThreePhaseEnergyEvaluator::ThreePhaseEnergyEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist) {
  //  my_key_ = plist_.get<std::string>("energy key", "energy");
  // my_key_ = plist_.get<std::string>("saturation liquid key", "saturation_liquid");

 Key domain_name = getDomain(my_key_);

 por_key_ = getKey(domain_name, "porosity");
 bpor_key_= getKey(domain_name, "base_porosity");
 sl_key_ = getKey(domain_name, "saturation_liquid");
 mdl_key_ = getKey(domain_name, "molar_density_liquid");
 iel_key_= getKey(domain_name, "internal_energy_liquid");
 sg_key_ = getKey(domain_name, "saturation_gas");
 mdg_key_ = getKey(domain_name, "molar_density_gas");
 ieg_key_ = getKey(domain_name, "internal_energy_gas");
 si_key_ = getKey(domain_name, "saturation_ice");
 mdi_key_ = getKey(domain_name, "molar_density_ice");
 iei_key_ = getKey(domain_name, "internal_energy_ice");
 dr_key_ = getKey(domain_name, "density_rock");
 ier_key_ = getKey(domain_name, "internal_energy_rock");
 cv_key_ = getKey(domain_name, "cell_volume");

 dependencies_.insert(std::string(por_key_));
 dependencies_.insert(std::string(bpor_key_));

 dependencies_.insert(std::string(sl_key_));
 dependencies_.insert(std::string(mdl_key_));
 dependencies_.insert(std::string(iel_key_));
 
 dependencies_.insert(std::string(sg_key_));
 dependencies_.insert(std::string(mdg_key_));
 dependencies_.insert(std::string(ieg_key_));
 
 dependencies_.insert(std::string(si_key_));
 dependencies_.insert(std::string(mdi_key_));
 dependencies_.insert(std::string(iei_key_));
 
 dependencies_.insert(std::string(dr_key_));
 dependencies_.insert(std::string(ier_key_));
 dependencies_.insert(std::string(cv_key_));
 
  //  check_derivative_ = true;
};

ThreePhaseEnergyEvaluator::ThreePhaseEnergyEvaluator(const ThreePhaseEnergyEvaluator& other) :
    SecondaryVariableFieldEvaluator(other) {};

Teuchos::RCP<FieldEvaluator>
ThreePhaseEnergyEvaluator::Clone() const {
  return Teuchos::rcp(new ThreePhaseEnergyEvaluator(*this));
};


void ThreePhaseEnergyEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result) {
  Key domain_name = getDomain(my_key_);
  por_key_ = getKey(domain_name, "porosity");
  bpor_key_= getKey(domain_name, "base_porosity");
  sl_key_ = getKey(domain_name, "saturation_liquid");
  mdl_key_ = getKey(domain_name, "molar_density_liquid");
  iel_key_= getKey(domain_name, "internal_energy_liquid");
  sg_key_ = getKey(domain_name, "saturation_gas");
  mdg_key_ = getKey(domain_name, "molar_density_gas");
  ieg_key_ = getKey(domain_name, "internal_energy_gas");
  si_key_ = getKey(domain_name, "saturation_ice");
  mdi_key_ = getKey(domain_name, "molar_density_ice");
  iei_key_ = getKey(domain_name, "internal_energy_ice");
  dr_key_ = getKey(domain_name, "density_rock");
  ier_key_ = getKey(domain_name, "internal_energy_rock");
  cv_key_ = getKey(domain_name, "cell_volume");
 
  const Epetra_MultiVector& s_l = *S->GetFieldData(sl_key_)->ViewComponent("cell",false);
  const Epetra_MultiVector& n_l = *S->GetFieldData(mdl_key_)->ViewComponent("cell",false);
  const Epetra_MultiVector& u_l = *S->GetFieldData(iel_key_)->ViewComponent("cell",false);

  const Epetra_MultiVector& s_g = *S->GetFieldData(sg_key_)->ViewComponent("cell",false);
  const Epetra_MultiVector& n_g = *S->GetFieldData(mdg_key_)->ViewComponent("cell",false);
  const Epetra_MultiVector& u_g = *S->GetFieldData(ieg_key_)->ViewComponent("cell",false);

  const Epetra_MultiVector& s_i = *S->GetFieldData(si_key_)->ViewComponent("cell",false);
  const Epetra_MultiVector& n_i = *S->GetFieldData(mdi_key_)->ViewComponent("cell",false);
  const Epetra_MultiVector& u_i = *S->GetFieldData(iei_key_)->ViewComponent("cell",false);

  const Epetra_MultiVector& phi = *S->GetFieldData(por_key_)->ViewComponent("cell",false);
  const Epetra_MultiVector& phib = *S->GetFieldData(bpor_key_)->ViewComponent("cell",false);
  const Epetra_MultiVector& u_rock = *S->GetFieldData(ier_key_)->ViewComponent("cell",false);
  const Epetra_MultiVector& rho_rock = *S->GetFieldData(dr_key_)->ViewComponent("cell",false);
  const Epetra_MultiVector& cell_volume = *S->GetFieldData(cv_key_)->ViewComponent("cell",false);
  Epetra_MultiVector& result_v = *result->ViewComponent("cell",false);


  int ncells = result->size("cell", false);
  for (int c=0; c!=ncells; ++c) {
    result_v[0][c] = phi[0][c] * (
        s_l[0][c] * n_l[0][c] * u_l[0][c]
        + s_i[0][c] * n_i[0][c] * u_i[0][c]
        + s_g[0][c] * n_g[0][c] * u_g[0][c])
        + (1.0 - phib[0][c]) * u_rock[0][c] * rho_rock[0][c];
    result_v[0][c] *= cell_volume[0][c];
  }
};


void ThreePhaseEnergyEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {
  
  Key domain_name = getDomain(my_key_);
  por_key_ = getKey(domain_name, "porosity");
  bpor_key_= getKey(domain_name, "base_porosity");
  sl_key_ = getKey(domain_name, "saturation_liquid");
  mdl_key_ = getKey(domain_name, "molar_density_liquid");
  iel_key_= getKey(domain_name, "internal_energy_liquid");
  sg_key_ = getKey(domain_name, "saturation_gas");
  mdg_key_ = getKey(domain_name, "molar_density_gas");
  ieg_key_ = getKey(domain_name, "internal_energy_gas");
  si_key_ = getKey(domain_name, "saturation_ice");
  mdi_key_ = getKey(domain_name, "molar_density_ice");
  iei_key_ = getKey(domain_name, "internal_energy_ice");
  dr_key_ = getKey(domain_name, "density_rock");
  ier_key_ = getKey(domain_name, "internal_energy_rock");
  cv_key_ = getKey(domain_name, "cell_volume");
  const Epetra_MultiVector& s_l = *S->GetFieldData(sl_key_)->ViewComponent("cell",false);
  const Epetra_MultiVector& n_l = *S->GetFieldData(mdl_key_)->ViewComponent("cell",false);
  const Epetra_MultiVector& u_l = *S->GetFieldData(iel_key_)->ViewComponent("cell",false);

  const Epetra_MultiVector& s_g = *S->GetFieldData(sg_key_)->ViewComponent("cell",false);
  const Epetra_MultiVector& n_g = *S->GetFieldData(mdg_key_)->ViewComponent("cell",false);
  const Epetra_MultiVector& u_g = *S->GetFieldData(ieg_key_)->ViewComponent("cell",false);

  const Epetra_MultiVector& s_i = *S->GetFieldData(si_key_)->ViewComponent("cell",false);
  const Epetra_MultiVector& n_i = *S->GetFieldData(mdi_key_)->ViewComponent("cell",false);
  const Epetra_MultiVector& u_i = *S->GetFieldData(iei_key_)->ViewComponent("cell",false);

  const Epetra_MultiVector& phi = *S->GetFieldData(por_key_)->ViewComponent("cell",false);
  const Epetra_MultiVector& phib = *S->GetFieldData(bpor_key_)->ViewComponent("cell",false);
  const Epetra_MultiVector& u_rock = *S->GetFieldData(ier_key_)->ViewComponent("cell",false);
  const Epetra_MultiVector& rho_rock = *S->GetFieldData(dr_key_)->ViewComponent("cell",false);
  const Epetra_MultiVector& cell_volume = *S->GetFieldData(cv_key_)->ViewComponent("cell",false);
  Epetra_MultiVector& result_v = *result->ViewComponent("cell",false);

  if (wrt_key == por_key_) {
    for (unsigned int c=0; c!=result->size("cell"); ++c) {
      result_v[0][c] = s_l[0][c]*n_l[0][c]*u_l[0][c]
          + s_i[0][c]*n_i[0][c]*u_i[0][c]
          + s_g[0][c]*n_g[0][c]*u_g[0][c];
    }
  } else if (wrt_key == bpor_key_) {
    for (unsigned int c=0; c!=result->size("cell"); ++c) {
      result_v[0][c] = - rho_rock[0][c] * u_rock[0][c];
    }

  } else if (wrt_key == sl_key_) {
    for (unsigned int c=0; c!=result->size("cell"); ++c) {
      result_v[0][c] = phi[0][c]*n_l[0][c]*u_l[0][c];
    }
  } else if (wrt_key == mdl_key_) {
    for (unsigned int c=0; c!=result->size("cell"); ++c) {
      result_v[0][c] = phi[0][c]*s_l[0][c]*u_l[0][c];
    }
  } else if (wrt_key == iel_key_) {
    for (unsigned int c=0; c!=result->size("cell"); ++c) {
      result_v[0][c] = phi[0][c]*s_l[0][c]*n_l[0][c];
    }

  } else if (wrt_key == sg_key_) {
    for (unsigned int c=0; c!=result->size("cell"); ++c) {
      result_v[0][c] = phi[0][c]*n_g[0][c]*u_g[0][c];
    }
  } else if (wrt_key == mdg_key_) {
    for (unsigned int c=0; c!=result->size("cell"); ++c) {
      result_v[0][c] = phi[0][c]*s_g[0][c]*u_g[0][c];
    }
  } else if (wrt_key == ieg_key_) {
    for (unsigned int c=0; c!=result->size("cell"); ++c) {
      result_v[0][c] = phi[0][c]*s_g[0][c]*n_g[0][c];
    }

  } else if (wrt_key == si_key_) {
    for (unsigned int c=0; c!=result->size("cell"); ++c) {
      result_v[0][c] = phi[0][c]*n_i[0][c]*u_i[0][c];
    }
  } else if (wrt_key == mdi_key_) {
    for (unsigned int c=0; c!=result->size("cell"); ++c) {
      result_v[0][c] = phi[0][c]*s_i[0][c]*u_i[0][c];
    }
  } else if (wrt_key == iei_key_) {
    for (unsigned int c=0; c!=result->size("cell"); ++c) {
      result_v[0][c] = phi[0][c]*s_i[0][c]*n_i[0][c];
    }

  } else if (wrt_key == ier_key_) {
    for (unsigned int c=0; c!=result->size("cell"); ++c) {
      result_v[0][c] = (1.0 - phib[0][c])*rho_rock[0][c];
    }
  } else if (wrt_key == dr_key_) {
    for (unsigned int c=0; c!=result->size("cell"); ++c) {
      result_v[0][c] = (1.0 - phib[0][c])*u_rock[0][c];
    }
  } else {
    ASSERT(0);
  }

  for (unsigned int c=0; c!=result->size("cell"); ++c) {
    result_v[0][c] *= cell_volume[0][c];
  }
};


} //namespace
} //namespace
