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

 dependencies_.insert(std::string(getKey(domain_name, "porosity")));
 dependencies_.insert(std::string(getKey(domain_name, "base_porosity")));

 dependencies_.insert(std::string(getKey(domain_name, "saturation_liquid")));
 dependencies_.insert(std::string(getKey(domain_name, "molar_density_liquid")));
 dependencies_.insert(std::string(getKey(domain_name, "internal_energy_liquid")));

 dependencies_.insert(std::string(getKey(domain_name, "saturation_gas")));
 dependencies_.insert(std::string(getKey(domain_name, "molar_density_gas")));
 dependencies_.insert(std::string(getKey(domain_name, "internal_energy_gas")));

 dependencies_.insert(std::string(getKey(domain_name, "saturation_ice")));
 dependencies_.insert(std::string(getKey(domain_name, "molar_density_ice")));
 dependencies_.insert(std::string(getKey(domain_name, "internal_energy_ice")));

 dependencies_.insert(std::string(getKey(domain_name, "density_rock")));
 dependencies_.insert(std::string(getKey(domain_name, "internal_energy_rock")));
  //  dependencies_.insert(std::string("cell_volume"));

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
  
  const Epetra_MultiVector& s_l = *S->GetFieldData("saturation_liquid")->ViewComponent("cell",false);
  const Epetra_MultiVector& n_l = *S->GetFieldData("molar_density_liquid")->ViewComponent("cell",false);
  const Epetra_MultiVector& u_l = *S->GetFieldData("internal_energy_liquid")->ViewComponent("cell",false);

  const Epetra_MultiVector& s_g = *S->GetFieldData("saturation_gas")->ViewComponent("cell",false);
  const Epetra_MultiVector& n_g = *S->GetFieldData("molar_density_gas")->ViewComponent("cell",false);
  const Epetra_MultiVector& u_g = *S->GetFieldData("internal_energy_gas")->ViewComponent("cell",false);

  const Epetra_MultiVector& s_i = *S->GetFieldData("saturation_ice")->ViewComponent("cell",false);
  const Epetra_MultiVector& n_i = *S->GetFieldData("molar_density_ice")->ViewComponent("cell",false);
  const Epetra_MultiVector& u_i = *S->GetFieldData("internal_energy_ice")->ViewComponent("cell",false);

  const Epetra_MultiVector& phi = *S->GetFieldData("porosity")->ViewComponent("cell",false);
  const Epetra_MultiVector& phib = *S->GetFieldData("base_porosity")->ViewComponent("cell",false);
  const Epetra_MultiVector& u_rock = *S->GetFieldData("internal_energy_rock")->ViewComponent("cell",false);
  const Epetra_MultiVector& rho_rock = *S->GetFieldData("density_rock")->ViewComponent("cell",false);
  const Epetra_MultiVector& cell_volume = *S->GetFieldData("cell_volume")->ViewComponent("cell",false);
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
  const Epetra_MultiVector& s_l = *S->GetFieldData("saturation_liquid")->ViewComponent("cell",false);
  const Epetra_MultiVector& n_l = *S->GetFieldData("molar_density_liquid")->ViewComponent("cell",false);
  const Epetra_MultiVector& u_l = *S->GetFieldData("internal_energy_liquid")->ViewComponent("cell",false);

  const Epetra_MultiVector& s_g = *S->GetFieldData("saturation_gas")->ViewComponent("cell",false);
  const Epetra_MultiVector& n_g = *S->GetFieldData("molar_density_gas")->ViewComponent("cell",false);
  const Epetra_MultiVector& u_g = *S->GetFieldData("internal_energy_gas")->ViewComponent("cell",false);

  const Epetra_MultiVector& s_i = *S->GetFieldData("saturation_ice")->ViewComponent("cell",false);
  const Epetra_MultiVector& n_i = *S->GetFieldData("molar_density_ice")->ViewComponent("cell",false);
  const Epetra_MultiVector& u_i = *S->GetFieldData("internal_energy_ice")->ViewComponent("cell",false);

  const Epetra_MultiVector& phi = *S->GetFieldData("porosity")->ViewComponent("cell",false);
  const Epetra_MultiVector& phib = *S->GetFieldData("base_porosity")->ViewComponent("cell",false);
  const Epetra_MultiVector& u_rock = *S->GetFieldData("internal_energy_rock")->ViewComponent("cell",false);
  const Epetra_MultiVector& rho_rock = *S->GetFieldData("density_rock")->ViewComponent("cell",false);
  const Epetra_MultiVector& cell_volume = *S->GetFieldData("cell_volume")->ViewComponent("cell",false);
  Epetra_MultiVector& result_v = *result->ViewComponent("cell",false);

  if (wrt_key == "porosity") {
    for (unsigned int c=0; c!=result->size("cell"); ++c) {
      result_v[0][c] = s_l[0][c]*n_l[0][c]*u_l[0][c]
          + s_i[0][c]*n_i[0][c]*u_i[0][c]
          + s_g[0][c]*n_g[0][c]*u_g[0][c];
    }
  } else if (wrt_key == "base_porosity") {
    for (unsigned int c=0; c!=result->size("cell"); ++c) {
      result_v[0][c] = - rho_rock[0][c] * u_rock[0][c];
    }

  } else if (wrt_key == "saturation_liquid") {
    for (unsigned int c=0; c!=result->size("cell"); ++c) {
      result_v[0][c] = phi[0][c]*n_l[0][c]*u_l[0][c];
    }
  } else if (wrt_key == "molar_density_liquid") {
    for (unsigned int c=0; c!=result->size("cell"); ++c) {
      result_v[0][c] = phi[0][c]*s_l[0][c]*u_l[0][c];
    }
  } else if (wrt_key == "internal_energy_liquid") {
    for (unsigned int c=0; c!=result->size("cell"); ++c) {
      result_v[0][c] = phi[0][c]*s_l[0][c]*n_l[0][c];
    }

  } else if (wrt_key == "saturation_gas") {
    for (unsigned int c=0; c!=result->size("cell"); ++c) {
      result_v[0][c] = phi[0][c]*n_g[0][c]*u_g[0][c];
    }
  } else if (wrt_key == "molar_density_gas") {
    for (unsigned int c=0; c!=result->size("cell"); ++c) {
      result_v[0][c] = phi[0][c]*s_g[0][c]*u_g[0][c];
    }
  } else if (wrt_key == "internal_energy_gas") {
    for (unsigned int c=0; c!=result->size("cell"); ++c) {
      result_v[0][c] = phi[0][c]*s_g[0][c]*n_g[0][c];
    }

  } else if (wrt_key == "saturation_ice") {
    for (unsigned int c=0; c!=result->size("cell"); ++c) {
      result_v[0][c] = phi[0][c]*n_i[0][c]*u_i[0][c];
    }
  } else if (wrt_key == "molar_density_ice") {
    for (unsigned int c=0; c!=result->size("cell"); ++c) {
      result_v[0][c] = phi[0][c]*s_i[0][c]*u_i[0][c];
    }
  } else if (wrt_key == "internal_energy_ice") {
    for (unsigned int c=0; c!=result->size("cell"); ++c) {
      result_v[0][c] = phi[0][c]*s_i[0][c]*n_i[0][c];
    }

  } else if (wrt_key == "internal_energy_rock") {
    for (unsigned int c=0; c!=result->size("cell"); ++c) {
      result_v[0][c] = (1.0 - phib[0][c])*rho_rock[0][c];
    }
  } else if (wrt_key == "density_rock") {
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
