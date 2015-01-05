/*
  This is the energy component of the ATS and Amanzi codes. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)

  FieldEvaluator for water content.

  Wrapping this conserved quantity as a field evaluator makes it easier to take
  derivatives, keep updated, and the like.  The equation for this is simply:

    WC = phi * (s_liquid * n_liquid + omega_gas * s_gas * n_gas)

  This is simply the conserved quantity in Richards equation.
*/

#include "twophase_energy_evaluator.hh"

namespace Amanzi {
namespace Energy {

/* ******************************************************************
* Constructor from ParameterList
****************************************************************** */
TwoPhaseEnergyEvaluator::TwoPhaseEnergyEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist) {
  my_key_ = plist_.get<std::string>("energy key", "energy");

  dependencies_.insert(std::string("porosity"));

  dependencies_.insert(std::string("saturation_liquid"));
  dependencies_.insert(std::string("molar_density_liquid"));
  dependencies_.insert(std::string("internal_energy_liquid"));

  // dependencies_.insert(std::string("saturation_gas"));
  dependencies_.insert(std::string("molar_density_gas"));
  dependencies_.insert(std::string("internal_energy_gas"));

  dependencies_.insert(std::string("internal_energy_rock"));
  dependencies_.insert(std::string("density_rock"));
};


/* ******************************************************************
* Copy constructor.
****************************************************************** */
TwoPhaseEnergyEvaluator::TwoPhaseEnergyEvaluator(const TwoPhaseEnergyEvaluator& other) :
    SecondaryVariableFieldEvaluator(other) {};


/* ******************************************************************
* Copy constructor.
****************************************************************** */
Teuchos::RCP<FieldEvaluator> TwoPhaseEnergyEvaluator::Clone() const
{
  return Teuchos::rcp(new TwoPhaseEnergyEvaluator(*this));
};


/* ******************************************************************
* Field evaluator.
****************************************************************** */
void TwoPhaseEnergyEvaluator::EvaluateField_(
    const Teuchos::Ptr<State>& S,
    const Teuchos::Ptr<CompositeVector>& result) 
{
  const Epetra_MultiVector& s_l = *S->GetFieldData("saturation_liquid")->ViewComponent("cell",false);
  const Epetra_MultiVector& n_l = *S->GetFieldData("molar_density_liquid")->ViewComponent("cell",false);
  const Epetra_MultiVector& u_l = *S->GetFieldData("internal_energy_liquid")->ViewComponent("cell",false);

  // const Epetra_MultiVector& s_g = *S->GetFieldData("saturation_gas")->ViewComponent("cell",false);
  const Epetra_MultiVector& n_g = *S->GetFieldData("molar_density_gas")->ViewComponent("cell",false);
  const Epetra_MultiVector& u_g = *S->GetFieldData("internal_energy_gas")->ViewComponent("cell",false);

  const Epetra_MultiVector& phi = *S->GetFieldData("porosity")->ViewComponent("cell",false);
  const Epetra_MultiVector& u_rock = *S->GetFieldData("internal_energy_rock")->ViewComponent("cell",false);
  const Epetra_MultiVector& rho_rock = *S->GetFieldData("density_rock")->ViewComponent("cell",false);
  Epetra_MultiVector& result_v = *result->ViewComponent("cell",false);

  Teuchos::RCP<const AmanziMesh::Mesh> mesh = S->GetMesh();

  int ncells = result->size("cell", false);
  for (int c=0; c!=ncells; ++c) {
    double s_g = 1.0 - s_l[0][c];
    result_v[0][c] = phi[0][c] * (
        s_l[0][c] * n_l[0][c] * u_l[0][c]
        + s_g * n_g[0][c] * u_g[0][c])
        + (1.0 - phi[0][c]) * u_rock[0][c] * rho_rock[0][c];
    result_v[0][c] *= mesh->cell_volume(c);
  }
};


/* ******************************************************************
* Field derivative evaluator.
****************************************************************** */
void TwoPhaseEnergyEvaluator::EvaluateFieldPartialDerivative_(
    const Teuchos::Ptr<State>& S,
    Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{
  const Epetra_MultiVector& s_l = *S->GetFieldData("saturation_liquid")->ViewComponent("cell",false);
  const Epetra_MultiVector& n_l = *S->GetFieldData("molar_density_liquid")->ViewComponent("cell",false);
  const Epetra_MultiVector& u_l = *S->GetFieldData("internal_energy_liquid")->ViewComponent("cell",false);

  const Epetra_MultiVector& s_g = *S->GetFieldData("saturation_gas")->ViewComponent("cell",false);
  const Epetra_MultiVector& n_g = *S->GetFieldData("molar_density_gas")->ViewComponent("cell",false);
  const Epetra_MultiVector& u_g = *S->GetFieldData("internal_energy_gas")->ViewComponent("cell",false);

  const Epetra_MultiVector& phi = *S->GetFieldData("porosity")->ViewComponent("cell",false);
  const Epetra_MultiVector& u_rock = *S->GetFieldData("internal_energy_rock")->ViewComponent("cell",false);
  const Epetra_MultiVector& rho_rock = *S->GetFieldData("density_rock")->ViewComponent("cell",false);
  const Epetra_MultiVector& cell_volume = *S->GetFieldData("cell_volume")->ViewComponent("cell",false);

  Epetra_MultiVector& result_v = *result->ViewComponent("cell",false);

  int ncells = result->size("cell",false);
  if (wrt_key == "porosity") {
    for (int c=0; c!=ncells; ++c) {
      result_v[0][c] = s_l[0][c]*n_l[0][c]*u_l[0][c]
          + s_g[0][c]*n_g[0][c]*u_g[0][c]
          - rho_rock[0][c] * u_rock[0][c];
    }

  } else if (wrt_key == "saturation_liquid") {
    for (int c=0; c!=ncells; ++c) {
      result_v[0][c] = phi[0][c]*n_l[0][c]*u_l[0][c];
    }
  } else if (wrt_key == "molar_density_liquid") {
    for (int c=0; c!=ncells; ++c) {
      result_v[0][c] = phi[0][c]*s_l[0][c]*u_l[0][c];
    }
  } else if (wrt_key == "internal_energy_liquid") {
    for (int c=0; c!=ncells; ++c) {
      result_v[0][c] = phi[0][c]*s_l[0][c]*n_l[0][c];
    }

  } else if (wrt_key == "saturation_gas") {
    for (int c=0; c!=ncells; ++c) {
      result_v[0][c] = phi[0][c]*n_g[0][c]*u_g[0][c];
    }
  } else if (wrt_key == "molar_density_gas") {
    for (int c=0; c!=ncells; ++c) {
      result_v[0][c] = phi[0][c]*s_g[0][c]*u_g[0][c];
    }
  } else if (wrt_key == "internal_energy_gas") {
    for (int c=0; c!=ncells; ++c) {
      result_v[0][c] = phi[0][c]*s_g[0][c]*n_g[0][c];
    }

  } else if (wrt_key == "internal_energy_rock") {
    for (int c=0; c!=ncells; ++c) {
      result_v[0][c] = (1.0 - phi[0][c])*rho_rock[0][c];
    }
  } else if (wrt_key == "density_rock") {
    for (int c=0; c!=ncells; ++c) {
      result_v[0][c] = (1.0 - phi[0][c])*u_rock[0][c];
    }
  } else {
    ASSERT(0);
  }

  for (int c=0; c!=ncells; ++c) {
    result_v[0][c] *= cell_volume[0][c];
  }
};

}  // namespace Energy
}  // namespace Amanzi
