/*
  Energy

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)

  FieldEvaluator for the total internal energy. Wrapping this conserved
  quantity as a field evaluator makes it easier to take derivatives, 
  keep updated, and the like. The equation for this is simply:

    IE = phi * (s_liquid * n_liquid * u_liquid + s_gas * n_gas * u_gas)
       + (1 - phi) * rho_rock * u_rock
*/

#include "TotalEnergyEvaluator.hh"

namespace Amanzi {
namespace Energy {

/* ******************************************************************
* Constructor from ParameterList
****************************************************************** */
TotalEnergyEvaluator::TotalEnergyEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist)
{
  my_key_ = plist_.get<std::string>("energy key", "energy");
  vapor_diffusion_ = plist_.get<bool>("vapor diffusion");

  dependencies_.insert(std::string("porosity"));

  dependencies_.insert(std::string("saturation_liquid"));
  dependencies_.insert(std::string("molar_density_liquid"));
  dependencies_.insert(std::string("internal_energy_liquid"));

  if (vapor_diffusion_) {
    dependencies_.insert(std::string("molar_density_gas"));
    dependencies_.insert(std::string("internal_energy_gas"));
  }

  dependencies_.insert(std::string("internal_energy_rock"));
  dependencies_.insert(std::string("particle_density"));
}


/* ******************************************************************
* Copy constructor.
****************************************************************** */
TotalEnergyEvaluator::TotalEnergyEvaluator(const TotalEnergyEvaluator& other) :
    SecondaryVariableFieldEvaluator(other) {};


/* ******************************************************************
* Copy constructor.
****************************************************************** */
Teuchos::RCP<FieldEvaluator> TotalEnergyEvaluator::Clone() const {
  return Teuchos::rcp(new TotalEnergyEvaluator(*this));
}


/* ******************************************************************
* Field evaluator.
****************************************************************** */
void TotalEnergyEvaluator::EvaluateField_(
    const Teuchos::Ptr<State>& S,
    const Teuchos::Ptr<CompositeVector>& result) 
{
  const auto& s_l = *S->GetFieldData("saturation_liquid")->ViewComponent("cell");
  const auto& n_l = *S->GetFieldData("molar_density_liquid")->ViewComponent("cell");
  const auto& u_l = *S->GetFieldData("internal_energy_liquid")->ViewComponent("cell");

  Teuchos::RCP<const Epetra_MultiVector> n_g, u_g;
  if (vapor_diffusion_) {
    n_g = S->GetFieldData("molar_density_gas")->ViewComponent("cell");
    u_g = S->GetFieldData("internal_energy_gas")->ViewComponent("cell");
  } 

  const auto& phi = *S->GetFieldData("porosity")->ViewComponent("cell");
  const auto& u_rock = *S->GetFieldData("internal_energy_rock")->ViewComponent("cell");
  const auto& rho_rock = *S->GetFieldData("particle_density")->ViewComponent("cell");

  auto& result_v = *result->ViewComponent("cell");
  int ncells = result->size("cell");

  for (int c = 0; c != ncells; ++c) {
    result_v[0][c] = phi[0][c] * s_l[0][c] * n_l[0][c] * u_l[0][c]
                   + (1.0 - phi[0][c]) * u_rock[0][c] * rho_rock[0][c];
    if (vapor_diffusion_) {
      double s_g = 1.0 - s_l[0][c];
      result_v[0][c] += phi[0][c] * s_g * (*n_g)[0][c] * (*u_g)[0][c];
    }
  }
}


/* ******************************************************************
* Field derivative evaluator.
****************************************************************** */
void TotalEnergyEvaluator::EvaluateFieldPartialDerivative_(
    const Teuchos::Ptr<State>& S,
    Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{
  const auto& s_l = *S->GetFieldData("saturation_liquid")->ViewComponent("cell");
  const auto& n_l = *S->GetFieldData("molar_density_liquid")->ViewComponent("cell");
  const auto& u_l = *S->GetFieldData("internal_energy_liquid")->ViewComponent("cell");

  Teuchos::RCP<const Epetra_MultiVector> n_g, u_g;
  if (vapor_diffusion_) {
    n_g = S->GetFieldData("molar_density_gas")->ViewComponent("cell");
    u_g = S->GetFieldData("internal_energy_gas")->ViewComponent("cell");
  }

  const auto& phi = *S->GetFieldData("porosity")->ViewComponent("cell");
  const auto& u_rock = *S->GetFieldData("internal_energy_rock")->ViewComponent("cell");
  const auto& rho_rock = *S->GetFieldData("particle_density")->ViewComponent("cell");

  auto& result_v = *result->ViewComponent("cell");
  int ncells = result->size("cell");

  if (wrt_key == "porosity") {
    for (int c = 0; c != ncells; ++c) {
      result_v[0][c] = s_l[0][c] * n_l[0][c] * u_l[0][c]
                     - rho_rock[0][c] * u_rock[0][c];
      if (vapor_diffusion_) {
        double s_g = 1.0 - s_l[0][c];
        result_v[0][c] += s_g * (*n_g)[0][c] * (*u_g)[0][c];
      }
    }
  } else if (wrt_key == "saturation_liquid") {
    for (int c = 0; c != ncells; ++c) {
      result_v[0][c] = phi[0][c] * n_l[0][c] * u_l[0][c];
      if (vapor_diffusion_) {
        result_v[0][c] -= (*n_g)[0][c] * (*u_g)[0][c];
      }
    }
  } else if (wrt_key == "molar_density_liquid") {
    for (int c = 0; c != ncells; ++c) {
      result_v[0][c] = phi[0][c] * s_l[0][c] * u_l[0][c];
    }
  } else if (wrt_key == "internal_energy_liquid") {
    for (int c = 0; c != ncells; ++c) {
      result_v[0][c] = phi[0][c] * s_l[0][c] * n_l[0][c];
    }

  } else if (wrt_key == "molar_density_gas") {
    for (int c = 0; c != ncells; ++c) {
      double s_g = 1.0 - s_l[0][c];
      result_v[0][c] = phi[0][c] * s_g * (*u_g)[0][c];
    }
  } else if (wrt_key == "internal_energy_gas") {
    for (int c = 0; c != ncells; ++c) {
      double s_g = 1.0 - s_l[0][c];
      result_v[0][c] = phi[0][c] * s_g * (*n_g)[0][c];
    }

  } else if (wrt_key == "internal_energy_rock") {
    for (int c = 0; c != ncells; ++c) {
      result_v[0][c] = (1.0 - phi[0][c]) * rho_rock[0][c];
    }
  } else if (wrt_key == "particle_density") {
    for (int c = 0; c != ncells; ++c) {
      result_v[0][c] = (1.0 - phi[0][c]) * u_rock[0][c];
    }
  } else {
    AMANZI_ASSERT(0);
  }
}

}  // namespace Energy
}  // namespace Amanzi
