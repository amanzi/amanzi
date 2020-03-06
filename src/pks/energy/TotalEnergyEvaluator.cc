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
  my_key_ = plist_.get<std::string>("energy key");
  auto prefix = Keys::getDomainPrefix(my_key_);

  vapor_diffusion_ = plist_.get<bool>("vapor diffusion");
  ie_rock_key_ = plist_.get<std::string>("internal energy rock key");
  ie_liquid_key_ = prefix + "internal_energy_liquid";
  ie_gas_key_ = prefix + "internal_energy_gas";

  mol_density_liquid_key_ = prefix + "molar_density_liquid";
  mol_density_gas_key_ = prefix + "molar_density_gas";

  particle_density_key_ = plist_.get<std::string>("particle density key");
  porosity_key_ = prefix + "porosity";
  sat_liquid_key_ = prefix + "saturation_liquid";

  dependencies_.insert(porosity_key_);
  dependencies_.insert(sat_liquid_key_);
  dependencies_.insert(mol_density_liquid_key_);
  dependencies_.insert(ie_liquid_key_);

  if (vapor_diffusion_) {
    dependencies_.insert(mol_density_gas_key_);
    dependencies_.insert(ie_gas_key_);
  }

  dependencies_.insert(ie_rock_key_);
  dependencies_.insert(particle_density_key_);
}


/* ******************************************************************
* Copy constructor.
****************************************************************** */
TotalEnergyEvaluator::TotalEnergyEvaluator(const TotalEnergyEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    vapor_diffusion_(other.vapor_diffusion_),
    particle_density_key_(other.particle_density_key_),
    porosity_key_(other.porosity_key_),
    sat_liquid_key_(other.sat_liquid_key_),
    ie_rock_key_(other.ie_rock_key_),
    ie_liquid_key_(other.ie_liquid_key_),
    ie_gas_key_(other.ie_gas_key_),
    mol_density_liquid_key_(other.mol_density_liquid_key_),
    mol_density_gas_key_(other.mol_density_gas_key_) {};


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
  const auto& n_l = *S->GetFieldData(mol_density_liquid_key_)->ViewComponent("cell");
  const auto& u_l = *S->GetFieldData(ie_liquid_key_)->ViewComponent("cell");

  Teuchos::RCP<const Epetra_MultiVector> n_g, u_g;
  if (vapor_diffusion_) {
    n_g = S->GetFieldData(mol_density_gas_key_)->ViewComponent("cell");
    u_g = S->GetFieldData(ie_gas_key_)->ViewComponent("cell");
  } 

  const auto& phi = *S->GetFieldData(porosity_key_)->ViewComponent("cell");
  const auto& u_rock = *S->GetFieldData(ie_rock_key_)->ViewComponent("cell");
  const auto& rho_rock = *S->GetFieldData(particle_density_key_)->ViewComponent("cell");

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
  const auto& s_l = *S->GetFieldData(sat_liquid_key_)->ViewComponent("cell");
  const auto& n_l = *S->GetFieldData(mol_density_liquid_key_)->ViewComponent("cell");
  const auto& u_l = *S->GetFieldData(ie_liquid_key_)->ViewComponent("cell");

  Teuchos::RCP<const Epetra_MultiVector> n_g, u_g;
  if (vapor_diffusion_) {
    n_g = S->GetFieldData(mol_density_gas_key_)->ViewComponent("cell");
    u_g = S->GetFieldData(ie_gas_key_)->ViewComponent("cell");
  }

  const auto& phi = *S->GetFieldData(porosity_key_)->ViewComponent("cell");
  const auto& u_rock = *S->GetFieldData(ie_rock_key_)->ViewComponent("cell");
  const auto& rho_rock = *S->GetFieldData(particle_density_key_)->ViewComponent("cell");

  auto& result_v = *result->ViewComponent("cell");
  int ncells = result->size("cell");

  if (wrt_key == porosity_key_) {
    for (int c = 0; c != ncells; ++c) {
      result_v[0][c] = s_l[0][c] * n_l[0][c] * u_l[0][c]
                     - rho_rock[0][c] * u_rock[0][c];
      if (vapor_diffusion_) {
        double s_g = 1.0 - s_l[0][c];
        result_v[0][c] += s_g * (*n_g)[0][c] * (*u_g)[0][c];
      }
    }
  } else if (wrt_key == sat_liquid_key_) {
    for (int c = 0; c != ncells; ++c) {
      result_v[0][c] = phi[0][c] * n_l[0][c] * u_l[0][c];
      if (vapor_diffusion_) {
        result_v[0][c] -= (*n_g)[0][c] * (*u_g)[0][c];
      }
    }
  } else if (wrt_key == mol_density_liquid_key_) {
    for (int c = 0; c != ncells; ++c) {
      result_v[0][c] = phi[0][c] * s_l[0][c] * u_l[0][c];
    }
  } else if (wrt_key == ie_liquid_key_) {
    for (int c = 0; c != ncells; ++c) {
      result_v[0][c] = phi[0][c] * s_l[0][c] * n_l[0][c];
    }

  } else if (wrt_key == mol_density_gas_key_) {
    for (int c = 0; c != ncells; ++c) {
      double s_g = 1.0 - s_l[0][c];
      result_v[0][c] = phi[0][c] * s_g * (*u_g)[0][c];
    }
  } else if (wrt_key == ie_gas_key_) {
    for (int c = 0; c != ncells; ++c) {
      double s_g = 1.0 - s_l[0][c];
      result_v[0][c] = phi[0][c] * s_g * (*n_g)[0][c];
    }

  } else if (wrt_key == ie_rock_key_) {
    for (int c = 0; c != ncells; ++c) {
      result_v[0][c] = (1.0 - phi[0][c]) * rho_rock[0][c];
    }
  } else if (wrt_key == particle_density_key_) {
    for (int c = 0; c != ncells; ++c) {
      result_v[0][c] = (1.0 - phi[0][c]) * u_rock[0][c];
    }
  } else {
    AMANZI_ASSERT(0);
  }
}

}  // namespace Energy
}  // namespace Amanzi
