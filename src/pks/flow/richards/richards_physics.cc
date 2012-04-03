/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  A base two-phase, thermal Richard's equation with water vapor.

  License: BSD
  Authors: Neil Carlson (version 1)
  Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
  Ethan Coon (ATS version) (ecoon@lanl.gov)
*/

#include "richards.hh"

namespace Amanzi {
namespace AmanziFlow {

void Richards::ApplyDiffusion_(const Teuchos::RCP<State>& S,
        const Teuchos::RCP<CompositeVector>& g) {
};

void Richards::AddAccumulation_(const Teuchos::RCP<CompositeVector>& g) {
};

/* ******************************************************************
 * Update secondary variables, calculated in various methods below.
 ****************************************************************** */
void Richards::UpdateSecondaryVariables_(const Teuchos::RCP<State>& S) {
  // get temp, pressure
  Teuchos::RCP<const CompositeVector> temp = S->GetFieldData("temperature");
  Teuchos::RCP<const CompositeVector> pres = S->GetFieldData("pressure");
  Teuchos::RCP<const double> p_atm = S->GetScalarData("atmospheric_pressure")

    Teuchos::RCP<CompositeVector> mol_frac_gas = S->GetFieldData("mol_frac_gas", "flow");
  Teuchos::RCP<CompositeVector> dens_gas = S->GetFieldData("density_gas", "flow");
  Teuchos::RCP<CompositeVector> mol_dens_gas = S->GetFieldData("molar_density_gas", "flow");
  Teuchos::RCP<CompositeVector> visc_liq = S->GetFieldData("viscosity_liquid", "flow");
  Teuchos::RCP<CompositeVector> dens_liq = S->GetFieldData("density_liquid", "flow");
  Teuchos::RCP<CompositeVector> mol_dens_liq = S->GetFieldData("molar_density_liquid", "flow");
  Teuchos::RCP<CompositeVector> sat_gas = S->GetFieldData("saturation_gas", "flow");
  Teuchos::RCP<CompositeVector> sat_liq = S->GetFieldData("saturation_liquid", "flow");
  Teuchos::RCP<CompositeVector> rel_perm = S->GetFieldData("relative_permeability", "flow");

  // calculate densities, viscosities using EOS
  DensityLiquid_(S, *temp, *pres, dens_liq, mol_dens_liq);
  ViscosityLiquid_(S, *temp, visc_liq);
  DensityGas_(S, *temp, *pres, *p_atm, mol_frac_gas, dens_gas, mol_dens_gas);

  // calculate saturations using WRM
  Saturation_(S, *pres, *p_atm, sat_liq);
  sat_gas->PutScalar(1.0);
  sat_gas->Update(-1.0, *sat_liq, 1.0);

  // calculate rel perm using WRM
  RelativePermeability_(S, *sat_liq, rel_perm);
};

void Richards::DensityLiquid_(const Teuchos::RCP<State>& S,
        const CompositeVector& temp, const CompositeVector& pres,
        const Teuchos::RCP<CompositeVector>& dens_liq,
        const Teuchos::RCP<CompositeVector>& mol_dens_liq) {

  int c_owned = S->mesh()->count_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c=0; c!=c_owned; ++c) {
    (*dens_liq)("cell",0,c) = eos_liquid_->MassDensity(temp("cell",0,c), pres("cell",0,c));
    (*mol_dens_liq)("cell",0,c) = eos_liquid_->
      MolarDensityFromMassDensity((*dens_liq)("cell",0,c));
  }
};

void Richards::DensityGas_(const Teuchos::RCP<State>& S,
                           const CompositeVector& temp,
                           const CompositeVector& pres, const double& p_atm,
                           const Teuchos::RCP<CompositeVector>& mol_frac_gas,
                           const Teuchos::RCP<CompositeVector>& dens_gas,
                           const Teuchos::RCP<CompositeVector>& mol_dens_gas) {

  int c_owned = S->mesh()->count_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c=0; c!=c_owned; ++c) {
    double p_sat = eos_gas_->SaturatedVaporPressure(temp("cell",0,c));
    (*mol_frac_gas)("cell",0,c) = p_sat/p_atm;
    (*mol_dens_liq)("cell",0,c) = eos_gas_->
      MolarDensity(temp("cell",0,c), pres("cell",0,c));
    (*dens_liq)("cell",0,c) = eos_gas_->
      MassDensityFromMolarDensity((*mol_dens_liq)("cell",0,c), (*mol_frac_gas)("cell",0,c));
  }
};

void Richards::ViscosityLiquid_(const Teuchos::RCP<State>& S,
        const CompositeVector& temp,
        const Teuchos::RCP<CompositeVector>& visc_liq) {
  int c_owned = S->mesh()->count_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c=0; c!=c_owned; ++c) {
    (*visc_liq)("cell",0,c) = eos_liquid_->Viscosity(temp("cell",0,c));
  }
};

void Richards::Saturation_(const Teuchos::RCP<State>& S,
                           const CompositeVector& pres, const double& p_atm,
                           const Teuchos::RCP<CompositeVector>& sat_liq) {
  int c_owned = S->mesh()->count_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c=0; c!=c_owned; ++c) {
    (*sat_liq)("cell",0,c) = wrm_->saturation(p_atm - pres("cell",0,c));
  }
};

void Richards::RelativePermeability_(const Teuchos::RCP<State>& S,
        const CompositeVector& pres, const double& p_atm,
        const Teuchos::RCP<CompositeVector>& rel_perm) {
  int c_owned = S->mesh()->count_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c=0; c!=c_owned; ++c) {
    (*rel_perm)("cell",0,c) = wrm_->k_relative(p_atm - pres("cell",0,c));
  }
};

} //namespace
} //namespace
