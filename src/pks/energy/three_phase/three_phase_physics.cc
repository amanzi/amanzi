/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Solves:

du/dt + v dot grad h = div Ke grad T .......... fix me!
------------------------------------------------------------------------- */

#include "three_phase.hh"

namespace Amanzi {
namespace Energy {

// -------------------------------------------------------------
// Accumulation of internal energy term du/dt
// -------------------------------------------------------------
void ThreePhase::AddAccumulation_(Teuchos::RCP<CompositeVector> g) {
  Teuchos::RCP<const CompositeVector> poro0 =
    S_inter_->GetFieldData("porosity");
  Teuchos::RCP<const CompositeVector> poro1 =
    S_next_->GetFieldData("porosity");

  Teuchos::RCP<const CompositeVector> density_gas0;
  Teuchos::RCP<const CompositeVector> density_gas1;
  if (iem_gas_->IsMolarBasis()) {
    density_gas0 = S_inter_->GetFieldData("molar_density_gas");
    density_gas1 = S_next_->GetFieldData("molar_density_gas");
  } else {
    density_gas0 = S_inter_->GetFieldData("density_gas");
    density_gas1 = S_next_->GetFieldData("density_gas");
  }

  Teuchos::RCP<const CompositeVector> density_liq0;
  Teuchos::RCP<const CompositeVector> density_liq1;
  if (iem_liquid_->IsMolarBasis()) {
    density_liq0 = S_inter_->GetFieldData("molar_density_liquid");
    density_liq1 = S_next_->GetFieldData("molar_density_liquid");
  } else {
    density_liq0 = S_inter_->GetFieldData("density_liquid");
    density_liq1 = S_next_->GetFieldData("density_liquid");
  }

  Teuchos::RCP<const CompositeVector> density_ice0;
  Teuchos::RCP<const CompositeVector> density_ice1;
  if (iem_ice_->IsMolarBasis()) {
    density_ice0 = S_inter_->GetFieldData("molar_density_ice");
    density_ice1 = S_next_->GetFieldData("molar_density_ice");
  } else {
    density_ice0 = S_inter_->GetFieldData("density_ice");
    density_ice1 = S_next_->GetFieldData("density_ice");
  }

  Teuchos::RCP<const CompositeVector> sat_gas0 =
    S_inter_->GetFieldData("saturation_gas");
  Teuchos::RCP<const CompositeVector> sat_gas1 =
    S_next_->GetFieldData("saturation_gas");

  Teuchos::RCP<const CompositeVector> sat_liq0 =
    S_inter_->GetFieldData("saturation_liquid");
  Teuchos::RCP<const CompositeVector> sat_liq1 =
    S_next_->GetFieldData("saturation_liquid");

  Teuchos::RCP<const CompositeVector> sat_ice0 =
    S_inter_->GetFieldData("saturation_ice");
  Teuchos::RCP<const CompositeVector> sat_ice1 =
    S_next_->GetFieldData("saturation_ice");

  Teuchos::RCP<const CompositeVector> int_energy_gas0 =
    S_inter_->GetFieldData("internal_energy_gas");
  Teuchos::RCP<const CompositeVector> int_energy_gas1 =
    S_next_->GetFieldData("internal_energy_gas");

  Teuchos::RCP<const CompositeVector> int_energy_liq0 =
    S_inter_->GetFieldData("internal_energy_liquid");
  Teuchos::RCP<const CompositeVector> int_energy_liq1 =
    S_next_->GetFieldData("internal_energy_liquid");

  Teuchos::RCP<const CompositeVector> int_energy_ice0 =
    S_inter_->GetFieldData("internal_energy_ice");
  Teuchos::RCP<const CompositeVector> int_energy_ice1 =
    S_next_->GetFieldData("internal_energy_ice");

  Teuchos::RCP<const CompositeVector> int_energy_rock0 =
    S_inter_->GetFieldData("internal_energy_rock");
  Teuchos::RCP<const CompositeVector> int_energy_rock1 =
    S_next_->GetFieldData("internal_energy_rock");

  Teuchos::RCP<const CompositeVector> cell_volume0 =
    S_inter_->GetFieldData("cell_volume");
  Teuchos::RCP<const CompositeVector> cell_volume1 =
    S_next_->GetFieldData("cell_volume");

  Teuchos::RCP<const double> density_rock =
    S_next_->GetScalarData("density_rock");

  double dt = S_next_->time() - S_inter_->time();

  // NOTE: gas, liquid, and ice are done in a ?? basis, but rock is done in a mass basis

  int c_owned = g->size("cell");
  for (int c=0; c != c_owned; ++c) {
    // calculate the energy density at the old and new times
    double edens_gas1 = (*density_gas1)("cell",c) * (*sat_gas1)("cell",c) *
      (*int_energy_gas1)("cell",c);
    double edens_liq1 = (*density_liq1)("cell",c) * (*sat_liq1)("cell",c) *
      (*int_energy_liq1)("cell",c);
    double edens_ice1 = (*density_ice1)("cell",c) * (*sat_ice1)("cell",c) *
      (*int_energy_ice1)("cell",c);
    double edens_rock1 = (*density_rock) * (*int_energy_rock1)("cell",c);
    double energy1 = ((*poro1)("cell",c) *
                      (edens_gas1 + edens_liq1 + edens_ice1) +
                      (1-(*poro1)("cell",c)) * (edens_rock1)) *
                     (*cell_volume1)("cell",c);

    double edens_gas0 = (*density_gas0)("cell",c) * (*sat_gas0)("cell",c) *
      (*int_energy_gas0)("cell",c);
    double edens_liq0 = (*density_liq0)("cell",c) * (*sat_liq0)("cell",c) *
      (*int_energy_liq0)("cell",c);
    double edens_ice0 = (*density_ice0)("cell",c) * (*sat_ice0)("cell",c) *
      (*int_energy_ice0)("cell",c);
    double edens_rock0 = (*density_rock) * (*int_energy_rock0)("cell",c);
    double energy0 = ((*poro0)("cell",c) *
                      (edens_gas0 + edens_liq0 + edens_ice0) +
                      (1-(*poro0)("cell",c)) * (edens_rock0)) *
                     (*cell_volume0)("cell",c);

    // add the time derivative of energy density to the residual
    (*g)("cell",c) += (energy1 - energy0)/dt;
  }
};


// -------------------------------------------------------------
// Advective term for transport of enthalpy, v dot grad h.
// -------------------------------------------------------------
void ThreePhase::AddAdvection_(const Teuchos::RCP<State> S,
        const Teuchos::RCP<CompositeVector> g, bool negate) {

  Teuchos::RCP<CompositeVector> field = advection_->field();
  field->PutScalar(0);

  // set the flux field as the darcy flux
  // NOTE: darcy_flux is a MOLAR flux by choice of the richards flow pk, i.e.
  // [flux] =  mol/(m^2*s)
  Teuchos::RCP<const CompositeVector> darcy_flux = S->GetFieldData("darcy_flux");
  advection_->set_flux(darcy_flux);

  // put the advected quantity in cells
  UpdateEnthalpyLiquid_(S);
  Teuchos::RCP<const CompositeVector> enthalpy_liq = S->GetFieldData("enthalpy_liquid");
  Teuchos::RCP<const CompositeVector> dens_liq = S->GetFieldData("density_liquid");
  Teuchos::RCP<const CompositeVector> n_liq = S->GetFieldData("molar_density_liquid");

  int c_owned = field->size("cell");
  if (iem_liquid_->IsMolarBasis()) {
    // this is clean:
    // if u is in units of J/mol, then the rate of change of energy [J/s] is given by:
    // [flux * u * face_area] = mol/(m^2*s) * J/mol * m^2 = J/s
    for (int c=0; c!=c_owned; ++c) {
      (*field)("cell",c) = (*enthalpy_liq)("cell",c);
    }
  } else {
    // this is ugly because the flux is in a molar basis:
    // if u is in units of J/kg, then the rate of change of energy [J/s] is given by:
    // [flux * rho/n * u * face_area] = mol/(m^2*s) * kg/mol * J/kg * m^2 = J/s
    for (int c=0; c!=c_owned; ++c) {
      (*field)("cell",c) = (*enthalpy_liq)("cell",c)
                        * (*dens_liq)("cell",c) / (*n_liq)("cell",c);
    }
  }

  // put the boundary fluxes in faces

  // Dirichlet data
  // NOTE this boundary flux is in enthalpy, while the Dirichlet BC is temperature.
  // MANY assumptions are coming in to play here... many more than I like.
  // h = n(T,p) * u_l(T) + p_l, and we are using:
  //  - p_l is assumed to be a mimetic discretization, and therefore has p on faces
  //  - u_l is calculated correctly using the boundary data
  //  - n(T,P) is the molar density of the inward _cell_, not calculated from
  //           T,p on the face as it ought to be, as the EOS is currently in
  //           the flow PK.
  Teuchos::RCP<const CompositeVector> pres = S->GetFieldData("pressure");

  for (Functions::BoundaryFunction::Iterator bc = bc_temperature_->begin();
       bc!=bc_temperature_->end(); ++bc) {
    int f = bc->first;
    double T = bc->second;
    double int_energy = iem_liquid_->InternalEnergy(T);

    AmanziMesh::Entity_ID_List cells;
    S->Mesh()->face_get_cells(f, AmanziMesh::OWNED, &cells);
    for (int i=0; i!=cells.size(); ++i) {
      int c = cells[i];
      if (c >= 0) { // only the inward cell is > 0
        double enthalpy = int_energy;
        if (!iem_liquid_->IsMolarBasis()) {
          enthalpy *= (*dens_liq)("cell",c) / (*n_liq)("cell",c);
        }
        enthalpy += (*pres)("face",f)/(*n_liq)("cell",c);
        (*field)("face",f) = enthalpy * fabs((*darcy_flux)("face",f));
      }
    }
  }

  // apply the advection operator and add to residual
  advection_->Apply(bc_flux_);
  if (negate) {
    for (int c=0; c!=c_owned; ++c) {
      (*g)("cell",c) -= (*field)("cell",c);
    }
  } else {
    for (int c=0; c!=c_owned; ++c) {
      (*g)("cell",c) += (*field)("cell",c);
    }
  }
};


// -------------------------------------------------------------
// Diffusion term, div K grad T
// -------------------------------------------------------------
void ThreePhase::ApplyDiffusion_(const Teuchos::RCP<State> S,
          const Teuchos::RCP<CompositeVector> g) {
  // compute the stiffness matrix at the new time
  Teuchos::RCP<const CompositeVector> temp =
    S->GetFieldData("temperature");

  // get conductivity
  UpdateThermalConductivity_(S);
  Teuchos::RCP<CompositeVector> thermal_conductivity =
    S->GetFieldData("thermal_conductivity", "energy");

  // calculate the div-grad operator, apply it to temperature, and add to residual
  matrix_->CreateMFDstiffnessMatrices(*thermal_conductivity);
  matrix_->CreateMFDrhsVectors();
  matrix_->ApplyBoundaryConditions(bc_markers_, bc_values_);
  matrix_->AssembleGlobalMatrices();
  matrix_->ComputeNegativeResidual(*temp, g);
};


// -------------------------------------------------------------
// Update variables, like internal energy, conductivity, etc
//
//    Note: UpdatePhysicalQuantity() methods take only a State as an
//          argument, while PhysicalQuantity() methods take the needed
//          vector quantities as arguments.  This division is on
//          purpose, as future versions of this code may use the
//          Update*() methods as generic "secondary variable PKs",
//          i.e. algebraic PKs.
// -------------------------------------------------------------
void ThreePhase::UpdateSecondaryVariables_(const Teuchos::RCP<State>& S) {
  // update secondary variables
  UpdateInternalEnergyGas_(S);
  UpdateInternalEnergyLiquid_(S);
  UpdateInternalEnergyIce_(S);
  UpdateInternalEnergyRock_(S);
};


// -------------------------------------------------------------
// Update the internal energy of gas in state S.
// -------------------------------------------------------------
void ThreePhase::UpdateInternalEnergyGas_(const Teuchos::RCP<State>& S) {
  Teuchos::RCP<const CompositeVector> temp = S->GetFieldData("temperature");
  Teuchos::RCP<const CompositeVector> mol_frac_gas = S->GetFieldData("mol_frac_gas");
  Teuchos::RCP<CompositeVector> int_energy_gas =
    S->GetFieldData("internal_energy_gas", "energy");

  InternalEnergyGas_(S, *temp, *mol_frac_gas, int_energy_gas);
};


// -------------------------------------------------------------
// Evaluate the internal energy of gas model.
// -------------------------------------------------------------
void ThreePhase::InternalEnergyGas_(const Teuchos::RCP<State>& S,
        const CompositeVector& temp,
        const CompositeVector& mol_frac_gas,
        const Teuchos::RCP<CompositeVector>& int_energy_gas) {
  // just a single model for now -- ignore blocks
  int c_owned = int_energy_gas->size("cell");
  for (int c=0; c != c_owned; ++c) {
    (*int_energy_gas)("cell",c) = iem_gas_->
      InternalEnergy(temp("cell",c), mol_frac_gas("cell",c));
  }
};


// -------------------------------------------------------------
// Update the internal energy of liquid in state S.
// -------------------------------------------------------------
void ThreePhase::UpdateInternalEnergyLiquid_(const Teuchos::RCP<State>& S) {
  Teuchos::RCP<const CompositeVector> temp = S->GetFieldData("temperature");
  Teuchos::RCP<CompositeVector> int_energy_liquid =
    S->GetFieldData("internal_energy_liquid", "energy");

  InternalEnergyLiquid_(S, *temp, int_energy_liquid);
};


// -------------------------------------------------------------
// Evaluate the internal energy of liquid model.
// -------------------------------------------------------------
void ThreePhase::InternalEnergyLiquid_(const Teuchos::RCP<State>& S,
        const CompositeVector& temp,
        const Teuchos::RCP<CompositeVector>& int_energy_liquid) {
  // just a single model for now -- ignore blocks
  int c_owned = int_energy_liquid->size("cell");
  for (int c=0; c != c_owned; ++c) {
    (*int_energy_liquid)("cell",c) = iem_liquid_->InternalEnergy(temp("cell",c));
  }
};


// -------------------------------------------------------------
// Update the internal energy of ice in state S.
// -------------------------------------------------------------
void ThreePhase::UpdateInternalEnergyIce_(const Teuchos::RCP<State>& S) {
  Teuchos::RCP<const CompositeVector> temp = S->GetFieldData("temperature");
  Teuchos::RCP<CompositeVector> int_energy_ice =
    S->GetFieldData("internal_energy_ice", "energy");

  InternalEnergyIce_(S, *temp, int_energy_ice);
};


// -------------------------------------------------------------
// Evaluate the internal energy of ice model.
// -------------------------------------------------------------
void ThreePhase::InternalEnergyIce_(const Teuchos::RCP<State>& S,
        const CompositeVector& temp,
        const Teuchos::RCP<CompositeVector>& int_energy_ice) {
  // just a single model for now -- ignore blocks
  int c_owned = int_energy_ice->size("cell");
  for (int c=0; c != c_owned; ++c) {
    (*int_energy_ice)("cell",c) = iem_ice_->InternalEnergy(temp("cell",c));
  }
};


// -------------------------------------------------------------
// Update the internal energy of rock in state S.
// -------------------------------------------------------------
void ThreePhase::UpdateInternalEnergyRock_(const Teuchos::RCP<State>& S) {
  Teuchos::RCP<const CompositeVector> temp = S->GetFieldData("temperature");
  Teuchos::RCP<CompositeVector> int_energy_rock =
    S->GetFieldData("internal_energy_rock", "energy");

  InternalEnergyRock_(S, *temp, int_energy_rock);
};


// -------------------------------------------------------------
// Evaluate the internal energy of rock model.
// -------------------------------------------------------------
void ThreePhase::InternalEnergyRock_(const Teuchos::RCP<State>& S,
        const CompositeVector& temp,
        const Teuchos::RCP<CompositeVector>& int_energy_rock) {
  // just a single model for now -- ignore blocks
  int c_owned = int_energy_rock->size("cell");
  for (int c=0; c != c_owned; ++c) {
    (*int_energy_rock)("cell",c) = iem_rock_->
      InternalEnergy(temp("cell",c));
  }
};


// -------------------------------------------------------------
// Update the enthalpy of liquid in state S.
// -------------------------------------------------------------
void ThreePhase::UpdateEnthalpyLiquid_(const Teuchos::RCP<State>& S) {
  Teuchos::RCP<const CompositeVector> pres = S->GetFieldData("pressure");

  Teuchos::RCP<const CompositeVector> dens_liq;
  if (iem_liquid_->IsMolarBasis()) {
    dens_liq = S->GetFieldData("molar_density_liquid");
  } else {
    dens_liq = S->GetFieldData("density_liquid");
  }

  Teuchos::RCP<const CompositeVector> int_energy_liquid =
    S->GetFieldData("internal_energy_liquid");

  Teuchos::RCP<CompositeVector> enthalpy_liq =
    S->GetFieldData("enthalpy_liquid", "energy");

  // update enthalpy of liquid
  EnthalpyLiquid_(S, *int_energy_liquid, *pres, *dens_liq, enthalpy_liq);
};


// -------------------------------------------------------------
// Evaluate the enthalpy of liquid, h = u + p/rho
// -------------------------------------------------------------
void ThreePhase::EnthalpyLiquid_(const Teuchos::RCP<State>& S,
        const CompositeVector& int_energy_liquid,
        const CompositeVector& pres, const CompositeVector& dens_liq,
        const Teuchos::RCP<CompositeVector>& enthalpy_liq) {

  // just a single model for now -- ignore blocks
  int c_owned = enthalpy_liq->size("cell");
  for (int c=0; c != c_owned; ++c) {
    (*enthalpy_liq)("cell",c) = int_energy_liquid("cell",c)
                              + pres("cell",c)/dens_liq("cell",c);
  }
};


// -------------------------------------------------------------
// Update the thermal conductivity in state S.
// -------------------------------------------------------------
void ThreePhase::UpdateThermalConductivity_(const Teuchos::RCP<State>& S) {
  Teuchos::RCP<const CompositeVector> poro =
    S->GetFieldData("porosity");
  Teuchos::RCP<const CompositeVector> sat_liq =
    S->GetFieldData("saturation_liquid");
  Teuchos::RCP<const CompositeVector> sat_ice =
    S->GetFieldData("saturation_ice");
  Teuchos::RCP<CompositeVector> thermal_conductivity =
    S->GetFieldData("thermal_conductivity", "energy");

  ThermalConductivity_(S, *poro, *sat_liq, *sat_ice, thermal_conductivity);
};


// -------------------------------------------------------------
// Evaluate the thermal conductivity model.
// -------------------------------------------------------------
void ThreePhase::ThermalConductivity_(const Teuchos::RCP<State>& S,
        const CompositeVector& porosity,
        const CompositeVector& sat_liq,
        const CompositeVector& sat_ice,
        const Teuchos::RCP<CompositeVector>& thermal_conductivity) {

  // just a single model for now -- ignore blocks
  int c_owned = thermal_conductivity->size("cell");
  for (int c=0; c != c_owned; ++c) {
    (*thermal_conductivity)("cell",c) = thermal_conductivity_model_->
      CalculateConductivity(porosity("cell",c), sat_liq("cell",c), sat_ice("cell",c));
  }
};

} //namespace Energy
} //namespace Amanzi
