/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
  A base three-phase, thermal Richard's equation with water vapor and ice.

  License: BSD
  Authors: Ethan Coon (ATS version) (ecoon@lanl.gov)
------------------------------------------------------------------------- */

#include "permafrost.hh"

namespace Amanzi {
namespace Flow {


// -------------------------------------------------------------
// Accumulation of water term du/dt
// -------------------------------------------------------------
void Permafrost::AddAccumulation_(const Teuchos::RCP<CompositeVector>& g) {
  Teuchos::RCP<const CompositeVector> poro0 =
    S_inter_->GetFieldData("porosity");
  Teuchos::RCP<const CompositeVector> poro1 =
    S_next_->GetFieldData("porosity");

  Teuchos::RCP<const CompositeVector> n_gas0 =
    S_inter_->GetFieldData("molar_density_gas");
  Teuchos::RCP<const CompositeVector> n_gas1 =
    S_next_->GetFieldData("molar_density_gas");

  Teuchos::RCP<const CompositeVector> mol_frac_gas0 =
    S_inter_->GetFieldData("mol_frac_gas");
  Teuchos::RCP<const CompositeVector> mol_frac_gas1 =
    S_next_->GetFieldData("mol_frac_gas");

  Teuchos::RCP<const CompositeVector> sat_gas0 =
    S_inter_->GetFieldData("saturation_gas");
  Teuchos::RCP<const CompositeVector> sat_gas1 =
    S_next_->GetFieldData("saturation_gas");

  Teuchos::RCP<const CompositeVector> n_liq0 =
    S_inter_->GetFieldData("molar_density_liquid");
  Teuchos::RCP<const CompositeVector> n_liq1 =
    S_next_->GetFieldData("molar_density_liquid");

  Teuchos::RCP<const CompositeVector> sat_liq0 =
    S_inter_->GetFieldData("saturation_liquid");
  Teuchos::RCP<const CompositeVector> sat_liq1 =
    S_next_->GetFieldData("saturation_liquid");

  Teuchos::RCP<const CompositeVector> n_ice0 =
    S_inter_->GetFieldData("molar_density_ice");
  Teuchos::RCP<const CompositeVector> n_ice1 =
    S_next_->GetFieldData("molar_density_ice");

  Teuchos::RCP<const CompositeVector> sat_ice0 =
    S_inter_->GetFieldData("saturation_ice");
  Teuchos::RCP<const CompositeVector> sat_ice1 =
    S_next_->GetFieldData("saturation_ice");

  Teuchos::RCP<const CompositeVector> cell_volume0 =
    S_inter_->GetFieldData("cell_volume");
  Teuchos::RCP<const CompositeVector> cell_volume1 =
    S_next_->GetFieldData("cell_volume");

  double dt = S_next_->time() - S_inter_->time();


  int c_owned = g->size("cell");
  for (int c=0; c!=c_owned; ++c) {
    // calculate water content of each phase and at each time
    double wc_liq1 = (*n_liq1)("cell",c) * (*sat_liq1)("cell",c);
    double wc_ice1 = (*n_ice1)("cell",c) * (*sat_ice1)("cell",c);
    double wc_gas1 = (*n_gas1)("cell",c) * (*sat_gas1)("cell",c)
                                        * (*mol_frac_gas1)("cell",c);
    double wc1 = (wc_liq1 + wc_gas1 + wc_ice1) * (*poro1)("cell",c)
                                        * (*cell_volume1)("cell",c);

    double wc_liq0 = (*n_liq0)("cell",c) * (*sat_liq0)("cell",c);
    double wc_ice0 = (*n_ice0)("cell",c) * (*sat_ice0)("cell",c);
    double wc_gas0 = (*n_gas0)("cell",c) * (*sat_gas0)("cell",c)
                                        * (*mol_frac_gas0)("cell",c);
    double wc0 = (wc_liq0 + wc_gas0 + wc_ice0) * (*poro0)("cell",c)
                                        * (*cell_volume0)("cell",c);

    // add the time derivative of total water content to the residual
    (*g)("cell",c) += (wc1 - wc0)/dt;
  }
};


// -----------------------------------------------------------------------------
// Update variables, like densities, saturations, etc, from constitutive models.
// -----------------------------------------------------------------------------
void Permafrost::UpdateSecondaryVariables_(const Teuchos::RCP<State>& S) {
  // calculate liquid properties
  UpdateDensityLiquid_(S);
  UpdateViscosityLiquid_(S);

  // calculate ice properties
  UpdateDensityIce_(S);

  // calculate molar fraction of vapor and density of gas
  UpdateDensityGas_(S);

  // calculate saturations using WRM
  UpdateSaturation_(S);

  // update abs perm if needed
  if (variable_abs_perm_) {
    // UpdateAbsolutePermeability_(S);
  }

  // calculate rel perm using WRM
  UpdateRelativePermeability_(S);
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
// -------------------------------------------------------------
// Update the density of ice in state S.
// -------------------------------------------------------------
void Permafrost::UpdateDensityIce_(const Teuchos::RCP<State>& S) {
  Teuchos::RCP<const CompositeVector> temp = S->GetFieldData("temperature");
  Teuchos::RCP<const CompositeVector> pres = S->GetFieldData("pressure");
  Teuchos::RCP<CompositeVector> rho_ice = S->GetFieldData("density_ice", "flow");
  Teuchos::RCP<CompositeVector> n_ice =
    S->GetFieldData("molar_density_ice", "flow");

  DensityIce_(S, *temp, *pres, rho_ice, n_ice);
};


// -------------------------------------------------------------
// Evaluate EOS of the ice phase.
// -------------------------------------------------------------
void Permafrost::DensityIce_(const Teuchos::RCP<State>& S,
        const CompositeVector& temp, const CompositeVector& pres,
        const Teuchos::RCP<CompositeVector>& rho_ice,
        const Teuchos::RCP<CompositeVector>& n_ice) {
  double Mw = eos_ice_->molar_mass();

  int c_owned = rho_ice->size("cell");
  for (int c=0; c!=c_owned; ++c) {
    double rho = eos_ice_->MassDensity(temp("cell",c), pres("cell",c));
    (*rho_ice)("cell",c) = rho;
    (*n_ice)("cell",c) = rho/Mw;
  }
};


// -------------------------------------------------------------
// Update saturation of all phases in state S.
// -------------------------------------------------------------
void Permafrost::UpdateSaturation_(const Teuchos::RCP<State>& S) {
  Teuchos::RCP<const CompositeVector> temp = S->GetFieldData("temperature");
  Teuchos::RCP<const CompositeVector> pres = S->GetFieldData("pressure");
  Teuchos::RCP<const double> p_atm = S->GetScalarData("atmospheric_pressure");

  Teuchos::RCP<CompositeVector> sat_gas =
    S->GetFieldData("saturation_gas", "flow");
  Teuchos::RCP<CompositeVector> sat_liq =
    S->GetFieldData("saturation_liquid", "flow");
  Teuchos::RCP<CompositeVector> sat_ice =
    S->GetFieldData("saturation_ice", "flow");

  Teuchos::RCP<const CompositeVector> n_ice = S->GetFieldData("molar_density_ice");
  Teuchos::RCP<const CompositeVector> rho_ice = S->GetFieldData("density_ice");

  // store frozen saturation S*(sig_gl/sig_il * pc_il) in sat_ice, this is 1/A
  FrozenSaturation_(S, *temp, *rho_ice, *n_ice, sat_ice);

  // store unfrozen saturation S*(pc_gl) in sat_gas, this is 1/B
  UnfrozenSaturation_(S, *pres, *p_atm, sat_gas);

  int c_owned = sat_liq->size("cell");
  for (int c=0; c!=c_owned; ++c) {
    // see notes, section 1.2.2
    double A = 1.0/(*sat_ice)("cell",c);
    double B = 1.0/(*sat_gas)("cell",c);

    double s_liq = 1.0/( A + B - 1.0);
    double s_gas = s_liq * (B - 1.0);

    (*sat_ice)("cell",c) = 1.0 - s_liq - s_gas;
    (*sat_liq)("cell",c) = s_liq;
    (*sat_gas)("cell",c) = s_gas;
  }
};


// -------------------------------------------------------------
// Call WRM to determine the unfrozen saturation, sat^*
// -------------------------------------------------------------
void Permafrost::UnfrozenSaturation_(const Teuchos::RCP<State>& S,
        const CompositeVector& pres, const double& p_atm,
        const Teuchos::RCP<CompositeVector>& sat_star) {
  // loop over region/wrm pairs
  for (std::vector< Teuchos::RCP<WRMRegionPair> >::iterator wrm=wrm_.begin();
       wrm!=wrm_.end(); ++wrm) {
    // get the owned cells in that region
    std::string region = (*wrm)->first;
    int ncells = sat_star->size("cell");
    std::vector<int> cells(ncells);
    S->Mesh()->get_set_entities(region, AmanziMesh::CELL, AmanziMesh::OWNED, &cells);

    // use the wrm to evaluate saturation on each cell in the region
    for (std::vector<int>::iterator c=cells.begin(); c!=cells.end(); ++c) {
      (*sat_star)("cell",*c) = (*wrm)->second->saturation(p_atm - pres("cell",*c));
    }
  }
};


// -------------------------------------------------------------
// Call WRM to calculate the unfrozen dsat^*/dp
// -------------------------------------------------------------
void Permafrost::DUnfrozenSaturationDp_(const Teuchos::RCP<State>& S,
        const CompositeVector& pres, const double& p_atm,
        const Teuchos::RCP<CompositeVector>& dsat_star) {
  // loop over region/wrm pairs
  for (std::vector< Teuchos::RCP<WRMRegionPair> >::iterator wrm=wrm_.begin();
       wrm!=wrm_.end(); ++wrm) {
    // get the owned cells in that region
    std::string region = (*wrm)->first;
    int ncells = dsat_star->size("cell");
    std::vector<int> cells(ncells);
    S->Mesh()->get_set_entities(region, AmanziMesh::CELL, AmanziMesh::OWNED, &cells);

    // use the wrm to evaluate saturation on each cell in the region
    for (std::vector<int>::iterator c=cells.begin(); c!=cells.end(); ++c) {
      (*dsat_star)("cell",*c) = -(*wrm)->second->d_saturation(p_atm - pres("cell",*c));
    }
  }
};


// -------------------------------------------------------------
// Call WRM to determine the frozen saturation, sat^*
// -------------------------------------------------------------
void Permafrost::FrozenSaturation_(const Teuchos::RCP<State>& S,
        const CompositeVector& temp,
        const CompositeVector& rho_ice,
        const CompositeVector& n_ice,
        const Teuchos::RCP<CompositeVector>& sat_star) {
  // This would have to be fixed for multiple pc models
  const CompositeVector* dens_ice;
  if (pc_ice_liq_model_->IsMolarBasis()) {
    dens_ice = &n_ice;
  } else {
    dens_ice = &rho_ice;
  }

  // loop over region/wrm pairs
  for (std::vector< Teuchos::RCP<WRMRegionPair> >::iterator wrm=wrm_.begin();
       wrm!=wrm_.end(); ++wrm) {
    // get the owned cells in that region
    std::string region = (*wrm)->first;
    int ncells = sat_star->size("cell");
    std::vector<int> cells(ncells);
    S->Mesh()->get_set_entities(region, AmanziMesh::CELL, AmanziMesh::OWNED, &cells);

    for (std::vector<int>::iterator c=cells.begin(); c!=cells.end(); ++c) {
      // use a pc model to evaluate the ice pressure and then rescale it to
      // the gas-liquid retention curve.  See notes eqns 6 and 9.
      double pc_rescaled = pc_ice_liq_model_->CapillaryPressure(temp("cell",*c), (*dens_ice)("cell",*c));

      // use the wrm to evaluate saturation on each cell in the region
      (*sat_star)("cell",*c) = (*wrm)->second->saturation(pc_rescaled);
    }
  }
};


// -------------------------------------------------------------
// Evaluate WRM for Krel
// -------------------------------------------------------------
void Permafrost::UpdateRelativePermeability_(const Teuchos::RCP<State>& S) {
  Teuchos::RCP<const CompositeVector> sat_liq = S->GetFieldData("saturation_liquid");
  Teuchos::RCP<CompositeVector> rel_perm = S->GetFieldData("relative_permeability", "flow");
  RelativePermeability_(S, *sat_liq, rel_perm);
};


// -------------------------------------------------------------
// Evaluate WRM for Krel
// -------------------------------------------------------------
void Permafrost::RelativePermeability_(const Teuchos::RCP<State>& S,
        const CompositeVector& sat_liq,
        const Teuchos::RCP<CompositeVector>& rel_perm) {

  AmanziMesh::Entity_ID_List faces;

  // loop over region/wrm pairs
  for (std::vector< Teuchos::RCP<WRMRegionPair> >::iterator wrm=wrm_.begin();
       wrm!=wrm_.end(); ++wrm) {
    // get the owned cells in that region
    std::string region = (*wrm)->first;
    int ncells = rel_perm->size("cell");
    std::vector<int> cells(ncells);
    S->Mesh()->get_set_entities(region, AmanziMesh::CELL, AmanziMesh::OWNED, &cells);

    // use the wrm to evaluate saturation on each cell in the region
    for (std::vector<int>::iterator c=cells.begin(); c!=cells.end(); ++c) {
      double pc = (*wrm)->second->capillaryPressure(sat_liq("cell",*c));
      (*rel_perm)("cell",*c) = (*wrm)->second->k_relative(pc);

      // and also on the BC faces
      // NOTE: this is 1st order on the boundary now, but we don't
      //   have sat on teh boundary, which would require a lot of
      //   extra work.
      rel_perm->mesh()->cell_get_faces(*c, &faces);
      for (int n=0; n!=faces.size(); ++n) {
        int f = faces[n];
        if (bc_markers_[f] != Operators::MFD_BC_NULL) {
          (*rel_perm)("face",f) = (*rel_perm)("cell", *c);
        }
      }
    }
  }
};


} //namespace
} //namespace
