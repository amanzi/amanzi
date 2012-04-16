/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  A base two-phase, thermal Richard's equation with water vapor.

  License: BSD
  Authors: Neil Carlson (version 1)
  Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
  Ethan Coon (ATS version) (ecoon@lanl.gov)
*/

#include "permafrost.hh"

namespace Amanzi {
namespace Flow {

void Permafrost::ApplyDiffusion_(const Teuchos::RCP<State>& S,
        const Teuchos::RCP<CompositeVector>& g) {

  // update the rel perm according to the scheme of choice
  UpdatePermeabilityData_(S);
  Teuchos::RCP<const CompositeVector> rel_perm_faces =
    S->GetFieldData("rel_perm_faces", "flow");

  // update the stiffness matrix
  matrix_->CreateMFDstiffnessMatrices(K_, rel_perm_faces);
  matrix_->CreateMFDrhsVectors();
  AddGravityFluxes_(S, matrix_);
  matrix_->ApplyBoundaryConditions(bc_markers_, bc_values_);
  matrix_->AssembleGlobalMatrices();

  // calculate the residual
  Teuchos::RCP<const CompositeVector> pres = S->GetFieldData("pressure");
  matrix_->ComputeNegativeResidual(*pres, g);
};


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


  int c_owned = S_->mesh()->count_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c=0; c!=c_owned; ++c) {
    // calculate water content of each phase and at each time
    double wc_liq1 = (*n_liq1)(c) * (*sat_liq1)(c);
    double wc_ice1 = (*n_ice1)(c) * (*sat_ice1)(c);
    double wc_gas1 = (*n_gas1)(c) * (*sat_gas1)(c) * (*mol_frac_gas1)(c);
    double wc1 = (wc_liq1 + wc_gas1 + wc_ice1) * (*poro1)(c) * (*cell_volume1)(c);

    double wc_liq0 = (*n_liq0)(c) * (*sat_liq0)(c);
    double wc_ice0 = (*n_ice0)(c) * (*sat_ice0)(c);
    double wc_gas0 = (*n_gas0)(c) * (*sat_gas0)(c) * (*mol_frac_gas0)(c);
    double wc0 = (wc_liq0 + wc_gas0 + wc_ice0) * (*poro0)(c) * (*cell_volume0)(c);

    // add the time derivative of total water content to the residual
    (*g)("cell",0,c) += (wc1 - wc0)/dt;
  }
};

/* ******************************************************************
 * Update secondary variables, calculated in various methods below.
 ****************************************************************** */
void Permafrost::UpdateSecondaryVariables_(const Teuchos::RCP<State>& S) {
  // get needed fields
  Teuchos::RCP<const CompositeVector> temp = S->GetFieldData("temperature");
  Teuchos::RCP<const CompositeVector> pres = S->GetFieldData("pressure");
  Teuchos::RCP<const double> p_atm = S->GetScalarData("atmospheric_pressure");

  Teuchos::RCP<CompositeVector> dens_liq = S->GetFieldData("density_liquid", "flow");
  Teuchos::RCP<CompositeVector> mol_dens_liq = S->GetFieldData("molar_density_liquid", "flow");
  Teuchos::RCP<CompositeVector> visc_liq = S->GetFieldData("viscosity_liquid", "flow");

  Teuchos::RCP<CompositeVector> dens_ice = S->GetFieldData("density_ice", "flow");
  Teuchos::RCP<CompositeVector> mol_dens_ice = S->GetFieldData("molar_density_ice", "flow");
  Teuchos::RCP<CompositeVector> visc_ice = S->GetFieldData("viscosity_ice", "flow");

  Teuchos::RCP<CompositeVector> dens_gas = S->GetFieldData("density_gas", "flow");
  Teuchos::RCP<CompositeVector> mol_dens_gas = S->GetFieldData("molar_density_gas", "flow");
  Teuchos::RCP<CompositeVector> mol_frac_gas = S->GetFieldData("mol_frac_gas", "flow");

  Teuchos::RCP<CompositeVector> sat_gas = S->GetFieldData("saturation_gas", "flow");
  Teuchos::RCP<CompositeVector> sat_liq = S->GetFieldData("saturation_liquid", "flow");
  Teuchos::RCP<CompositeVector> sat_ice = S->GetFieldData("saturation_ice", "flow");
  Teuchos::RCP<CompositeVector> rel_perm = S->GetFieldData("relative_permeability", "flow");

  Teuchos::RCP<CompositeVector> flux = S->GetFieldData("darcy_flux", "flow");

  // calculate liquid properties
  DensityLiquid_(S, *temp, *pres, dens_liq, mol_dens_liq);
  ViscosityLiquid_(S, *temp, visc_liq);

  // calculate liquid properties
  DensityIce_(S, *temp, *pres, dens_ice, mol_dens_ice);

  // calculate molar fraction of vapor and density of gas
  DensityGas_(S, *temp, *pres, *p_atm, mol_frac_gas, dens_gas, mol_dens_gas);

  // calculate saturations using WRM
  UpdateSaturation_(S);

  // update abs perm if needed
  if (variable_abs_perm_) {
    // Teuchos::RCP<const CompositeVector> phi = S->GetFieldData("porosity");
    // Teuchos::RCP<CompositeVector> abs_perm = S->GetFieldData("permeability", "flow");
    // AbsolutePermeability_(*phi, abs_perm)
  }

  // calculate rel perm using WRM
  RelativePermeability_(S, *pres, *p_atm, rel_perm);
};

void Permafrost::DensityLiquid_(const Teuchos::RCP<State>& S,
        const CompositeVector& temp, const CompositeVector& pres,
        const Teuchos::RCP<CompositeVector>& dens_liq,
        const Teuchos::RCP<CompositeVector>& mol_dens_liq) {

  int c_owned = S->mesh()->count_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  double Mw = eos_liquid_->molar_mass();
  for (int c=0; c!=c_owned; ++c) {
    double rho = eos_liquid_->MassDensity(temp("cell",0,c), pres("cell",0,c));
    (*dens_liq)("cell",0,c) = rho;
    (*mol_dens_liq)("cell",0,c) = rho/Mw;
  }
};

void Permafrost::ViscosityLiquid_(const Teuchos::RCP<State>& S,
        const CompositeVector& temp,
        const Teuchos::RCP<CompositeVector>& visc_liq) {
  int c_owned = S->mesh()->count_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c=0; c!=c_owned; ++c) {
    (*visc_liq)("cell",0,c) = eos_liquid_->Viscosity(temp("cell",0,c));
  }
};

void Permafrost::DensityGas_(const Teuchos::RCP<State>& S,
                           const CompositeVector& temp,
                           const CompositeVector& pres, const double& p_atm,
                           const Teuchos::RCP<CompositeVector>& mol_frac_gas,
                           const Teuchos::RCP<CompositeVector>& dens_gas,
                           const Teuchos::RCP<CompositeVector>& mol_dens_gas) {

  int c_owned = S->mesh()->count_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c=0; c!=c_owned; ++c) {
    double p_sat = eos_gas_->SaturatedVaporPressure(temp("cell",0,c));
    (*mol_frac_gas)("cell",0,c) = p_sat/p_atm;
    double Mv = eos_gas_->molar_mass((*mol_frac_gas)("cell",0,c));
    double n = eos_gas_->MolarDensity(temp("cell",0,c), pres("cell",0,c));
    (*dens_gas)("cell",0,c) = Mv*n;
    (*mol_dens_gas)("cell",0,c) = n;
  }
};

void Permafrost::UpdateSaturation_(const Teuchos::RCP<State>& S) {
  Teuchos::RCP<const CompositeVector> temp = S->GetFieldData("temperature");
  Teuchos::RCP<const CompositeVector> pres = S->GetFieldData("pressure");
  Teuchos::RCP<const double> p_atm = S->GetScalarData("atmospheric_pressure");

  Teuchos::RCP<CompositeVector> sat_gas =
    S->GetFieldData("saturation_gas", "falow");
  Teuchos::RCP<CompositeVector> sat_liq =
    S->GetFieldData("saturation_liquid", "flow");
  Teuchos::RCP<CompositeVector> sat_ice =
    S->GetFieldData("saturation_ice", "flow");

  Teuchos::RCP<const CompositeVector> n_ice = S->GetFieldData("molar_density_ice");
  Teuchos::RCP<const CompositeVector> rho_ice = S->GetFieldData("density_ice");

  // store frozen saturation S*_il(pc_il) in sat_ice, this is 1/A
  FrozenSaturation_(S, *temp, *rho_ice, *n_ice, sat_ice);

  // store unfrozen saturation S*_gl(pc_gl) in sat_liq, this is 1/B
  UnfrozenSaturation_(S, *pres, *p_atm, sat_liq);

  int c_owned = S->mesh()->count_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c=0; c!=c_owned; ++c) {
    // see notes, eqn 10?
    double A = 1.0/(*sat_ice)("cell",0,c);
    double B = 1.0/(*sat_liq)("cell",0,c);

    double s_liq = 1.0/( A + B - 1.0);
    double s_gas = s_liq * (B - 1.0);

    (*sat_ice)("cell",0,c) = 1.0 - s_liq - s_gas;
    (*sat_liq)("cell",0,c) = s_liq;
    (*sat_gas)("cell",0,c) = s_gas;
  }
};

void Permafrost::UnfrozenSaturation_(const Teuchos::RCP<State>& S,
        const CompositeVector& pres, const double& p_atm,
        const Teuchos::RCP<CompositeVector>& sat_star) {
  // loop over region/wrm pairs
  for (std::vector< Teuchos::RCP<WRMRegionPair> >::iterator wrm=wrm_liq_gas_.begin();
       wrm!=wrm_liq_gas_.end(); ++wrm) {
    // get the owned cells in that region
    std::string region = (*wrm)->first;
    int ncells = S->mesh()->get_set_size(region, AmanziMesh::CELL, AmanziMesh::OWNED);
    std::vector<unsigned int> cells(ncells);
    S->mesh()->get_set_entities(region, AmanziMesh::CELL, AmanziMesh::OWNED, &cells);

    // use the wrm to evaluate saturation on each cell in the region
    for (std::vector<unsigned int>::iterator c=cells.begin(); c!=cells.end(); ++c) {
      (*sat_star)("cell",0,*c) = (*wrm)->second->saturation(p_atm - pres("cell",0,*c));
    }
  }
};

void Permafrost::DUnfrozenSaturationDp_(const Teuchos::RCP<State>& S,
        const CompositeVector& pres, const double& p_atm,
        const Teuchos::RCP<CompositeVector>& dsat_star) {
  // loop over region/wrm pairs
  for (std::vector< Teuchos::RCP<WRMRegionPair> >::iterator wrm=wrm_liq_gas_.begin();
       wrm!=wrm_liq_gas_.end(); ++wrm) {
    // get the owned cells in that region
    std::string region = (*wrm)->first;
    int ncells = S->mesh()->get_set_size(region, AmanziMesh::CELL, AmanziMesh::OWNED);
    std::vector<unsigned int> cells(ncells);
    S->mesh()->get_set_entities(region, AmanziMesh::CELL, AmanziMesh::OWNED, &cells);

    // use the wrm to evaluate saturation on each cell in the region
    for (std::vector<unsigned int>::iterator c=cells.begin(); c!=cells.end(); ++c) {
      (*dsat_star)("cell",0,*c) = -(*wrm)->second->d_saturation(p_atm - pres("cell",0,*c));
    }
  }
};

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
  for (std::vector< Teuchos::RCP<WRMRegionPair> >::iterator wrm=wrm_ice_liq_.begin();
       wrm!=wrm_ice_liq_.end(); ++wrm) {
    // get the owned cells in that region
    std::string region = (*wrm)->first;
    int ncells = S->mesh()->get_set_size(region, AmanziMesh::CELL, AmanziMesh::OWNED);
    std::vector<unsigned int> cells(ncells);
    S->mesh()->get_set_entities(region, AmanziMesh::CELL, AmanziMesh::OWNED, &cells);

    for (std::vector<unsigned int>::iterator c=cells.begin(); c!=cells.end(); ++c) {
      // use a pc model to evaluate the ice pressure
      double pc_ice_liq = pc_ice_liq_model_->CapillaryPressure(temp("cell",0,*c), (*dens_ice)("cell",0,*c));

      // use the wrm to evaluate saturation on each cell in the region
      (*sat_star)("cell",0,*c) = (*wrm)->second->saturation(pc_ice_liq);
    }
  }
};

void Permafrost::RelativePermeability_(const Teuchos::RCP<State>& S,
        const CompositeVector& pres, const double& p_atm,
        const Teuchos::RCP<CompositeVector>& rel_perm) {

  // loop over region/wrm pairs
  for (std::vector< Teuchos::RCP<WRMRegionPair> >::iterator wrm=wrm_.begin();
       wrm!=wrm_.end(); ++wrm) {
    // get the owned cells in that region
    std::string region = (*wrm)->first;
    int ncells = S->mesh()->get_set_size(region, AmanziMesh::CELL, AmanziMesh::OWNED);
    std::vector<unsigned int> cells(ncells);
    S->mesh()->get_set_entities(region, AmanziMesh::CELL, AmanziMesh::OWNED, &cells);

    // use the wrm to evaluate saturation on each cell in the region
    for (std::vector<unsigned int>::iterator c=cells.begin(); c!=cells.end(); ++c) {
      (*rel_perm)("cell",0,*c) = (*wrm)->second->k_relative(p_atm - pres("cell",0,*c));
    }
  }
};


/* ******************************************************************
 * Converts absolute perm to tensor
 ****************************************************************** */
void Permafrost::SetAbsolutePermeabilityTensor_(const Teuchos::RCP<State>& S) {
  // currently assumes isotropic perm, should be updated
  Teuchos::RCP<const CompositeVector> perm = S->GetFieldData("permeability");
  for (int c=0; c!=K_.size(); ++c) {
    K_[c](0, 0) = (*perm)(c);
  }
};

/* ******************************************************************
 * Routine updates elemental discretization matrices and must be
 * called before applying boundary conditions and global assembling.
 ****************************************************************** */
void Permafrost::AddGravityFluxes_(const Teuchos::RCP<State>& S,
        const Teuchos::RCP<Operators::MatrixMFD>& matrix) {

  Teuchos::RCP<const CompositeVector> rho = S->GetFieldData("density_liquid");
  Teuchos::RCP<const Epetra_Vector> g_vec = S->GetConstantVectorData("gravity");
  Teuchos::RCP<const CompositeVector> Krel = S->GetFieldData("rel_perm_faces");

  AmanziGeometry::Point gravity(g_vec->MyLength());
  for (int i=0; i!=g_vec->MyLength(); ++i) gravity[i] = (*g_vec)[i];

  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  int c_owned = S->mesh()->count_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c=0; c!=c_owned; ++c) {
    S->mesh()->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    Epetra_SerialDenseVector& Ff = matrix->Ff_cells()[c];
    double& Fc = matrix->Fc_cells()[c];

    for (int n=0; n!=nfaces; ++n) {
      int f = faces[n];
      const AmanziGeometry::Point& normal = S->mesh()->face_normal(f);

      double outward_flux = ((K_[c] * gravity) * normal)
        * dirs[n] * (*Krel)("face",0,f) * (*rho)("cell",0,c);
      Ff[n] += outward_flux;
      Fc -= outward_flux;  // Nonzero-sum contribution when not upwinding
    }
  }
};


/* ******************************************************************
 * Updates global Darcy vector calculated by a discretization method.
 ****************************************************************** */
void Permafrost::AddGravityFluxesToVector_(const Teuchos::RCP<State>& S,
        const Teuchos::RCP<CompositeVector>& darcy_flux) {

  Teuchos::RCP<const CompositeVector> rho = S->GetFieldData("density_liquid");
  Teuchos::RCP<const Epetra_Vector> g_vec = S->GetConstantVectorData("gravity");
  Teuchos::RCP<const CompositeVector> Krel = S->GetFieldData("rel_perm_faces");

  AmanziGeometry::Point gravity(g_vec->MyLength());
  for (int i=0; i!=g_vec->MyLength(); ++i) gravity[i] = (*g_vec)[i];

  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  int f_used = S->mesh()->count_entities(AmanziMesh::FACE, AmanziMesh::USED);
  int f_owned = S->mesh()->count_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  std::vector<bool> done(f_used, false);

  int c_owned = S->mesh()->count_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c=0; c!=c_owned; ++c) {
    S->mesh()->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    for (int n=0; n!=nfaces; ++n) {
      int f = faces[n];
      const AmanziGeometry::Point& normal = S->mesh()->face_normal(f);

      if (f<f_owned && !done[f]) {
        (*darcy_flux)(f) += ((K_[c] * gravity) * normal)
          * (*Krel)("face",0,f) * (*rho)("cell",0,c);
        done[f] = true;
      }
    }
  }
};

} //namespace
} //namespace
