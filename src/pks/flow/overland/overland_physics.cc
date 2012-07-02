/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  A base two-phase, thermal Richard's equation with water vapor.

  License: BSD
  Authors: Neil Carlson (version 1)
  Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
  Ethan Coon (ATS version) (ecoon@lanl.gov)
*/

#include "overland.hh"

namespace Amanzi {
namespace Flow {

#if 0
}}
#endif

// -------------------------------------------------------------
// Diffusion term, div K grad T
// -------------------------------------------------------------
void OverlandFlow::ApplyDiffusion_(const Teuchos::RCP<State>& S,
                                   const Teuchos::RCP<CompositeVector>& g) {

  // update the rel perm according to the scheme of choice
  UpdatePermeabilityData_(S);
  Teuchos::RCP<const CompositeVector> upwind_conductivity =
    S->GetFieldData("upwind_overland_conductivity", "overland_flow");

  // update the stiffness matrix
  matrix_->CreateMFDstiffnessMatrices(*upwind_conductivity);
  matrix_->CreateMFDrhsVectors();

  //  std::cout << "BC in res: " << bc_values_[3] << ", " << bc_values_[5] << std::endl;
  matrix_->ApplyBoundaryConditions(bc_markers_, bc_values_);
  matrix_->AssembleGlobalMatrices();

  // calculate the residual
  Teuchos::RCP<const CompositeVector> pres_elev = S->GetFieldData("pres_elev");
  matrix_->ComputeNegativeResidual(*pres_elev, g);
};


// -------------------------------------------------------------
// Accumulation of internal energy term du/dt
// -------------------------------------------------------------
void OverlandFlow::AddAccumulation_(const Teuchos::RCP<CompositeVector>& g) {
  const CompositeVector & pres0        = *(S_inter_->GetFieldData("overland_pressure"));
  const CompositeVector & pres1        = *(S_next_ ->GetFieldData("overland_pressure"));
  const CompositeVector & cell_volume0 = *(S_inter_->GetFieldData("surface_cell_volume"));
  const CompositeVector & cell_volume1 = *(S_next_ ->GetFieldData("surface_cell_volume"));

  double dt = S_next_->time() - S_inter_->time();

  int c_owned = g->size("cell");
  for (int c=0; c!=c_owned; ++c) {
    (*g)("cell",0,c) += (cell_volume1("cell",c)*pres1("cell",c)
                         - cell_volume0("cell",c)*pres0("cell",c))/dt ;
  }
};


// -------------------------------------------------------------
// Source term
// -------------------------------------------------------------
void OverlandFlow::AddLoadValue_(const Teuchos::RCP<State>& S,
        const Teuchos::RCP<CompositeVector>& g) {

  Teuchos::RCP <const CompositeVector> cell_volume = S->GetFieldData("surface_cell_volume");
  Teuchos::RCP <const CompositeVector> rain = S->GetFieldData("rainfall_rate");

  int c_owned = g->size("cell");
  for (int c=0; c!=c_owned; ++c) {
    (*g)("cell",c) -= (*rain)("cell",c) * (*cell_volume)("cell",c);
  }
};


// -----------------------------------------------------------------------------
// Update variables, like p + z and rel perm.
// -----------------------------------------------------------------------------
void OverlandFlow::UpdateSecondaryVariables_(const Teuchos::RCP<State>& S) {
  // add elevation
  AddElevation_(S) ;

  // compute non-linear coefficient
  UpdateConductivity_(S);

  // update rainfall
  UpdateLoadValue_(S);
};


// -----------------------------------------------------------------------------
// Evaluate model for elevation + pressure
// -----------------------------------------------------------------------------
void OverlandFlow::AddElevation_(const Teuchos::RCP<State>& S) {

  const CompositeVector & pres      = *(S->GetFieldData("overland_pressure"));
  const CompositeVector & elevation = *(S->GetFieldData("elevation"));
  CompositeVector       & pres_elev = *(S->GetFieldData("pres_elev","overland_flow"));

  // add elevation to pres_elev, cell dofs
  int c_owned = pres_elev.size("cell");
  for (int c=0; c!=c_owned; ++c) {
    pres_elev("cell",c) = pres("cell",c) + elevation("cell",c) ;
  }

  // add elevation to pres_elev, face dofs
  int f_owned = pres_elev.size("face");
  for (int f=0; f!=f_owned; ++f) {
    pres_elev("face",f) = pres("face",f) + elevation("face",f) ;
  }
}


// -----------------------------------------------------------------------------
// Update the rainfall rate.
// -----------------------------------------------------------------------------
void OverlandFlow::UpdateLoadValue_(const Teuchos::RCP<State>& S) {
  Teuchos::RCP<CompositeVector> rain_rate =
    S->GetFieldData("rainfall_rate", "overland_flow");
  rain_rate_function_->Compute(S->time(), rain_rate.ptr());
};


// -----------------------------------------------------------------------------
// Update the conductivity.
// -----------------------------------------------------------------------------
void OverlandFlow::UpdateConductivity_(const Teuchos::RCP<State>& S) {
  Teuchos::RCP<const CompositeVector> pres = S->GetFieldData("overland_pressure");
  Teuchos::RCP<const CompositeVector> mann = S->GetFieldData("manning_coef");
  Teuchos::RCP<const CompositeVector> slope = S->GetFieldData("slope_magnitude");

  Teuchos::RCP<CompositeVector> cond =
    S->GetFieldData("overland_conductivity", "overland_flow");

  Conductivity_(S, *pres, *mann, *slope, cond);
};


// -----------------------------------------------------------------------------
// Update the "relative permeability".
// -----------------------------------------------------------------------------
void OverlandFlow::Conductivity_(const Teuchos::RCP<State>& S,
                                 const CompositeVector& pressure,
                                 const CompositeVector& manning_coef,
                                 const CompositeVector& slope_mag,
                                 const Teuchos::RCP<CompositeVector>& conductivity) {
  double eps = 1.e-14;
  double exponent = 1.0 + manning_exp_;

  int ncells = conductivity->size("cell");
  for (int c=0; c!=ncells; ++c ) {
    double scaling = manning_coef("cell",c) * std::sqrt(slope_mag("cell",c) + eps);
    (*conductivity)("cell",c) = std::pow(pressure("cell",c), exponent) / scaling;
  }
};

} //namespace
} //namespace
