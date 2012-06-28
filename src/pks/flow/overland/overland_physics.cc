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
  const CompositeVector & cell_volume0 = *(S_inter_->GetFieldData("cell_volume"));
  const CompositeVector & cell_volume1 = *(S_next_ ->GetFieldData("cell_volume"));

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
  int c_owned = g->size("cell");
  const CompositeVector& cell_volume = *(S->GetFieldData("cell_volume"));
  for (int c=0; c!=c_owned; ++c) {
    (*g)("cell",c) -= rhs_load_value() * cell_volume("cell",c) ;
  }
}


// -----------------------------------------------------------------------------
// Update variables, like p + z and rel perm.
// -----------------------------------------------------------------------------
void OverlandFlow::UpdateSecondaryVariables_(const Teuchos::RCP<State>& S) {
  // get needed fields
  const CompositeVector & pres = *(S->GetFieldData ("overland_pressure"));
  Teuchos::RCP<CompositeVector> cond =
    S->GetFieldData("overland_conductivity", "overland_flow");

  // add elevation
  AddElevation_(S) ;

  // compute non-linear coefficient
  RelativePermeability_(S, pres, cond);
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


/* ******************************************************************
 * OVERLAND Update secondary variables, calculated in various methods below.
 ****************************************************************** */
#if 0
void OverlandFlow::RelativePermeability_( const Teuchos::RCP<State>& S,
                                          const CompositeVector & pres, 
                                          const Teuchos::RCP<CompositeVector>& rel_perm ) {

  double manning_coeff = manning[0] ;
  double slope_coeff   = slope_x[0]+1.e-14 ;
  double scaling       = manning_coeff * std::sqrt(slope_coeff) ;

  int ncells = S->mesh()->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c=0; c!=ncells; ++c ) {
    (*rel_perm)("cell",0,c) = pow( pres("cell",0,c), 5./3. )/scaling ;
  }
};
#else
void OverlandFlow::RelativePermeability_( const Teuchos::RCP<State>& S,
                                          const CompositeVector & pres, 
                                          const Teuchos::RCP<CompositeVector>& rel_perm ) {

  int ncells = rel_perm->size("cell");
  for (int c=0; c!=ncells; ++c ) {
    // get cell center coords
    Amanzi::AmanziGeometry::Point pc = S->Mesh()->cell_centroid(c);
    // get coefficients
    int izn = TestTwoZoneFlag_( pc[0], pc[1] ) ;
    //int izn = 0 ; // 1D-pblm, single zone
    double manning_coeff = manning[izn] ;
    double slope_coeff = std::sqrt(std::pow(slope_x[izn],2)+std::pow(slope_y[izn],2))+1.e-14 ;
    double scaling = manning_coeff * std::sqrt(slope_coeff) ;
    // compute the krel term
    (*rel_perm)("cell",c) = std::pow( pres("cell",c), 5./3. )/scaling ;
  }
};
#endif

} //namespace
} //namespace
