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

void OverlandFlow::ApplyDiffusion_(const Teuchos::RCP<State>& S,
                                   const Teuchos::RCP<CompositeVector>& g) {

  // update the rel perm according to the scheme of choice
  UpdatePermeabilityData_(S);
  Teuchos::RCP<const CompositeVector> upwind_conductivity =
    S->GetFieldData("upwind_overland_conductivity", "overland_flow");

  // update the stiffness matrix
  matrix_->CreateMFDstiffnessMatrices(*upwind_conductivity);
  matrix_->CreateMFDrhsVectors();

  std::cout << "BC in res: " << bc_values_[3] << ", " << bc_values_[5] << std::endl;
  matrix_->ApplyBoundaryConditions(bc_markers_, bc_values_);
  matrix_->AssembleGlobalMatrices();

  // calculate the residual
  const CompositeVector& pres_elev = *(S->GetFieldData("pres_elev"));
  matrix_->ComputeNegativeResidual( pres_elev, g );

  // std::cout << "  res0 (after diff, pre fix): " << (*g)("cell",0,0) << " " << (*g)("face",0,3) << std::endl;
  // std::cout << "  res1 (after diff, pre fix): " << (*g)("cell",0,9) << " " << (*g)("face",0,29) << std::endl;

  // // 
  // BoundaryFunction::Iterator bc;
  // const CompositeVector & elevation = *(S->GetFieldData("elevation"));
  // for (bc=bc_pressure_->begin(); bc!=bc_pressure_->end(); ++bc) {
  //   int f = bc->first;
  //   (*g)("face",0,f) += elevation("face",0,f);
  // }
};

void OverlandFlow::AddElevation_(const Teuchos::RCP<State>& S) {

  const CompositeVector & pres      = *(S->GetFieldData("overland_pressure"));
  const CompositeVector & elevation = *(S->GetFieldData("elevation"));
  CompositeVector       & pres_elev = *(S->GetFieldData("pres_elev","overland_flow"));

  // add elevation to pres_elev, cell dofs
  int c_owned = S->mesh()->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c=0; c!=c_owned; ++c) {
    pres_elev("cell",0,c) = pres("cell",0,c) + elevation("cell",0,c) ;
  }

  // add elevation to pres_elev, face dofs
  int f_owned = S->mesh()->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  for (int f=0; f!=f_owned; ++f) {
    pres_elev("face",0,f) = pres("face",0,f) + elevation("face",0,f) ;
  }
}

void OverlandFlow::AddAccumulation_(const Teuchos::RCP<CompositeVector>& g) {
  const CompositeVector & pres0        = *(S_inter_->GetFieldData("overland_pressure"));
  const CompositeVector & pres1        = *(S_next_ ->GetFieldData("overland_pressure"));
  const CompositeVector & cell_volume0 = *(S_inter_->GetFieldData("cell_volume"));
  const CompositeVector & cell_volume1 = *(S_next_ ->GetFieldData("cell_volume"));

  double dt = S_next_->time() - S_inter_->time();

  int c_owned = S_->mesh()->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c=0; c!=c_owned; ++c) {
    (*g)("cell",0,c) += (cell_volume1("cell",0,c)*pres1("cell",0,c)
                         - cell_volume0("cell",0,c)*pres0("cell",0,c))/dt ;
  }
};

void OverlandFlow::AddLoadValue_(const Teuchos::RCP<State>& S,
                                 const Teuchos::RCP<CompositeVector>& g) {
  int c_owned = S->mesh()->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  const CompositeVector& cell_volume = *(S->GetFieldData("cell_volume"));
  for (int c=0; c!=c_owned; ++c) {
    (*g)("cell",0,c) -= rhs_load_value(S->time()) * cell_volume("cell",0,c) ;
  }
}

/* ******************************************************************
 * Update secondary variables, calculated in various methods below.
 ****************************************************************** */
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

/* ******************************************************************
 * OVERLAND Update secondary variables, calculated in various methods below.
 ****************************************************************** */
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

/* ******************************************************************
 * Converts absolute perm to tensor
 ****************************************************************** */
void OverlandFlow::SetAbsolutePermeabilityTensor_(const Teuchos::RCP<State>& S) {
  // perm is not needed here
  for (int c=0; c!=K_.size(); ++c) {
    K_[c](0, 0) = 1.0;
  }
};

} //namespace
} //namespace
