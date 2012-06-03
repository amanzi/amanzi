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
  Teuchos::RCP<const CompositeVector> rel_perm_faces =
    S->GetFieldData("rel_perm_faces", "overland_flow");
  
  // update the stiffness matrix
  matrix_->CreateMFDstiffnessMatrices( K_, rel_perm_faces );
  matrix_->CreateMFDrhsVectors();
  
  // remove
  //AddGravityFluxes_(S, matrix_);
  
  std::cout << "BC in res: " << bc_values_[3] << ", " << bc_values_[5] << std::endl;
  matrix_->ApplyBoundaryConditions(bc_markers_, bc_values_);
  matrix_->AssembleGlobalMatrices();
  
  // calculate the residual
  const CompositeVector & pres_elev = *(S->GetFieldData("pres_elev"));
  matrix_->ComputeNegativeResidual( pres_elev, g );

  // 
  BoundaryFunction::Iterator bc;
  const CompositeVector & elevation = *(S->GetFieldData("elevation"));
  for (bc=bc_pressure_->begin(); bc!=bc_pressure_->end(); ++bc) {
    int f = bc->first;
    (*g)("face",0,f) += elevation("face",0,f);
  }
};

void OverlandFlow::AddElevation_(const Teuchos::RCP<State>& S) {

  const CompositeVector & pres      = *(S->GetFieldData("overland_pressure"));
  const CompositeVector & elevation = *(S->GetFieldData("elevation"));
  CompositeVector       & pres_elev = *(S->GetFieldData("pres_elev","overland_flow"));
  
  // add elevation to pres_elev, cell dofs
  int c_owned = S->mesh()->count_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c=0; c!=c_owned; ++c) {
    pres_elev("cell",0,c) = pres("cell",0,c) + elevation("cell",0,c) ;
  }  

  // add elevation to pres_elev, face dofs
  int f_owned = S->mesh()->count_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
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
  
  int c_owned = S_->mesh()->count_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c=0; c!=c_owned; ++c) {
     (*g)("cell",0,c) += (pres1("cell",0,c)-pres0("cell",0,c))/dt ;
  }
};

void OverlandFlow::AddLoadValue_(const Teuchos::RCP<State>& S,
                                 const Teuchos::RCP<CompositeVector>& g) {
  int c_owned = S->mesh()->count_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  const CompositeVector & cell_volume = *(S->GetFieldData("cell_volume"));
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
  Teuchos::RCP<CompositeVector> rel_perm = S->GetFieldData("relative_permeability", "overland_flow");
  
  // add elevation
  AddElevation_(S) ;

  // compute non-linear coefficient
  RelativePermeability_( S, pres, rel_perm );
};


// void OverlandFlow::Saturation_( const Teuchos::RCP<State> & S,
//                                 const CompositeVector     & pres, 
//                                 const double                p_atm,
//                                 CompositeVector           & sat_liq ) {
//   // loop over region/wrm pairs
//   for (std::vector< Teuchos::RCP<WRMRegionPair> >::iterator wrm=wrm_.begin();
//        wrm!=wrm_.end(); ++wrm) {
//     // get the owned cells in that region
//     std::string region = (*wrm)->first;
//     int ncells = S->mesh()->get_set_size(region, AmanziMesh::CELL, AmanziMesh::OWNED);
//     std::vector<unsigned int> cells(ncells);
//     S->mesh()->get_set_entities(region, AmanziMesh::CELL, AmanziMesh::OWNED, &cells);

//     // use the wrm to evaluate saturation on each cell in the region
//     for (std::vector<unsigned int>::iterator c=cells.begin(); c!=cells.end(); ++c) {
//       sat_liq("cell",0,*c) = (*wrm)->second->saturation(p_atm - pres("cell",0,*c));
//     }
//   }
// };

// void OverlandFlow::DSaturationDp_(const Teuchos::RCP<State>& S,
//                                   const CompositeVector& pres, 
//                                   const double& p_atm,
//                                   const Teuchos::RCP<CompositeVector>& dsat_liq) {
//   // loop over region/wrm pairs
//   for (std::vector< Teuchos::RCP<WRMRegionPair> >::iterator wrm=wrm_.begin();
//        wrm!=wrm_.end(); ++wrm) {
//     // get the owned cells in that region
//     std::string region = (*wrm)->first;
//     int ncells = S->mesh()->get_set_size(region, AmanziMesh::CELL, AmanziMesh::OWNED);
//     std::vector<unsigned int> cells(ncells);
//     S->mesh()->get_set_entities(region, AmanziMesh::CELL, AmanziMesh::OWNED, &cells);
    
//     // use the wrm to evaluate saturation on each cell in the region
//     for (std::vector<unsigned int>::iterator c=cells.begin(); c!=cells.end(); ++c) {
//       (*dsat_liq)("cell",0,*c) = -(*wrm)->second->d_saturation(p_atm - pres("cell",0,*c));
//     }
//   }
// };

// void OverlandFlow::RelativePermeability_( const Teuchos::RCP<State>& S,
//                                           const CompositeVector    & pres, 
//                                           const double             & p_atm,
//                                           CompositeVector          & rel_perm ) {
  
//   // loop over region/wrm pairs
//   for (std::vector< Teuchos::RCP<WRMRegionPair> >::iterator wrm=wrm_.begin();
//        wrm!=wrm_.end(); ++wrm) {
//     // get the owned cells in that region
//     std::string region = (*wrm)->first;
//     int ncells = S->mesh()->get_set_size(region, AmanziMesh::CELL, AmanziMesh::OWNED);
//     std::vector<unsigned int> cells(ncells);
//     S->mesh()->get_set_entities(region, AmanziMesh::CELL, AmanziMesh::OWNED, &cells);

//     // use the wrm to evaluate saturation on each cell in the region
//     for (std::vector<unsigned int>::iterator c=cells.begin(); c!=cells.end(); ++c) {
//       rel_perm("cell",0,*c) = (*wrm)->second->k_relative(p_atm - pres("cell",0,*c));
//     }
//   }
// };

/* ******************************************************************
 * OVERLAND Update secondary variables, calculated in various methods below.
 ****************************************************************** */

void OverlandFlow::RelativePermeability_( const Teuchos::RCP<State>& S,
                                          const CompositeVector & pres, 
                                          const Teuchos::RCP<CompositeVector>& rel_perm ) {

  double manning_coeff = manning[0] ;
  double slope_coeff   = slope_x[0] ;
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
  // currently assumes isotropic perm, should be updated
  for (int c=0; c!=K_.size(); ++c) {
    K_[c](0, 0) = 1.0;
  }
};

// /* ******************************************************************
//  * Routine updates elemental discretization matrices and must be
//  * called before applying boundary conditions and global assembling.
//  ****************************************************************** */
// void OverlandFlow::AddGravityFluxes_(const Teuchos::RCP<State>& S,
//                                      const Teuchos::RCP<Operators::MatrixMFD>& matrix) {
  
//   double rho = *(S->GetScalarData("density_liquid"));
//   Teuchos::RCP<const Epetra_Vector>   g_vec = S->GetConstantVectorData("gravity");
//   Teuchos::RCP<const CompositeVector> Krel  = S->GetFieldData("rel_perm_faces");
  
//   AmanziGeometry::Point gravity(g_vec->MyLength());
//   for (int i=0; i!=g_vec->MyLength(); ++i) gravity[i] = (*g_vec)[i];
  
//   AmanziMesh::Entity_ID_List faces;
//   std::vector<int> dirs;

//   int c_owned = S->mesh()->count_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
//   for (int c=0; c!=c_owned; ++c) {
//     S->mesh()->cell_get_faces_and_dirs(c, &faces, &dirs);
//     int nfaces = faces.size();
    
//     Epetra_SerialDenseVector& Ff = matrix->Ff_cells()[c];
//     double& Fc = matrix->Fc_cells()[c];

//     for (int n=0; n!=nfaces; ++n) {
//       int f = faces[n];
//       const AmanziGeometry::Point& normal = S->mesh()->face_normal(f);

//       double outward_flux = ((K_[c] * gravity) * normal)
//         * dirs[n] * (*Krel)("face",0,f) * rho ;
//       Ff[n] += outward_flux;
//       Fc -= outward_flux;  // Nonzero-sum contribution when not upwinding
//     }
//   }
// };


// /* ******************************************************************
//  * Updates global Darcy vector calculated by a discretization method.
//  ****************************************************************** */
// void OverlandFlow::AddGravityFluxesToVector_(const Teuchos::RCP<State>& S,
//                                              const Teuchos::RCP<CompositeVector>& darcy_flux) {

//   double rho = *(S->GetScalarData("density_liquid"));
//   Teuchos::RCP<const Epetra_Vector>   g_vec = S->GetConstantVectorData("gravity");
//   Teuchos::RCP<const CompositeVector> Krel  = S->GetFieldData("rel_perm_faces");
  
//   AmanziGeometry::Point gravity(g_vec->MyLength());
//   for (int i=0; i!=g_vec->MyLength(); ++i) gravity[i] = (*g_vec)[i];

//   AmanziMesh::Entity_ID_List faces;
//   std::vector<int> dirs;

//   int f_used = S->mesh()->count_entities(AmanziMesh::FACE, AmanziMesh::USED);
//   int f_owned = S->mesh()->count_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
//   std::vector<bool> done(f_used, false);

//   int c_owned = S->mesh()->count_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
//   for (int c=0; c!=c_owned; ++c) {
//     S->mesh()->cell_get_faces_and_dirs(c, &faces, &dirs);
//     int nfaces = faces.size();

//     for (int n=0; n!=nfaces; ++n) {
//       int f = faces[n];
//       const AmanziGeometry::Point& normal = S->mesh()->face_normal(f);

//       if (f<f_owned && !done[f]) {
//         (*darcy_flux)(f) += ((K_[c] * gravity) * normal)
//           * (*Krel)("face",0,f) * rho;
//         done[f] = true;
//       }
//     }
//   }
// };

} //namespace
} //namespace
