/*
  Operators 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
          Ethan Coon (ecoon@lanl.gov)
*/

#include <vector>

#include "WhetStoneDefs.hh"

#include "OperatorDefs.hh"
#include "Operator_Cell.hh"
#include "Op_Face_Cell.hh"
#include "Op_SurfaceFace_SurfaceCell.hh"
#include "PDE_AdvectionUpwindFracture.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Advection requires a velocity field.
****************************************************************** */
void PDE_AdvectionUpwindFracture::Setup(const CompositeVector& u)
{
  IdentifyUpwindCells_(u);
}

  
/* ******************************************************************
* A simple first-order transport method.
* Advection operator is of the form: div (u C), where u is the given
* velocity field and C is the advected field.
****************************************************************** */
void PDE_AdvectionUpwindFracture::UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& flux)
{
  std::vector<WhetStone::DenseMatrix>& matrix = local_op_->matrices;
  std::vector<WhetStone::DenseMatrix>& matrix_shadow = local_op_->matrices_shadow;

  AmanziMesh::Entity_ID_List cells;
  int c;
  double u;
 
  for (int f = 0; f < nfaces_owned; ++f) {
   
    mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
    int ncells = cells.size();
    WhetStone::DenseMatrix Aface(ncells, ncells);
    Aface.PutScalar(0.0);

    double flux_in(0.0);
    std::vector<int> upwind_loc(upwind_flux_[f].size());

    
    for (int n=0; n<upwind_flux_[f].size(); n++){
      c = upwind_cells_[f][n];
      u = upwind_flux_[f][n];      
      for (int j=0; j<cells.size(); j++){
        if (cells[j]==c){
          upwind_loc[n] = j;
          Aface(j,j) = u;
          break;
        }
      }      
    }

    for (int n = 0; n < downwind_cells_[f].size(); ++n) {
      flux_in -= downwind_flux_[f][n];
    }
    if (flux_in == 0.0) flux_in = 1e-12;
    
    for (int n = 0; n < downwind_cells_[f].size(); ++n) {
      int c = downwind_cells_[f][n];
      u = downwind_flux_[f][n];

      
      if (c < ncells_owned) {
        double tmp = u / flux_in;
        for (int m=0; m<upwind_flux_[f].size(); m++){
          double v = upwind_flux_[f][m];
          for (int j=0; j<cells.size(); j++){
            if (cells[j]==c){
              Aface(j, upwind_loc[m]) = (u / flux_in)*v;
              break;
            }
          }
        }
      }
    }

  }

  
}



/* *******************************************************************
* Apply boundary condition to the local matrices
*
* Advection only problem.
* Recommended options: primary=true, eliminate=false, essential_eqn=true
*  - must deal with Dirichlet BC on inflow boundary
*  - Dirichlet on outflow boundary is ill-posed
*  - Neumann on inflow boundary is typically not used, since it is 
*    equivalent to Dirichlet BC. We perform implicit conversion to 
*    Dirichlet BC.
*
* Advection-diffusion problem.
* Recommended options: primary=false, eliminate=true, essential_eqn=true
*  - Dirichlet BC is treated as usual
*  - Neuman on inflow boundary: If diffusion takes care of the total
*    flux, then TOTAL_FLUX model must be used. If diffusion deals
*    with the diffusive flux only (NEUMANN model), value of the 
*    advective flux is in general not available and negative value
*    is added to matrix diagonal. The discrete system may lose SPD 
*    property.
*  - Neuman on outflow boundary: If diffusion takes care of the total
*    flux, then TOTAL_FLUX model must be used. Otherwise, do nothing.
*
* FIXME: So far we support the case bc_test = bc_trial
******************************************************************* */
// void PDE_AdvectionUpwindFracture::ApplyBCs(bool primary, bool eliminate, bool essential_eqn)
// {
//   std::vector<WhetStone::DenseMatrix>& matrix = local_op_->matrices;
//   std::vector<WhetStone::DenseMatrix>& matrix_shadow = local_op_->matrices_shadow;

//   Epetra_MultiVector& rhs_cell = *global_op_->rhs()->ViewComponent("cell");

//   const std::vector<int>& bc_model = bcs_trial_[0]->bc_model();
//   const std::vector<double>& bc_value = bcs_trial_[0]->bc_value();

//   for (int f = 0; f < nfaces_owned; f++) {
//     int c1 = (*upwind_cell_)[f];
//     int c2 = (*downwind_cell_)[f];
//     if (bc_model[f] == OPERATOR_BC_DIRICHLET) {
//       if (c2 < 0) {
//         // pass, the upwind cell is internal to the domain, so all is good
//       } else if (c1 < 0) {
//         // downwind cell is internal to the domain
//         rhs_cell[0][c2] += matrix[f](0, 0) * bc_value[f];
//         matrix[f] = 0.0;
//       }
//     } 

//     // treat as essential inflow BC for pure advection
//     else if (bc_model[f] == OPERATOR_BC_NEUMANN && primary) {
//       if (c1 < 0) {
//         rhs_cell[0][c2] += mesh_->face_area(f) * bc_value[f];
//         matrix[f] = 0.0;
//       }
//     }
//     // leave in matrix for composite operator
//     else if (bc_model[f] == OPERATOR_BC_NEUMANN && ! primary) {
//       if (c1 < 0)
//         matrix[f] *= -1.0;
//     }
//     // total flux was processed by another operator -> remove here
//     else if (bc_model[f] == OPERATOR_BC_TOTAL_FLUX && ! primary) {
//       matrix[f] = 0.0;
//     }
//     // do not know what to do
//     else if (bc_model[f] != OPERATOR_BC_NONE) {
//       AMANZI_ASSERT(false);
//     } 
//   }
// }



/* *******************************************************************
* Identify flux direction based on orientation of the face normal 
* and sign of the  Darcy velocity.                               
******************************************************************* */
void PDE_AdvectionUpwindFracture::IdentifyUpwindCells_(const CompositeVector& u)
{

  upwind_cells_.clear();
  downwind_cells_.clear();

  upwind_cells_.resize(nfaces_wghost);
  downwind_cells_.resize(nfaces_wghost);

  upwind_flux_.clear();
  downwind_flux_.clear();

  upwind_flux_.resize(nfaces_wghost);
  downwind_flux_.resize(nfaces_wghost);

  AmanziMesh::Entity_ID_List faces, cells;
  std::vector<int> dirs;

  const Epetra_Map& cmap = mesh_->cell_map(true);

 // desinged for $domain-darcy_flux_fracture

  u.ScatterMasterToGhosted();

  const Epetra_MultiVector& flux_c = *u.ViewComponent("cell");
  
  for (int c = 0; c < ncells_wghost; c++) {
    mesh_->cell_get_faces(c, &faces);
    
    for (int i = 0; i < faces.size(); i++) {
      int f = faces[i];
      double u = flux_c[i][c];
      if (u >= 0.0) {
        upwind_cells_[f].push_back(c);
        upwind_flux_[f].push_back(u);
      } else {
        downwind_cells_[f].push_back(c);
        downwind_flux_[f].push_back(u);
      }      
    }
  }
  
}

}  // namespace Operators
}  // namespace Amanzi
