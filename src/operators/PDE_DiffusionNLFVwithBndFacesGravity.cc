/*
  Operators 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatskiy (dasvyat@lanl.gov)
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <vector>

// Amanzi
#include "DenseMatrix.hh"

// Operators
#include "Op_Face_Cell.hh"
#include "PDE_DiffusionNLFVwithBndFacesGravity.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Populate face-based matrices.
****************************************************************** */
void PDE_DiffusionNLFVwithBndFacesGravity::UpdateMatrices(
    const Teuchos::Ptr<const CompositeVector>& flux,
    const Teuchos::Ptr<const CompositeVector>& u)
{
  // affine map of u. It is equivalent to calculating hydraulic head.
  Teuchos::RCP<CompositeVector> hh = Teuchos::rcp(new CompositeVector(*u));
  Epetra_MultiVector& hh_c = *hh->viewComponent("cell");
  Epetra_MultiVector& hh_bnd = *hh->viewComponent("boundary_face");
  const Epetra_MultiVector& u_c = *u->viewComponent("cell");
  const Epetra_MultiVector& u_bnd = *u->viewComponent("boundary_face");
  const Epetra_MultiVector* rho_c = is_scalar_ ? nullptr :
                              rho_cv_->viewComponent("cell", true).get();

  for (int c = 0; c < ncells_owned; ++c) {
    double rho_g = (is_scalar_ ? rho_ : (*rho_c)[0][c]) * fabs(g_[dim_ - 1]);
    double zc = (mesh_->getCellCentroid(c))[dim_ - 1];
    hh_c[0][c] = u_c[0][c] + rho_g * zc;
    AmanziMesh::Entity_ID_List faces;
    mesh_->cell_get_faces(c, &faces);
    for (auto f : faces) {
      int bf = mesh_->exterior_face_map(false).LID(mesh_->face_map(false).GID(f));
      if (bf >= 0) {
        double zf = (mesh_->getFaceCentroid(f))[dim_ - 1];
        hh_bnd[0][bf] = u_bnd[0][bf] + rho_g*zf;
      }
    }
  }

  if (!is_scalar_) {
    rho_cv_->scatterMasterToGhosted("cell");
  }

  PDE_DiffusionNLFVwithBndFaces::UpdateMatrices(flux, hh.ptr());

  // add gravity fluxes to the right-hand side.
  global_op_->rhs()->PutScalarGhosted(0.0);

  const std::vector<int>& bc_model = bcs_trial_[0]->bc_model();
  Epetra_MultiVector& rhs_cell = *global_op_->rhs()->viewComponent("cell", true);
  Epetra_MultiVector& rhs_bnd = *global_op_->rhs()->viewComponent("boundary_face", true);
  Epetra_MultiVector& weight = *stencil_data_->viewComponent("weight", true);

  AmanziMesh::Entity_ID_List cells;

  for (int f = 0; f < nfaces_owned; ++f) {
    WhetStone::DenseMatrix& Aface = local_op_->matrices[f];
    
    mesh_->face_get_cells(f, AmanziMesh::Parallel_kind::ALL, &cells);
    int ncells = cells.size();

    if (ncells == 2) {
      WhetStone::DenseVector v(ncells), av(ncells);
      for (int n = 0; n < ncells; n++) {
        int c = cells[n];
        double rho_g = (is_scalar_ ? rho_ : (*rho_c)[0][c]) * fabs(g_[dim_ - 1]);
        double zc = (mesh_->getCellCentroid(c))[dim_ - 1];
        v(n) = zc * rho_g;
      }

      Aface.Multiply(v, av, false);

      for (int n = 0; n < ncells; n++) {
        rhs_cell[0][cells[n]] -= av(n);
      }
    } else if ((bc_model[f] == OPERATOR_BC_DIRICHLET)||(bc_model[f] == OPERATOR_BC_NEUMANN)) {
      int c = cells[0];
      double rho_g = (is_scalar_ ? rho_ : (*rho_c)[0][c]) * fabs(g_[dim_ - 1]);
      double zf = (mesh_->getFaceCentroid(f))[dim_ - 1];
      double zc = (mesh_->getCellCentroid(c))[dim_ - 1];
      double gravity_flux = 0.;
      gravity_flux = Aface(0, 0)*zc + Aface(0, 1)*zf;

      //rhs_cell[0][c] -= Aface(0, 0) * (zc - zf) * rho_g;
      rhs_cell[0][c] -=  (Aface(0, 0)*zc + Aface(0, 1)*zf )* rho_g;
      
      int bf = mesh_->exterior_face_map(false).LID(mesh_->face_map(false).GID(f));
      rhs_bnd[0][bf] -=  (Aface(1, 0)*zc + Aface(1, 1)*zf ) * rho_g;

    }
  }

  global_op_->rhs()->gatherGhostedToMaster();

}





/* ******************************************************************
* Calculate flux using cell-centered data.
* **************************************************************** */
void PDE_DiffusionNLFVwithBndFacesGravity::UpdateFlux(const Teuchos::Ptr<const CompositeVector>& u,
                                              const Teuchos::Ptr<CompositeVector>& flux) 
{
  const std::vector<int>& bc_model = bcs_trial_[0]->bc_model();

  // Map field u for the local system. For Richards's equation, this
  // is equivalent to calculating the hydraulic head.
  Teuchos::RCP<CompositeVector> hh = Teuchos::rcp(new CompositeVector(*u));
  Epetra_MultiVector& hh_c = *hh->viewComponent("cell");
  Epetra_MultiVector& hh_bnd = *hh->viewComponent("boundary_face");
  const Epetra_MultiVector& u_c = *u->viewComponent("cell"); 
  const Epetra_MultiVector& u_bnd = *u->viewComponent("boundary_face"); 

  const Epetra_MultiVector* rho_c = is_scalar_ ? nullptr : rho_cv_->viewComponent("cell", false).get();
  
  for (int c = 0; c < ncells_owned; ++c) {
    double rho_g = (is_scalar_ ? rho_ : (*rho_c)[0][c]) * fabs(g_[dim_ - 1]);
    double zc = (mesh_->getCellCentroid(c))[dim_ - 1];
    hh_c[0][c] = u_c[0][c] + rho_g * zc;
    AmanziMesh::Entity_ID_List faces;
    mesh_->cell_get_faces(c, &faces);
    for (auto f : faces) {
      int bf = mesh_->exterior_face_map(false).LID(mesh_->face_map(false).GID(f));
      if (bf >= 0) {
        double zf = (mesh_->getFaceCentroid(f))[dim_ - 1];
        hh_bnd[0][bf] = u_bnd[0][bf] + rho_g*zf;
      }
    }
  }
  PDE_DiffusionNLFVwithBndFaces::UpdateFlux(hh.ptr(), flux);
}


}  // namespace Operators
}  // namespace Amanzi

