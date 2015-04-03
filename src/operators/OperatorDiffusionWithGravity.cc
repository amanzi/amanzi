/*
  This is the operators component of the Amanzi code. 

  Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <vector>

#include "OperatorDefs.hh"
#include "OperatorDiffusionWithGravity.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Add a gravity term to the diffusion operator.
****************************************************************** */
void OperatorDiffusionWithGravity::UpdateMatrices(
    const Teuchos::Ptr<const CompositeVector>& flux,
    const Teuchos::Ptr<const CompositeVector>& u)
{
  OperatorDiffusionMFD::UpdateMatrices(flux, u);
  AddGravityToRHS_();
}


/* ******************************************************************
* Add a gravity term to the RHS of the operator
****************************************************************** */
void OperatorDiffusionWithGravity::AddGravityToRHS_()
{
  // vector or scalar rho?
  const Epetra_MultiVector* rho_c;
  if (rho_cv_ != Teuchos::null) rho_c = &*rho_cv_->ViewComponent("cell", false);

  if (global_op_->rhs()->HasComponent("face")) {
    int dim = mesh_->space_dimension();

    // preparing little-k data
    Teuchos::RCP<const Epetra_MultiVector> k_cell = Teuchos::null;
    Teuchos::RCP<const Epetra_MultiVector> k_face = Teuchos::null;
    if (k_ != Teuchos::null) {
      if (k_->HasComponent("cell")) k_cell = k_->ViewComponent("cell");
      if (k_->HasComponent("face")) k_face = k_->ViewComponent("face", true);
    }

    AmanziMesh::Entity_ID_List faces;
    std::vector<int> dirs;

    Epetra_MultiVector& rhs_cell = *global_op_->rhs()->ViewComponent("cell");
    Epetra_MultiVector& rhs_face = *global_op_->rhs()->ViewComponent("face", true);
    for (int f = nfaces_owned; f < nfaces_wghost; f++) rhs_face[0][f] = 0.0;

    WhetStone::Tensor Kc(mesh_->space_dimension(), 1);
    Kc(0, 0) = 1.0;

    for (int c = 0; c < ncells_owned; c++) {
      mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
      int nfaces = faces.size();
      double zc = (mesh_->cell_centroid(c))[dim - 1];

      // building blocks for the gravity term
      double rho = (rho_cv_ == Teuchos::null) ? rho_ : (*rho_c)[0][c];
      WhetStone::DenseMatrix& Wff = Wff_cells_[c];

      // Update terms due to nonlinear coefficient
      double kc(1.0);
      std::vector<double> kf(nfaces, 1.0);
      if (little_k_ == OPERATOR_LITTLE_K_ARTIFICIAL_DIFFUSION) {
        kc = (*k_cell)[0][c];
        for (int n = 0; n < nfaces; n++) kf[n] = kc;
      } else if (little_k_ == OPERATOR_LITTLE_K_DIVK) {
        kc = (*k_cell)[0][c];
        for (int n = 0; n < nfaces; n++) kf[n] = (*k_face)[0][faces[n]];
      } else if (little_k_ == OPERATOR_LITTLE_K_STANDARD && k_cell != Teuchos::null) {
        kc = (*k_cell)[0][c];
        for (int n = 0; n < nfaces; n++) kf[n] = kc;
      } else if(little_k_ == OPERATOR_LITTLE_K_UPWIND) {
        for (int n = 0; n < nfaces; n++) kf[n] = (*k_face)[0][faces[n]];
      }

      if (little_k_ != OPERATOR_LITTLE_K_DIVK) {
        if (K_.get()) Kc = (*K_)[c];
        AmanziGeometry::Point Kcg(Kc * g_);

        for (int n = 0; n < nfaces; n++) {
          int f = faces[n];
          const AmanziGeometry::Point& normal = mesh_->face_normal(f);
          double tmp, zf = (mesh_->face_centroid(f))[dim - 1];

          if (gravity_special_projection_) {
            const AmanziGeometry::Point& xcc = GravitySpecialDirection_(f);
            double sign = normal * xcc;

            tmp = (Kcg * xcc) * rho * kf[n] * dirs[n];
            tmp *= copysign(norm(normal) / norm(xcc), sign);
          } else {
            tmp = (Kcg * normal) * rho * kf[n] * dirs[n];
          }

          if (little_k_ == OPERATOR_LITTLE_K_ARTIFICIAL_DIFFUSION) {
            double alpha = (*k_face)[0][f] - kc;
            if (alpha > 0) {
              alpha *= Wff(n, n) * rho * norm(g_);
              tmp -= alpha * (zf - zc);
            }
          }

          rhs_face[0][f] += tmp; 
          rhs_cell[0][c] -= tmp; 
        }
      }

      if (little_k_ == OPERATOR_LITTLE_K_DIVK) {
        WhetStone::DenseVector v(nfaces), av(nfaces);
        for (int n = 0; n < nfaces; n++) {
          int f = faces[n];
          double zf = (mesh_->face_centroid(f))[dim - 1];
          v(n) = -(zf - zc) * kf[n] * rho * norm(g_) / kc;
        }

        Wff.Multiply(v, av, false);

        for (int n = 0; n < nfaces; n++) {
          int f = faces[n];
          double tmp = av(n) * kf[n];

          rhs_face[0][f] += tmp; 
          rhs_cell[0][c] -= tmp; 
        }
      }
    }
    global_op_->rhs()->GatherGhostedToMaster("face", Epetra_CombineMode(Add));
  }
}


/* ******************************************************************
* WARNING: Since gravity flux is not continuous, we derive it only once
* (using flag) and in exactly the same manner as in other routines.
* **************************************************************** */
void OperatorDiffusionWithGravity::UpdateFlux(
    const CompositeVector& u, CompositeVector& flux)
{
  // Calculate diffusive part of the flux.
  OperatorDiffusionMFD::UpdateFlux(u, flux);

  // vector or scalar rho?
  const Epetra_MultiVector* rho_c;
  if (rho_cv_ != Teuchos::null) rho_c = &*rho_cv_->ViewComponent("cell", false);

  // preparing little-k data
  Teuchos::RCP<const Epetra_MultiVector> k_cell = Teuchos::null;
  Teuchos::RCP<const Epetra_MultiVector> k_face = Teuchos::null;
  if (k_ != Teuchos::null) {
    if (k_->HasComponent("cell")) k_cell = k_->ViewComponent("cell");
    if (k_->HasComponent("face")) k_face = k_->ViewComponent("face", true);
  }

  int dim = mesh_->space_dimension();
  Epetra_MultiVector& flux_data = *flux.ViewComponent("face", true);

  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;
  std::vector<int> flag(nfaces_wghost, 0);

  WhetStone::Tensor Kc(dim, 1);
  Kc(0, 0) = 1.0;

  for (int c = 0; c < ncells_owned; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();
    double zc = mesh_->cell_centroid(c)[dim - 1];

    // building blocks for the gravity term
    double rho = (rho_cv_ == Teuchos::null) ? rho_ : (*rho_c)[0][c];
    WhetStone::DenseMatrix& Wff = Wff_cells_[c];

    // Update terms due to nonlinear coefficient
    double kc(1.0);
    std::vector<double> kf(nfaces, 1.0);
    if (little_k_ == OPERATOR_LITTLE_K_ARTIFICIAL_DIFFUSION) {
      kc = (*k_cell)[0][c];
      for (int n = 0; n < nfaces; n++) kf[n] = kc;
    } else if (little_k_ == OPERATOR_LITTLE_K_DIVK) {
      kc = (*k_cell)[0][c];
      for (int n = 0; n < nfaces; n++) kf[n] = (*k_face)[0][faces[n]];
    } else if (little_k_ == OPERATOR_LITTLE_K_STANDARD && k_cell != Teuchos::null) {
      kc = (*k_cell)[0][c];
      for (int n = 0; n < nfaces; n++) kf[n] = kc;
    } else if(little_k_ == OPERATOR_LITTLE_K_UPWIND) {
      for (int n = 0; n < nfaces; n++) kf[n] = (*k_face)[0][faces[n]];
    }

    if (little_k_ != OPERATOR_LITTLE_K_DIVK) {
      if (K_.get()) Kc = (*K_)[c];
      AmanziGeometry::Point Kcg(Kc * g_);

      for (int n = 0; n < nfaces; n++) {
        int f = faces[n];
        if (f < nfaces_owned && !flag[f]) {
          const AmanziGeometry::Point& normal = mesh_->face_normal(f);

          if (gravity_special_projection_) {
            const AmanziGeometry::Point& xcc = GravitySpecialDirection_(f);
            double sign = normal * xcc;
            double tmp = copysign(norm(normal) / norm(xcc), sign);
            flux_data[0][f] += (Kcg * xcc) * rho * kf[n] * tmp;
          } else {
            flux_data[0][f] += (Kcg * normal) * rho * kf[n];
          }
            
          if (little_k_ == OPERATOR_LITTLE_K_ARTIFICIAL_DIFFUSION) {
            double alpha = (*k_face)[0][f] - kc;
            if (alpha > 0) {
              int dir;
              const AmanziGeometry::Point& exterior = mesh_->face_normal(f, false, c, &dir);
              double zf = mesh_->face_centroid(f)[dim - 1];
              alpha *= Wff(n, n) * rho * norm(g_);
              flux_data[0][f] -= alpha * (zf - zc) * dir;
            }
          }

          flag[f] = 1;
        }
      }
    }

    if (little_k_ == OPERATOR_LITTLE_K_DIVK) {
      WhetStone::DenseVector v(nfaces), av(nfaces);
      for (int n = 0; n < nfaces; n++) {
        int f = faces[n];
        double zf = (mesh_->face_centroid(f))[dim - 1];
        v(n) = -(zf - zc) * kf[n] * rho * norm(g_) / kc;
      }

      Wff.Multiply(v, av, false);

      for (int n = 0; n < nfaces; n++) {
        int dir, f = faces[n];
        if (f < nfaces_owned && !flag[f]) {
          const AmanziGeometry::Point& normal = mesh_->face_normal(f, false, c, &dir);
            
          double tmp = av(n) * kf[n] * dir;
          flux_data[0][f] += tmp;

          flag[f] = 1;
        }
      }
    }
  }
}


/* ******************************************************************
* Compute non-normalized unsigned direction to the next cell needed
* to project gravity vector in the MFD-TPFA discretization method.
****************************************************************** */
AmanziGeometry::Point OperatorDiffusionWithGravity::GravitySpecialDirection_(int f) const
{
  AmanziMesh::Entity_ID_List cells;
  mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
  int ncells = cells.size();

  if (ncells == 2) {
    return mesh_->cell_centroid(cells[1]) - mesh_->cell_centroid(cells[0]);
  } else {
    return mesh_->face_centroid(f) - mesh_->cell_centroid(cells[0]);
  }
}

}  // namespace Operators
}  // namespace Amanzi
