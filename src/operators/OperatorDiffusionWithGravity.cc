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
void OperatorDiffusionWithGravity::UpdateMatrices(Teuchos::RCP<const CompositeVector> flux,
                                                  Teuchos::RCP<const CompositeVector> u)
{
  // add the diffusion matrices
  OperatorDiffusion::UpdateMatrices(flux, u);

  // add the gravity terms
  AmanziGeometry::Point rho_g(g_);
  rho_g *= rho_ * rho_ / mu_;

  if (rhs_->HasComponent("face")) {
    int dim = mesh_->space_dimension();

    // preparing upwind data
    Teuchos::RCP<const Epetra_MultiVector> k_cell = Teuchos::null;
    Teuchos::RCP<const Epetra_MultiVector> k_face = Teuchos::null;
    if (k_ != Teuchos::null) k_cell = k_->ViewComponent("cell");
    if (upwind_ == OPERATOR_UPWIND_WITH_FLUX || upwind_ == OPERATOR_UPWIND_AMANZI) {
      k_face = k_->ViewComponent("face", true);
    }

    AmanziMesh::Entity_ID_List faces;
    std::vector<int> dirs;

    Epetra_MultiVector& rhs_cell = *rhs_->ViewComponent("cell");
    Epetra_MultiVector& rhs_face = *rhs_->ViewComponent("face", true);
    for (int f = nfaces_owned; f < nfaces_wghost; f++) rhs_face[0][f] = 0.0;

    for (int c = 0; c < ncells_owned; c++) {
      mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
      int nfaces = faces.size();

      // Update terms due to nonlinear coefficient
      double kc(1.0);
      std::vector<double> kf(nfaces, 1.0);
      if (upwind_ == OPERATOR_UPWIND_AMANZI) {
        kc = (*k_cell)[0][c];
        for (int n = 0; n < nfaces; n++) kf[n] = std::max(kc, (*k_face)[0][faces[n]]);
      } else if (upwind_ == OPERATOR_UPWIND_NONE && k_cell != Teuchos::null) {
        kc = (*k_cell)[0][c];
        for (int n = 0; n < nfaces; n++) kf[n] = kc;
      } else if(upwind_ == OPERATOR_UPWIND_WITH_FLUX) {
        for (int n = 0; n < nfaces; n++) kf[n] = (*k_face)[0][faces[n]];
      }

      WhetStone::Tensor& Kc = (*K_)[c]; 
      for (int n = 0; n < nfaces; n++) {
        int f = faces[n];
        const AmanziGeometry::Point& normal = mesh_->face_normal(f);
        double tmp = ((Kc * rho_g) * normal) * kf[n] * dirs[n];
        rhs_face[0][f] += tmp; 
        rhs_cell[0][c] -= tmp; 
      }
    }

    rhs_->GatherGhostedToMaster("face", Epetra_CombineMode(Add));
  }
}


/* ******************************************************************
* WARNING: Since gavity flux is not continuous, we derive it only once
* (using flag) and in exactly the same manner as in other routines.
* **************************************************************** */
void OperatorDiffusionWithGravity::UpdateFlux(
    const CompositeVector& u, CompositeVector& flux)
{
  // Calculate diffusive part of the flux.
  OperatorDiffusion::UpdateFlux(u, flux);

  // preparing upwind data
  Teuchos::RCP<const Epetra_MultiVector> k_cell = Teuchos::null;
  Teuchos::RCP<const Epetra_MultiVector> k_face = Teuchos::null;
  if (k_ != Teuchos::null) k_cell = k_->ViewComponent("cell");
  if (upwind_ == OPERATOR_UPWIND_WITH_FLUX || upwind_ == OPERATOR_UPWIND_AMANZI) {
    k_face = k_->ViewComponent("face", true);
  }

  // Add gravity part to the flux.
  Epetra_MultiVector& flux_data = *flux.ViewComponent("face", true);
  AmanziGeometry::Point rho_g(g_);
  rho_g *= rho_ * rho_ / mu_;

  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;
  std::vector<int> flag(nfaces_wghost, 0);

  for (int c = 0; c < ncells_owned; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    // Update terms due to nonlinear coefficient
    double kc(1.0);
    std::vector<double> kf(nfaces, 1.0);
    if (upwind_ == OPERATOR_UPWIND_AMANZI) {
      // for (int n = 0; n < nfaces; n++) kf[n] = (*k_face)[0][faces[n]];
      for (int n = 0; n < nfaces; n++) kf[n] = std::max((*k_cell)[0][c], (*k_face)[0][faces[n]]);
    } else if (upwind_ == OPERATOR_UPWIND_NONE && k_cell != Teuchos::null) {
      kc = (*k_cell)[0][c];
      for (int n = 0; n < nfaces; n++) kf[n] = kc;
    } else if(upwind_ == OPERATOR_UPWIND_WITH_FLUX) {
      for (int n = 0; n < nfaces; n++) kf[n] = (*k_face)[0][faces[n]];
    }

    WhetStone::Tensor& Kc = (*K_)[c];
    AmanziGeometry::Point Kg(Kc * rho_g);
    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];
      if (f < nfaces_owned && !flag[f]) {
        const AmanziGeometry::Point& normal = mesh_->face_normal(f);
        flux_data[0][f] += (Kg * normal) * kf[n];
        flag[f] = 1;
      }
    }
  }
}

}  // namespace Operators
}  // namespace Amanzi
