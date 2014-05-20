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
void OperatorDiffusionWithGravity::UpdateMatrices(Teuchos::RCP<const CompositeVector> flux)
{
  // add the diffusion matrices
  OperatorDiffusion::UpdateMatrices(flux);

  // add the gravity terms
  AmanziGeometry::Point rho_g(g_);
  rho_g *= rho_ * rho_ / mu_;

  if (rhs_->HasComponent("face")) {
    Epetra_MultiVector& rhs = *rhs_->ViewComponent("face", true);

    AmanziMesh::Entity_ID_List faces;
    std::vector<int> dirs;

    for (int f = nfaces_owned; f < nfaces_wghost; f++) rhs[0][f] = 0.0;

    for (int c = 0; c < ncells_owned; c++) {
      mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
      int nfaces = faces.size();

      WhetStone::Tensor& Kc = (*K_)[c]; 
      for (int n = 0; n < nfaces; n++) {
        int f = faces[n];
        const AmanziGeometry::Point& normal = mesh_->face_normal(f);
        rhs[0][f] += ((Kc * rho_g) * normal) * dirs[n]; 
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

    WhetStone::Tensor& Kc = (*K_)[c];
    AmanziGeometry::Point Kg(Kc * rho_g);
    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];
      if (f < nfaces_owned && !flag[f]) {
        const AmanziGeometry::Point& normal = mesh_->face_normal(f);
        flux_data[0][f] += (Kg * normal);
        flag[f] = 1;
      }
    }
  }
}

}  // namespace Operators
}  // namespace Amanzi
