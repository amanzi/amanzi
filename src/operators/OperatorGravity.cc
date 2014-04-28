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
#include "OperatorGravity.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Add a gravity term to the operator
****************************************************************** */
void OperatorGravity::UpdateMatrices(
    const std::vector<WhetStone::Tensor>& K, const AmanziGeometry::Point& rho_g)
{
  if (rhs_->HasComponent("face")) {
    Epetra_MultiVector& rhs = *rhs_->ViewComponent("face", true);

    AmanziMesh::Entity_ID_List faces;
    std::vector<int> dirs;

    for (int c = 0; c < ncells_owned; c++) {
      mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
      int nfaces = faces.size();

      for (int n = 0; n < nfaces; n++) {
        int f = faces[n];
        const AmanziGeometry::Point& normal = mesh_->face_normal(f);
        rhs[0][f] += ((K[c] * rho_g) * normal) * dirs[n]; 
      }
    }
  }
}


/* ******************************************************************
* WARNING: Since gavity flux is not continuous, we derive it only once
* (using flag) and in exactly the same manner as in other routines.
* **************************************************************** */
void OperatorGravity::UpdateFlux(
    const std::vector<WhetStone::Tensor>& K, const AmanziGeometry::Point& rho_g, 
    CompositeVector& flux, double scalar)
{
  // Initialize the flux in the case of additive operators.
  if (scalar == 0.0) flux.PutScalar(0.0);
  Epetra_MultiVector& flux_data = *flux.ViewComponent("face", true);

  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;
  std::vector<int> flag(nfaces_wghost, 0);

  for (int c = 0; c < ncells_owned; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];
      if (f < nfaces_owned && !flag[f]) {
        const AmanziGeometry::Point& normal = mesh_->face_normal(f);
        flux_data[0][f] += ((K[c] * rho_g) * normal) * dirs[n]; 
        flag[f] = 1;
      }
    }
  }
}

}  // namespace Operators
}  // namespace Amanzi
