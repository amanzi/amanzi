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

}  // namespace Operators
}  // namespace Amanzi
