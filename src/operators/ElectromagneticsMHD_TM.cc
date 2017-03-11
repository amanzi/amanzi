/*
  Operators 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <vector>

// TPLs
#include "Epetra_Vector.h"

// Amanzi::Operators
#include "ElectromagneticsMHD_TM.hh"
#include "Op.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* System modification before solving the problem.
* **************************************************************** */
void ElectromagneticsMHD_TM::ModifyMatrices(
   CompositeVector& E, CompositeVector& B, double dt)
{
  const Epetra_MultiVector& Bf = *B.ViewComponent("face", true);
  Epetra_MultiVector& rhs_v = *global_op_->rhs()->ViewComponent("node", true);

  std::vector<int> dirs;
  AmanziMesh::Entity_ID_List faces, nodes;

  for (int c = 0; c < ncells_owned; ++c) {
    WhetStone::DenseMatrix& Acell = local_op_->matrices[c];
    Acell.Scale(dt / 2);

    const WhetStone::DenseMatrix& Mcell = mass_op_[c];
    const WhetStone::DenseMatrix& Ccell = curl_op_[c];

    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    mesh_->cell_get_nodes(c, &nodes);

    int nfaces = faces.size();
    int nnodes = nodes.size();

    WhetStone::DenseVector v1(nfaces), v2(nfaces), v3(nnodes);

    for (int n = 0; n < nfaces; ++n) {
      int f = faces[n];
      v1(n) = Bf[0][f] * dirs[n] * mesh_->face_area(f);
    }

    Mcell.Multiply(v1, v2, false);
    Ccell.Multiply(v2, v3, true);

    for (int n = 0; n < nnodes; ++n) {
      int v = nodes[n];
      rhs_v[0][v] += v3(n);
    }
  }
}


/* ******************************************************************
* Solution postprocessing
* **************************************************************** */
void ElectromagneticsMHD_TM::ModifyFields(
   CompositeVector& E, CompositeVector& B, double dt)
{
  B.ScatterMasterToGhosted("face");

  Epetra_MultiVector& Ev = *E.ViewComponent("node", true);
  Epetra_MultiVector& Bf = *B.ViewComponent("face", false);
  
  std::vector<int> dirs;
  AmanziMesh::Entity_ID_List faces, nodes;

  std::vector<bool> fflag(nedges_wghost, false);

  for (int c = 0; c < ncells_owned; ++c) {
    const WhetStone::DenseMatrix& Ccell = curl_op_[c];

    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    mesh_->cell_get_nodes(c, &nodes);

    int nfaces = faces.size();
    int nnodes = nodes.size();

    WhetStone::DenseVector v1(nnodes), v2(nfaces);

    for (int n = 0; n < nnodes; ++n) {
      int v = nodes[n];
      v1(n) = Ev[0][v];
    }

    Ccell.Multiply(v1, v2, false);

    for (int n = 0; n < nfaces; ++n) {
      int f = faces[n];
      if (!fflag[f]) {
        Bf[0][f] -= dt * v2(n) * dirs[n] / mesh_->face_area(f);
        fflag[f] = true;
      }
    }
  }
}

}  // namespace Operators
}  // namespace Amanzi
