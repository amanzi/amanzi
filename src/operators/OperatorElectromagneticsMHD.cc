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

// Amanzi
#include "errors.hh"
#include "MatrixFE.hh"
#include "mfd3d_electromagnetics.hh"
#include "PreconditionerFactory.hh"
#include "SuperMap.hh"
#include "WhetStoneDefs.hh"

// Operators
#include "Op.hh"
#include "Op_Cell_Edge.hh"
#include "Op_Cell_Face.hh"
#include "OperatorDefs.hh"
#include "Operator_Edge.hh"

#include "OperatorElectromagneticsMHD.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Calculate elemental matrices.
****************************************************************** */
void OperatorElectromagneticsMHD::UpdateMatrices()
{
  std::vector<int> dirs;
  AmanziMesh::Entity_ID_List edges, faces;
  WhetStone::MFD3D_Electromagnetics mfd(mesh_);

  WhetStone::Tensor Kc(mesh_->space_dimension(), 1);
  Kc(0, 0) = 1.0;
  
  for (int c = 0; c < ncells_owned; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    mesh_->cell_get_edges(c, &edges);
    int nedges = edges.size();

    WhetStone::DenseMatrix Mcell(nfaces, nfaces);
    WhetStone::DenseMatrix Ccell(nfaces, nedges);
    WhetStone::DenseMatrix Acell(nedges, nedges);

    if (K_.get()) Kc = (*K_)[c];
    mfd.StiffnessMatrix(c, Kc, Acell, Mcell, Ccell);

    local_op_->matrices[c] = Acell;
    mass_op_[c] = Mcell;
    curl_op_[c] = Ccell;
  }
}


/* ******************************************************************
* System modification before solving the problem.
* **************************************************************** */
void OperatorElectromagneticsMHD::ModifyMatrices(CompositeVector& E,
                                                 CompositeVector& B)
{
  const Epetra_MultiVector& Bf = *B.ViewComponent("face", true);

  std::vector<int> dirs;
  AmanziMesh::Entity_ID_List faces, edges;

  for (int c = 0; c < ncells_owned; ++c) {
    const WhetStone::DenseMatrix& Mcell = mass_op_[c];
    const WhetStone::DenseMatrix& Ccell = curl_op_[c];

    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    mesh_->cell_get_edges(c, &edges);

    int nfaces = faces.size();
    int nedges = edges.size();

    WhetStone::DenseVector v(nfaces), mv(nfaces);

    for (int n = 0; n < nfaces; ++n) {
      int f = faces[n];
      v(n) = Bf[0][f];
    }

    Mcell.Multiply(v, mv, false);
  }
}


/* ******************************************************************
* Solution postprocessing
* **************************************************************** */
void OperatorElectromagneticsMHD::ModifyFields(CompositeVector& E,
                                               CompositeVector& B)
{
  B.PutScalar(0.0);
}


/* ******************************************************************
* Put here stuff that has to be done in constructor.
****************************************************************** */
void OperatorElectromagneticsMHD::InitElectromagneticsMHD_()
{
  mass_op_.resize(ncells_owned);
  curl_op_.resize(ncells_owned);
}

}  // namespace Operators
}  // namespace Amanzi
