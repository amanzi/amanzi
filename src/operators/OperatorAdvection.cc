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
#include "OperatorAdvection.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Initialization
****************************************************************** */
void OperatorAdvection::InitOperator(const CompositeVector& u)
{
  IdentifyUpwindCells_(u);
}

  
/* ******************************************************************
* Add a simple first-order transport method 
****************************************************************** */
void OperatorAdvection::UpdateMatrices(const CompositeVector& u)
{
  // find the block of matrices
  bool flag(false);
  int m, nblocks = matrix_blocks_type_.size();
  for (int n = 0; n < nblocks; n++) {
    int type = matrix_blocks_type_[n];
    if (type == OPERATOR_STENCIL_TYPE_FACE_TPFA) {
      m = n;
      flag = true;
      break;
    }
  }

  if (flag == false) { 
    m = nblocks++;
    matrix_blocks_type_.push_back(OPERATOR_STENCIL_TYPE_FACE_TPFA);
    matrix_blocks_.push_back(Teuchos::rcp(new std::vector<WhetStone::DenseMatrix>));
  }
  std::vector<WhetStone::DenseMatrix>& matrix = *matrix_blocks_[m];

  // apply preconditioner inversion
  AmanziMesh::Entity_ID_List cells;
  const Epetra_MultiVector& uf = *u.ViewComponent("face");

  for (int f = 0; f < nfaces_owned; f++) {  // loop over master and slave faces
    int c1 = (*upwind_cell_)[f];
    int c2 = (*downwind_cell_)[f];

    mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
    int ncells = cells.size();
    WhetStone::DenseMatrix Aface(ncells, ncells);
    Aface.PutScalar(0.0);

    double umod = fabs(uf[0][f]);
    if (c1 < 0) {
      Aface(0, 0) = -umod;
    } else if (c2 < 0) {
      Aface(0, 0) = umod;
    } else {
      int i = (cells[0] == c1) ? 0 : 1;
      Aface(i, i) = umod;
      Aface(1 - i, i) = -umod;
    }

    if (flag) {
      matrix[f] += Aface;
    } else {
      matrix.push_back(Aface);
    }
  }
}


/* *******************************************************************
* Identify flux direction based on orientation of the face normal 
* and sign of the  Darcy velocity.                               
******************************************************************* */
void OperatorAdvection::IdentifyUpwindCells_(const CompositeVector& u)
{
  const Epetra_MultiVector& uf = *u.ViewComponent("face", true);

  const Epetra_Map& fmap_wghost = mesh_->face_map(true);
  upwind_cell_ = Teuchos::rcp(new Epetra_IntVector(fmap_wghost));
  downwind_cell_ = Teuchos::rcp(new Epetra_IntVector(fmap_wghost));

  for (int f = 0; f < nfaces_wghost; f++) {
    (*upwind_cell_)[f] = -1;  // negative value indicates boundary
    (*downwind_cell_)[f] = -1;
  }

  AmanziMesh::Entity_ID_List faces;
  std::vector<int> fdirs;

  for (int c = 0; c < ncells_wghost; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &fdirs);

    for (int i = 0; i < faces.size(); i++) {
      int f = faces[i];
      if (uf[0][f] * fdirs[i] >= 0) {
        (*upwind_cell_)[f] = c;
      } else {
        (*downwind_cell_)[f] = c;
      }
    }
  }
}

}  // namespace Operators
}  // namespace Amanzi
