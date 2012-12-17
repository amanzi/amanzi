/*
This is the flow component of the Amanzi code. 

Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided in the top-level COPYRIGHT file.

Authors: Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
*/

#include <vector>

#include "Epetra_FECrsGraph.h"

#include "Flow_constants.hpp"
#include "Matrix_MFD_PLambda.hpp"

namespace Amanzi {
namespace AmanziFlow {

/* ******************************************************************
* Initialize Trilinos matrices. It must be called only once. 
* If matrix is non-symmetric, we generate transpose of the matrix 
* block Afc to reuse cf_graph; otherwise, pointer Afc = Acf.   
****************************************************************** */
void Matrix_MFD_PLambda::SymbolicAssembleGlobalMatrices(const Epetra_Map& super_map)
{
  const Epetra_Map& cmap = mesh_->cell_map(false);
  const Epetra_Map& fmap = mesh_->face_map(false);
  const Epetra_Map& fmap_wghost = mesh_->face_map(true);

  int avg_entries_row = (mesh_->space_dimension() == 2) ? FLOW_QUAD_FACES : FLOW_HEX_FACES;
  Epetra_FECrsGraph graph(Copy, super_map, 2*avg_entries_row);

  AmanziMesh::Entity_ID_List faces, cells;
  std::vector<int> dirs;
  int dof_LID[FLOW_MAX_FACES + 1];  // Contigious memory is required.
  int dof_GID[FLOW_MAX_FACES + 1];

  // diffusion part
  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c = 0; c < ncells; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();
    int ndof = nfaces + 1;

    for (int n = 0; n < nfaces; n++) dof_LID[n] = ncells + faces[n];
    dof_LID[nfaces] = c;
    for (int n = 0; n < ndof; n++) dof_GID[n] = super_map.GID(dof_LID[n]);
  
    graph.InsertGlobalIndices(ndof, dof_GID, ndof, dof_GID);
  }

  // advection part
  int nfaces = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  for (int f = 0; f < nfaces; f++) {
    mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
    int ncells = cells.size();

    for (int n = 0; n < ncells; n++) dof_GID[n] = cells[n];

    graph.InsertGlobalIndices(ncells, dof_GID, ncells, dof_GID);
  }
  graph.GlobalAssemble();  // Symbolic graph is complete.

  // create global matrices
  APLambda_ = Teuchos::rcp(new Epetra_FECrsMatrix(Copy, graph));
  APLambda_->GlobalAssemble();

  rhs_ = Teuchos::rcp(new Epetra_Vector(super_map));
  rhs_cells_ = Teuchos::rcp(FS->CreateCellView(*rhs_));
  rhs_faces_ = Teuchos::rcp(FS->CreateFaceView(*rhs_));

  // set a preconditioner
  Sff_ = APLambda_;
}


/* ******************************************************************
* Assemble elemental mass matrices into one global matrix. 
* We need an auxiliary GHOST-based vector to assemble the RHS.
****************************************************************** */
void Matrix_MFD_PLambda::AssembleGlobalMatrices()
{
  APLambda_->PutScalar(0.0);
  const Epetra_Map& Amap = APLambda_->RowMap();

  int dof_LID[FLOW_MAX_FACES + 1];
  int dof_GID[FLOW_MAX_FACES + 1];

  // diffusion part
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;
  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);

  for (int c = 0; c < ncells; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();
    int ndof = nfaces + 1;

    for (int n = 0; n < nfaces; n++) dof_LID[n] = ncells + faces[n];
    dof_LID[nfaces] = c;
    for (int n = 0; n < ndof; n++) dof_GID[n] = Amap.GID(dof_LID[n]);

    Teuchos::SerialDenseMatrix<int, double> BPLambda(ndof, ndof);
    for (int n = 0; n < nfaces; n++)
      for (int m = 0; m < nfaces; m++) BPLambda(m, n) = Aff_cells_[c](m, n);

    for (int n = 0; n < nfaces; n++) {
      BPLambda(nfaces, n) = Acf_cells_[c][n];
      BPLambda(n, nfaces) = Afc_cells_[c][n];
    }
    BPLambda(nfaces, nfaces) = Acc_cells_[c];

    APLambda_->SumIntoGlobalValues(ndof, dof_GID, BPLambda.values());
  }

  // advection part
  int nfaces = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);

  if (Acc_faces_.size() > 0) {
    AmanziMesh::Entity_ID_List cells;
    for (int f = 0; f < nfaces; f++) {
      mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
      int ncells = cells.size();

      for (int n = 0; n < ncells; n++) dof_GID[n] = cells[n];

      Teuchos::SerialDenseMatrix<int, double>& Bcc = Acc_faces_[f];
      APLambda_->SumIntoGlobalValues(ncells, dof_GID, Bcc.values());
    }
  }
  APLambda_->GlobalAssemble();

  // We repeat some of the loops for code clarity.
  const Epetra_Map& fmap_wghost = mesh_->face_map(true);
  Epetra_Vector rhs_faces_wghost(fmap_wghost);

  for (int c = 0; c < ncells; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int mfaces = faces.size();

    (*rhs_cells_)[c] = Fc_cells_[c];
    for (int n = 0; n < mfaces; n++) {
      int f = faces[n];
      rhs_faces_wghost[f] += Ff_cells_[c][n];
    }
  }
  FS->CombineGhostFace2MasterFace(rhs_faces_wghost, Add);

  for (int f = 0; f < nfaces; f++) (*rhs_faces_)[f] = rhs_faces_wghost[f];
}


/* ******************************************************************
* Parallel matvec product A * X.                                              
****************************************************************** */
int Matrix_MFD_PLambda::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  APLambda_->Multiply(false, X, Y);
  return 0;
}


/* ******************************************************************
* Parallel matvec product inv(Prec) * X. 
* For some unknwon reasons, interface with Hypre requires to 
* re-define arguments of ApplyInverse().  
****************************************************************** */
int Matrix_MFD_PLambda::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  if (method_ == FLOW_PRECONDITIONER_TRILINOS_ML) {
    MLprec->ApplyInverse(X, Y);
  } else if (method_ == FLOW_PRECONDITIONER_HYPRE_AMG) { 
#ifdef HAVE_HYPRE
    Epetra_MultiVector XX(X);
    Epetra_MultiVector YY(Y);
    IfpHypre_Sff_->ApplyInverse(XX, YY);
    Y = YY;
#endif
  } else if (method_ == FLOW_PRECONDITIONER_TRILINOS_BLOCK_ILU) {
    ifp_prec_->ApplyInverse(X, Y);
  }
  return 0;
}

}  // namespace AmanziFlow
}  // namespace Amanzi




