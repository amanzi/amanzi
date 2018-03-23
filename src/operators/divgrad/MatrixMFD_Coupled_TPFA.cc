/*
  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)

 */

#include "Epetra_FECrsGraph.h"
#include "EpetraExt_RowMatrixOut.h"

#include "errors.hh"
#include "MatrixMFD.hh"

#include "MatrixMFD_Coupled_TPFA.hh"

namespace Amanzi {
namespace Operators {

MatrixMFD_Coupled_TPFA::MatrixMFD_Coupled_TPFA(Teuchos::ParameterList& plist,
        const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) :
  MatrixMFD_Coupled(plist,mesh) {}

MatrixMFD_Coupled_TPFA::MatrixMFD_Coupled_TPFA(const MatrixMFD_Coupled_TPFA& other) :
  MatrixMFD_Coupled(other) {}

void MatrixMFD_Coupled_TPFA::SetSubBlocks(const Teuchos::RCP<MatrixMFD>& blockA,
				     const Teuchos::RCP<MatrixMFD>& blockB) {
  MatrixMFD_Coupled::SetSubBlocks(blockA,blockB);
  blockA_TPFA_ = Teuchos::rcp_dynamic_cast<MatrixMFD_TPFA>(blockA);
  ASSERT(blockA_TPFA_ != Teuchos::null);
  blockB_TPFA_ = Teuchos::rcp_dynamic_cast<MatrixMFD_TPFA>(blockB);
  ASSERT(blockB_TPFA_ != Teuchos::null);
}

int MatrixMFD_Coupled_TPFA::ApplyInverse(const TreeVector& X,
        TreeVector& Y) const {
  if (!assembled_schur_) {
    AssembleSchur_();
    UpdatePreconditioner_();
  }

  if (S_pc_ == Teuchos::null) {
    Errors::Message msg("MatrixMFD::ApplyInverse called but no preconditioner sublist was provided");
    Exceptions::amanzi_throw(msg);
  }

  // pull X data
  Teuchos::RCP<const CompositeVector> XA = X.SubVector(0)->Data();
  Teuchos::RCP<const CompositeVector> XB = X.SubVector(1)->Data();
  const Epetra_MultiVector& XA_c = *XA->ViewComponent("cell", false);
  const Epetra_MultiVector& XB_c = *XB->ViewComponent("cell", false);

  Teuchos::RCP<CompositeVector> YA = Y.SubVector(0)->Data();
  Teuchos::RCP<CompositeVector> YB = Y.SubVector(1)->Data();
  Epetra_MultiVector& YA_c = *YA->ViewComponent("cell", false);
  Epetra_MultiVector& YB_c = *YB->ViewComponent("cell", false);

  // Temporary cell and face vectors.
  Epetra_MultiVector Xc(*double_cmap_, 1);
  Epetra_MultiVector Yc(*double_cmap_, 1);

  int ierr(0);

  int ncells = XA_c.MyLength();
  double val = 0.;

  for (int c=0; c!=ncells; ++c) {
    Xc[0][2*c] = XA_c[0][c];
    Xc[0][2*c+1] = XB_c[0][c];
  }

  // Solve the Schur complement system Spp * Yc = Xc.
  //  Xc.Print(std::cout);
  //S_pc_->ApplyInverse(Xc, Yc);
  ierr = S_pc_->ApplyInverse(Xc, Yc);
  ASSERT(!ierr);

  for (int c=0; c!=ncells; ++c) {
    YA_c[0][c] = Yc[0][2*c];
    YB_c[0][c] = Yc[0][2*c+1];
  }

  if (ierr) {
    Errors::Message msg("MatrixMFD_TPFA::ApplyInverse has failed in calculating y = inv(A)*x.");
    Exceptions::amanzi_throw(msg);
  }

  // Update the faces
  if (YA->HasComponent("face")) {
    const Epetra_MultiVector& XA_f = *XA->ViewComponent("face", false);
    Epetra_MultiVector& YA_f = *YA->ViewComponent("face", false);
    const Epetra_MultiVector& DffA_f = *blockA_TPFA_->Dff()->ViewComponent("face",false);

    blockA_TPFA_->ApplyAfc(*YA, *YA,0.);
    YA_f.Update(1., XA_f, -1.);

    int nfaces = YA_f.MyLength();
    for (int f=0; f!=nfaces; ++f) {
      YA_f[0][f] /= DffA_f[0][f];
    }
  }

  if (YB->HasComponent("face")) {
    const Epetra_MultiVector& XB_f = *XB->ViewComponent("face", false);
    Epetra_MultiVector& YB_f = *YB->ViewComponent("face", false);
    const Epetra_MultiVector& DffB_f = *blockB_TPFA_->Dff()->ViewComponent("face",false);

    blockB_TPFA_->ApplyAfc(*YB, *YB,0.);
    YB_f.Update(1., XB_f, -1.);

    int nfaces = YB_f.MyLength();
    for (int f=0; f!=nfaces; ++f) {
      YB_f[0][f] /= DffB_f[0][f];
    }
  }

  return ierr;
}


void MatrixMFD_Coupled_TPFA::AssembleSchur_() const {
  int ierr(0);

  const Epetra_BlockMap& cmap = mesh_->cell_map(false);
  const Epetra_BlockMap& fmap = mesh_->face_map(false);
  const Epetra_BlockMap& fmap_wghost = mesh_->face_map(true);

  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  ASSERT(Ccc_->MyLength() == ncells);
  ASSERT(Dcc_->MyLength() == ncells);

  // Get the assorted sub-blocks
  const Epetra_FECrsMatrix& App = *blockA_TPFA_->TPFA();
  const Epetra_FECrsMatrix& Bpp = *blockB_TPFA_->TPFA();

  // workspace
  Epetra_SerialDenseMatrix values(2, 2);
  AmanziMesh::Entity_ID_List faces;

  double *valuesA;
  valuesA = new double[9];
  double *valuesB;
  valuesB = new double[9];
  int *indicesA;
  indicesA = new int[9];
  int *indicesB;
  indicesB = new int[9];

  int nentriesA = 0;
  int nentriesB = 0;

  // Assemble
  for (int c=0; c!=ncells; ++c){
    int cell_GID = cmap.GID(c);

    // Access the row in App,Bpp
    ierr = App.ExtractGlobalRowCopy(cell_GID, 9, nentriesA, valuesA, indicesA);
    ierr |= Bpp.ExtractGlobalRowCopy(cell_GID, 9, nentriesB, valuesB, indicesB);
    ASSERT(!ierr);
    ASSERT(nentriesA == nentriesB);

    ierr = P2c2c_->BeginReplaceGlobalValues(cell_GID, nentriesA, indicesA);
    ASSERT(!ierr);

    for (int n=0; n!=nentriesA; ++n) {
      ASSERT(indicesA[n] == indicesB[n]);
      values(0,0) = valuesA[n];
      values(1,1) = valuesB[n];

      if (indicesA[n] == cell_GID) {
	values(0,1) = (*Ccc_)[0][c] * scaling_;
	values(1,0) = (*Dcc_)[0][c] * scaling_;
      } else {
	values(0,1) = 0.;
	values(1,0) = 0.;
      }

      ierr |= P2c2c_->SubmitBlockEntry(values.A(), values.LDA(), values.M(), values.N());
      ASSERT(!ierr);
    }

    ierr = P2c2c_->EndSubmitEntries();
    ASSERT(!ierr);
  }

  ierr |= P2c2c_->GlobalAssemble();
  ASSERT(!ierr);

  P2f2f_ = P2c2c_; // alias for use by PC

  // DEBUG dump
  if (dump_schur_) {
    std::stringstream filename_s;
    filename_s << "schur_MatrixMFD_Coupled_TPFA_" << 0 << ".txt";
    EpetraExt::RowMatrixToMatlabFile(filename_s.str().c_str(), *P2f2f_);
  }

}


void MatrixMFD_Coupled_TPFA::SymbolicAssembleGlobalMatrices() {

  // get the standard matrices
  MatrixMFD_Coupled::SymbolicAssembleGlobalMatrices();

  // also create a cell-cell matrix for the TPFA
  const Epetra_Map& cmap = mesh_->cell_map(false);
  const Epetra_Map& cmap_wghost = mesh_->cell_map(true);
  const Epetra_Map& fmap_wghost = mesh_->face_map(true);

  int avg_entries_row = (mesh_->space_dimension() == 2) ? MFD_QUAD_FACES : MFD_HEX_FACES;
  Epetra_FECrsGraph double_pp_graph(Copy, *double_cmap_, 2*(avg_entries_row + 1), false);

  AmanziMesh::Entity_ID_List cells;
  int cells_GID[2];

  int ierr = 0;
  int nfaces = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  for (int f = 0; f < nfaces; f++) {
    mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
    int ncells = cells.size();

    for (int n = 0; n < ncells; n++) cells_GID[n] = cmap_wghost.GID(cells[n]);

    ierr |= double_pp_graph.InsertGlobalIndices(ncells, cells_GID, ncells, cells_GID);
    ASSERT(!ierr);
  }
  ierr |= double_pp_graph.GlobalAssemble();  // Symbolic graph is complete.
  ASSERT(!ierr);

  P2c2c_ = Teuchos::rcp(new Epetra_FEVbrMatrix(Copy, double_pp_graph, false));
  ierr |= P2c2c_->GlobalAssemble();
  ASSERT(!ierr);
}



} // namespace
} // namespace
