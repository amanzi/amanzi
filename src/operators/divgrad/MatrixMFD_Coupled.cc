/*
  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)

  MatrixMFD_Coupled takes two MatrixMFD objects, along with the cell coupling
  terms, and forms a coupled system that is 2x cell + 2x face sqare.

  MatrixMFD provides for a system:

      [  Acc   Acf  ]
      [  Afc   Aff  ]

  where Acc is size ncells x ncells and Aff is size nfaces x nfaces.  The key
  to solving this system is that Acc is diagonal, allowing us to precondition
  with the Schur complement:

      [  Acc   Acf  ]
      [   0    Sff  ],   Sff = Aff - Afc * Acc^-1 * Acf

  For many of the systems needed in ATS, we have two, coupled equations whose
  preconditioner looks something like:

      [  Acc   Acf    Ccc    0  ]  ( y_Ac )    ( x_Ac )
      [  Afc   Aff     0     0  ]  ( y_Af ) =  ( x_Af )
      [  Dcc    0     Bcc   Bcf ]  ( y_Bc )    ( x_Bc )
      [   0     0     Bfc   Bff ]  ( y_Bf )    ( x_Bf )

  where Acc, Bcc, Ccc, and Dcc are all diagonal.  Effectively these four
  diagonal matrices form a block-diagonal matrix where each cell corresponds
  to a 2x2 dense matrix on the diagonal.  The same Schur complement trick can
  then be applied:

     S_(2f x 2f) = [ Aff   0  ]  _  [ Afc  0  ][ Acc  Ccc ]^-1 [ Acf  0  ]
                   [  0   Bff ]     [  0  Bfc ][ Dcc  Bcc ]    [  0  Bcf ]

  This class forms this matrix.  Note that the rhs is then:

     rhs_S = ( x_Af )  _  [ Afc  0  ][ Acc  Ccc ]^-1  ( x_Ac )
             ( x_Bf )     [  0  Bfc ][ Dcc  Bcc ]     ( x_Bc )

 */

#include "Epetra_FECrsGraph.h"
#include "EpetraExt_RowMatrixOut.h"

#include "errors.hh"
#include "MatrixMFD.hh"

#include "MatrixMFD_Coupled.hh"

namespace Amanzi {
namespace Operators {

MatrixMFD_Coupled::MatrixMFD_Coupled(Teuchos::ParameterList& plist,
        const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) :
    plist_(plist),
    mesh_(mesh),
    assembled_schur_(false),
    is_schur_created_(false),
    is_operator_created_(false) {
  InitializeFromPList_();
}


MatrixMFD_Coupled::MatrixMFD_Coupled(const MatrixMFD_Coupled& other) :
    plist_(other.plist_),
    mesh_(other.mesh_),
    assembled_schur_(false),
    is_schur_created_(false) {
  InitializeFromPList_();
}


void MatrixMFD_Coupled::InitializeFromPList_() {
  // preconditioner
  if (plist_.isSublist("preconditioner")) {
    Teuchos::ParameterList pc_list = plist_.sublist("preconditioner");
    AmanziPreconditioners::PreconditionerFactory pc_fac;
    S_pc_ = pc_fac.Create(pc_list);
  }

  // verbose object
  vo_ = Teuchos::rcp(new VerboseObject("MatrixMFD", plist_));

  // dump
  dump_schur_ = plist_.get<bool>("dump Schur complement", false);
}


void MatrixMFD_Coupled::SetSubBlocks(const Teuchos::RCP<MatrixMFD>& blockA,
        const Teuchos::RCP<MatrixMFD>& blockB) {

  blockA_ = blockA;
  blockB_ = blockB;

  // set up the space
  if (space_ == Teuchos::null) {
    Teuchos::RCP<const CompositeVectorSpace> spaceA = Teuchos::rcpFromRef(blockA->DomainMap());
    Teuchos::RCP<const CompositeVectorSpace> spaceB = Teuchos::rcpFromRef(blockB->DomainMap());
    Teuchos::RCP<const TreeVectorSpace> spaceA_TV = Teuchos::rcp(new TreeVectorSpace(spaceA));
    Teuchos::RCP<const TreeVectorSpace> spaceB_TV = Teuchos::rcp(new TreeVectorSpace(spaceB));
    space_ = Teuchos::rcp(new TreeVectorSpace());
    space_->PushBack(spaceA_TV);
    space_->PushBack(spaceB_TV);    
  }

  MarkLocalMatricesAsChanged_();
}

int MatrixMFD_Coupled::ApplyInverse(const TreeVector& X,
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
  const Epetra_MultiVector& XA_f = *XA->ViewComponent("face", false);
  const Epetra_MultiVector& XB_f = *XB->ViewComponent("face", false);

  // Temporary cell and face vectors.
  Epetra_MultiVector Xf(*double_fmap_, 1);
  Epetra_MultiVector Yf(*double_fmap_, 1);

  Teuchos::RCP<CompositeVector> A = 
    Teuchos::rcp(new CompositeVector(*XA, INIT_MODE_ZERO));
  Teuchos::RCP<CompositeVector> B = 
    Teuchos::rcp(new CompositeVector(*XB, INIT_MODE_ZERO));
  Teuchos::RCP<const CompositeVector> A_const = A;
  Teuchos::RCP<const CompositeVector> B_const = B;

  int ierr(0);

  int ncells = XA_c.MyLength();
  int nfaces = XA_f.MyLength();
  double val = 0.;

  // Forward Eliminate
  // Wc <-- A2c2c_Inv * (x_Ac, x_Bc)^T
  {
    Epetra_MultiVector& Ac = *A->ViewComponent("cell",false);
    Epetra_MultiVector& Bc = *B->ViewComponent("cell",false);
    for (int c=0; c!=ncells; ++c){
      Ac[0][c] = A2c2c_cells_Inv_[c](0,0)*XA_c[0][c] 
	+ A2c2c_cells_Inv_[c](0,1)*XB_c[0][c];
      Bc[0][c] = A2c2c_cells_Inv_[c](1,0)*XA_c[0][c] 
	+ A2c2c_cells_Inv_[c](1,1)*XB_c[0][c];
    }
  }

  // Wf <-- Afc * Wc
  ierr |= blockA_->ApplyAfc(*A,*A,0.);
  ierr |= blockB_->ApplyAfc(*B,*B,0.);
  AMANZI_ASSERT(!ierr);

  // Xf <-- (x_Af, x_Ac)^T - Afc * A2c2c_Inv * (x_Ac, x_Bc)^T, the rhs.
  {
    const Epetra_MultiVector& Af = *A_const->ViewComponent("face",false);
    const Epetra_MultiVector& Bf = *B_const->ViewComponent("face",false);

    for (int f=0; f!=nfaces; ++f) {
      Xf[0][2*f] = XA_f[0][f] - Af[0][f]; 
      Xf[0][2*f+1] = XB_f[0][f] - Bf[0][f]; 
    }
  }

  // Apply Schur Inverse,  Yf = Schur^-1 * Xf
  ierr = S_pc_->ApplyInverse(Xf, Yf);
  AMANZI_ASSERT(!ierr);

  // copy back into subblock
  Teuchos::RCP<CompositeVector> YA = Y.SubVector(0)->Data();
  Teuchos::RCP<CompositeVector> YB = Y.SubVector(1)->Data();
  {
    Epetra_MultiVector& YA_f = *YA->ViewComponent("face", false);
    Epetra_MultiVector& YB_f = *YB->ViewComponent("face", false);

    for (int f=0; f!=nfaces; ++f) {
      YA_f[0][f] = Yf[0][2*f];
      YB_f[0][f] = Yf[0][2*f+1];
    }
  }

  // Backward Substitution, Yc = inv( A2c2c) [  (x_Ac,x_Bc)^T - A2c2f * Yf ]
  // Yc <-- A2c2f * Yf
  ierr |= blockA_->ApplyAcf(*YA,*YA,0.);
  ierr |= adv_block_->ApplyAcf(*YA,*YB,0.);
  ierr |= blockB_->ApplyAcf(*YB,*YB,1.);
  //  ierr |= blockB_->ApplyAcf(*YB,*YB,0.);
  AMANZI_ASSERT(!ierr);

  // Yc -= (x_Ac,x_Bc)^T
  {
    Epetra_MultiVector& YA_c = *YA->ViewComponent("cell", false);
    Epetra_MultiVector& YB_c = *YB->ViewComponent("cell", false);
    for (int c=0; c!=ncells; ++c){
      double tmpA = XA_c[0][c] - YA_c[0][c];
      double tmpB = XB_c[0][c] - YB_c[0][c];

      YA_c[0][c] = A2c2c_cells_Inv_[c](0,0)*tmpA + A2c2c_cells_Inv_[c](0,1)*tmpB;
      YB_c[0][c] = A2c2c_cells_Inv_[c](1,0)*tmpA + A2c2c_cells_Inv_[c](1,1)*tmpB;
    }
  }

  return ierr;
}


/* ******************************************************************
 * Parallel matvec product Y <-- A * X.
 ****************************************************************** */
int MatrixMFD_Coupled::Apply(const TreeVector& X,
        TreeVector& Y) const {
  AMANZI_ASSERT(0); // this is currently screwed up with the inclusion of advective terms

  Teuchos::RCP<const CompositeVector> XA = X.SubVector(0)->Data();
  Teuchos::RCP<const CompositeVector> XB = X.SubVector(1)->Data();
  Teuchos::RCP<CompositeVector> YA = Y.SubVector(0)->Data();
  Teuchos::RCP<CompositeVector> YB = Y.SubVector(1)->Data();

  int ierr = blockA_->Apply(*XA, *YA);
  ierr |= blockB_->Apply(*XB, *YB);

  // add in the off-diagonals
  ierr |= YA->ViewComponent("cell",false)->Multiply(scaling_, *Ccc_,
          *XB->ViewComponent("cell",false), 1.);
  ierr |= YB->ViewComponent("cell",false)->Multiply(scaling_, *Dcc_,
          *XA->ViewComponent("cell",false), 1.);
  AMANZI_ASSERT(!ierr);
  return ierr;
}


void MatrixMFD_Coupled::AssembleAff_() const {
  int ierr(0);

  const Epetra_BlockMap& cmap = mesh_->cell_map(false);
  const Epetra_BlockMap& fmap = mesh_->face_map(false);
  const Epetra_BlockMap& fmap_wghost = mesh_->face_map(true);

  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

  // Get the assorted sub-blocks
  std::vector<Teuchos::SerialDenseMatrix<int, double> >& Aff = blockA_->Aff_cells();
  std::vector<Teuchos::SerialDenseMatrix<int, double> >& Bff = blockB_->Aff_cells();

  // workspace
  Epetra_SerialDenseMatrix values(2, 2);
  AmanziMesh::Entity_ID_List faces;
  const int MFD_MAX_FACES = 14;
  int faces_LID[MFD_MAX_FACES];  // Contigious memory is required.
  int faces_GID[MFD_MAX_FACES];

  if (is_operator_created_) A2f2f_->PutScalar(0.0);

  // Assemble
  for (int c=0; c!=ncells; ++c){
    int cell_GID = cmap.GID(c);
    mesh_->cell_get_faces(c, &faces);
    int nfaces = faces.size();
    int nentries = nfaces; // not sure if this is required, but may be passed by ref

    Epetra_SerialDenseMatrix S2f2f(2*nfaces, 2*nfaces);

    // get IDs of faces
    for (int i=0; i!=nfaces; ++i) {
      faces_LID[i] = faces[i];
      faces_GID[i] = fmap_wghost.GID(faces_LID[i]);
    }

    for (int i=0; i!=nfaces; ++i) {
      for (int j=0; j!=nfaces; ++j) {
        S2f2f(i, j) = Aff[c](i, j);
        S2f2f(nfaces + i, nfaces + j) = Bff[c](i, j);
      }
    }

    // -- Assemble Schur complement
    for (int i=0; i!=nfaces; ++i) {
      ierr = A2f2f_->BeginSumIntoGlobalValues(faces_GID[i], nentries, faces_GID);
      AMANZI_ASSERT(!ierr);

      for (int j=0; j!=nfaces; ++j){
        values(0,0) = S2f2f(i,j);
        values(0,1) = S2f2f(i,j + nfaces);
        values(1,0) = S2f2f(i + nfaces,j);
        values(1,1) = S2f2f(i+ nfaces,j+ nfaces);

        //ierr = A2f2f_->SubmitBlockEntry(values);  // Bug in Trilinos 10.10 FeVbrMatrix
        ierr = A2f2f_->SubmitBlockEntry(values.A(), values.LDA(),
                values.M(), values.N());
        AMANZI_ASSERT(!ierr);
      }

      ierr = A2f2f_->EndSubmitEntries();
      AMANZI_ASSERT(!ierr);
    }
  }

  // Finish assembly
  ierr = A2f2f_->GlobalAssemble();
  is_operator_created_ = true;
  assembled_operator_ = true;
  AMANZI_ASSERT(!ierr);
}


/* ******************************************************************
 * Assemble Schur complement from elemental matrices.
 ****************************************************************** */
void MatrixMFD_Coupled::AssembleSchur_() const {
  int ierr(0);

  const Epetra_BlockMap& cmap = mesh_->cell_map(false);
  const Epetra_BlockMap& fmap = mesh_->face_map(false);
  const Epetra_BlockMap& fmap_wghost = mesh_->face_map(true);

  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  AMANZI_ASSERT(Ccc_->MyLength() == ncells);
  AMANZI_ASSERT(Dcc_->MyLength() == ncells);

  // Get the assorted sub-blocks
  std::vector<Teuchos::SerialDenseMatrix<int, double> >& Aff = blockA_->Aff_cells();
  std::vector<double>& Acc = blockA_->Acc_cells();
  std::vector<Epetra_SerialDenseVector>& Afc = blockA_->Afc_cells();
  std::vector<Epetra_SerialDenseVector>& Acf = blockA_->Acf_cells();

  std::vector<Teuchos::SerialDenseMatrix<int, double> >& Bff = blockB_->Aff_cells();
  std::vector<double>& Bcc = blockB_->Acc_cells();
  std::vector<Epetra_SerialDenseVector>& Bfc = blockB_->Afc_cells();
  std::vector<Epetra_SerialDenseVector>& Bcf = blockB_->Acf_cells();

  std::vector<double>& Gcc = adv_block_->Acc_cells();
  std::vector<Epetra_SerialDenseVector>& Gcf = adv_block_->Acf_cells();
  
  // workspace
  Teuchos::SerialDenseMatrix<int, double> cell_inv(2, 2);
  Epetra_SerialDenseMatrix values(2, 2);
  AmanziMesh::Entity_ID_List faces;
  const int MFD_MAX_FACES = 14;
  int faces_LID[MFD_MAX_FACES];  // Contigious memory is required.
  int faces_GID[MFD_MAX_FACES];

  // allocate global space, if necessary
  if (A2c2c_cells_Inv_.size() != ncells) {
    A2c2c_cells_Inv_.resize(static_cast<size_t>(ncells));
  }

  if (is_schur_created_) P2f2f_->PutScalar(0.0);

  // Assemble
  for (int c=0; c!=ncells; ++c){
    int cell_GID = cmap.GID(c);
    mesh_->cell_get_faces(c, &faces);
    int nfaces = faces.size();
    int nentries = nfaces; // not sure if this is required, but may be passed by ref

    // space for the local matrices
    Epetra_SerialDenseMatrix S2f2f(2*nfaces, 2*nfaces);

    // get IDs of faces
    for (int i=0; i!=nfaces; ++i) {
      faces_LID[i] = faces[i];
      faces_GID[i] = fmap_wghost.GID(faces_LID[i]);
    }

    // Invert the cell block
    double Dcc = (*Dcc_)[0][c]*scaling_ + Gcc[c];
    double Ccc = (*Ccc_)[0][c]*scaling_;
    double det_cell = Acc[c] * Bcc[c] - Ccc * Dcc;
    if (std::abs(det_cell) > 1.e-30) {
      cell_inv(0, 0) = Bcc[c] / det_cell;
      cell_inv(1, 1) = Acc[c] / det_cell;
      cell_inv(0, 1) = -Ccc / det_cell;
      cell_inv(1, 0) = -Dcc / det_cell;
    } else {
      std::cout << "MatrixMFD_Coupled: Division by zero: determinant of the cell block is zero" << std::endl;
      //      AMANZI_ASSERT(0);
      Exceptions::amanzi_throw(Errors::CutTimeStep());
    }
    A2c2c_cells_Inv_[c] = cell_inv;

    // -- Assemble Schur complement
    for (int i=0; i!=nfaces; ++i) {
      ierr = P2f2f_->BeginSumIntoGlobalValues(faces_GID[i], nentries, faces_GID);
      AMANZI_ASSERT(!ierr);

      for (int j=0; j!=nfaces; ++j){
        values(0,0) = Aff[c](i, j) - Afc[c](i) * (cell_inv(0,0)*Acf[c](j) + cell_inv(0,1)*Gcf[c](j));
        values(0,1) = - Afc[c](i)*cell_inv(0,1)*Bcf[c](j);
        values(1,0) = - Bfc[c](i) * (cell_inv(1,0)*Acf[c](j) + cell_inv(1,1)*Gcf[c](j));
        values(1,1) = Bff[c](i, j) - Bfc[c](i)*cell_inv(1,1)*Bcf[c](j);
        ierr = P2f2f_->SubmitBlockEntry(values.A(), values.LDA(),
                values.M(), values.N());
        AMANZI_ASSERT(!ierr);
      }

      ierr = P2f2f_->EndSubmitEntries();
      AMANZI_ASSERT(!ierr);
    }
  }

  // Finish assembly
  ierr = P2f2f_->GlobalAssemble();
  AMANZI_ASSERT(!ierr);
  assembled_schur_ = true;
  is_schur_created_ = true;

  // DEBUG dump
  if (dump_schur_) {
    std::stringstream filename_s;
    filename_s << "schur_MatrixMFD_Coupled_" << 0 << ".txt";
    EpetraExt::RowMatrixToMatlabFile(filename_s.str().c_str(), *P2f2f_);
  }

}


/* ******************************************************************
 * Initialize global matrices.
 *
 * This likely should only be called once.
 ****************************************************************** */
void MatrixMFD_Coupled::SymbolicAssembleGlobalMatrices() {
  int ierr(0);
  const Epetra_BlockMap& cmap = mesh_->cell_map(false);
  const Epetra_BlockMap& fmap = mesh_->face_map(false);
  const Epetra_BlockMap& fmap_wghost = mesh_->face_map(true);

  // Make the double maps
  double_fmap_ = Teuchos::rcp(new Epetra_BlockMap(fmap.NumGlobalPoints(),
          fmap.NumMyPoints(), fmap.MyGlobalElements(), 2,
          fmap.IndexBase(), fmap.Comm() ));

  double_cmap_ = Teuchos::rcp(new Epetra_BlockMap(cmap.NumGlobalPoints(),
          cmap.NumMyPoints(), cmap.MyGlobalElements(), 2,
          cmap.IndexBase(), cmap.Comm() ));

  double_fmap_wghost_ = Teuchos::rcp(
      new Epetra_BlockMap(fmap_wghost.NumGlobalPoints(),
                          fmap_wghost.NumMyPoints(),
                          fmap_wghost.MyGlobalElements(), 2,
                          fmap_wghost.IndexBase(), fmap_wghost.Comm() ));

  // Make the matrix graphs
  int avg_entries_row = 6;
  Teuchos::RCP<Epetra_CrsGraph> cf_graph =
      Teuchos::rcp(new Epetra_CrsGraph(Copy, *double_cmap_, *double_fmap_wghost_,
              avg_entries_row, false));
  Teuchos::RCP<Epetra_FECrsGraph> ff_graph =
      Teuchos::rcp(new Epetra_FECrsGraph(Copy, *double_fmap_, 2*avg_entries_row - 1,
              false));

  // Fill the graphs from the symbolic patterns in matrix_mfd from one of the
  // sub-block matrices.

  // IT IS ASSUMED that the two sub-block matrices fill in identical ways!
  blockA_->FillMatrixGraphs_(cf_graph.ptr(), ff_graph.ptr());

  // Assemble the graphs
  ierr = cf_graph->FillComplete(*double_fmap_, *double_cmap_);
  AMANZI_ASSERT(!ierr);
  ierr = ff_graph->GlobalAssemble();  // Symbolic graph is complete.
  AMANZI_ASSERT(!ierr);

  // Create the matrices
  A2f2c_ = Teuchos::rcp(new Epetra_VbrMatrix(Copy, *cf_graph)); // stored in transpose
  A2c2f_ = Teuchos::rcp(new Epetra_VbrMatrix(Copy, *cf_graph));
  P2f2f_ = Teuchos::rcp(new Epetra_FEVbrMatrix(Copy, *ff_graph, false));
  //  A2f2f_ = Teuchos::rcp(new Epetra_FEVbrMatrix(Copy, *ff_graph, false));
  ierr = P2f2f_->GlobalAssemble();
  //  ierr = A2f2f_->GlobalAssemble();
  AMANZI_ASSERT(!ierr);
}



/* ******************************************************************
 * Rebuild preconditioner.
 ****************************************************************** */
void MatrixMFD_Coupled::UpdatePreconditioner_() const {
  if (S_pc_ == Teuchos::null) {
    Errors::Message msg("MatrixMFD::UpdatePreconditioner called but no preconditioner sublist was provided");
    Exceptions::amanzi_throw(msg);
  }
  S_pc_->Destroy();
  S_pc_->Update(P2f2f_);
}


/* ******************************************************************
 * Solve the bottom row of the block system for lambda, given p.
 ****************************************************************** */
void MatrixMFD_Coupled::UpdateConsistentFaceCorrection(const TreeVector& u,
        const Teuchos::Ptr<TreeVector>& Pu) {
  blockA_->UpdateConsistentFaceCorrection(*u.SubVector(0)->Data(),
          Pu->SubVector(0)->Data().ptr());
  blockB_->UpdateConsistentFaceCorrection(*u.SubVector(1)->Data(),
          Pu->SubVector(1)->Data().ptr());
}

} // namespace
} // namespace
