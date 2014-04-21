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
    is_matrix_constructed_(false) {
  InitializeFromPList_();
}


MatrixMFD_Coupled::MatrixMFD_Coupled(const MatrixMFD_Coupled& other) :
    plist_(other.plist_),
    mesh_(other.mesh_),
    is_matrix_constructed_(false) {
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
  Teuchos::RCP<const CompositeVectorSpace> spaceA = Teuchos::rcpFromRef(blockA->DomainMap());
  Teuchos::RCP<const CompositeVectorSpace> spaceB = Teuchos::rcpFromRef(blockB->DomainMap());
  Teuchos::RCP<const TreeVectorSpace> spaceA_TV = Teuchos::rcp(new TreeVectorSpace(spaceA));
  Teuchos::RCP<const TreeVectorSpace> spaceB_TV = Teuchos::rcp(new TreeVectorSpace(spaceB));
  space_ = Teuchos::rcp(new TreeVectorSpace());
  space_->PushBack(spaceA_TV);
  space_->PushBack(spaceB_TV);    
}

int MatrixMFD_Coupled::ApplyInverse(const TreeVector& X,
        TreeVector& Y) const {
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
  Epetra_MultiVector Xc(*double_cmap_, 1);
  Epetra_MultiVector Xf(*double_fmap_, 1);

  Epetra_MultiVector Yc(*double_cmap_, 1);
  Epetra_MultiVector Yf(*double_fmap_, 1);

  int ierr(0);

  int ncells = XA_c.MyLength();
  int nfaces = XA_f.MyLength();
  double val = 0.;

  // Forward Eliminate
  // Xc <-- - A2c2c_Inv * (x_Ac, x_Bc)^T
  for (int c=0; c!=ncells; ++c){
    double a_c = XA_c[0][c];
    double b_c = XB_c[0][c];
    val = -(A2c2c_cells_Inv_[c](0,0)*a_c + A2c2c_cells_Inv_[c](0,1)*b_c);
    ierr = Xc.ReplaceMyValue(c, 0, 0, val);
    ASSERT(!ierr);

    val = -(A2c2c_cells_Inv_[c](1,0)*a_c + A2c2c_cells_Inv_[c](1,1)*b_c);
    ierr = Xc.ReplaceMyValue(c, 1, 0, val);
    ASSERT(!ierr);
  }

  // Xf <-- A2f2c * Xc
  ierr = A2f2c_->Multiply(true, Xc, Xf);
  ASSERT(!ierr);

  // Xf <-- (x_Af, x_Ac)^T + Xf
  // This gives Xf = (x_Af, x_Ac)^T  - A2f2c * A2c2c_Inv * (x_Ac, x_Bc)^T, the rhs.
  for (int f=0; f!=nfaces; ++f){
    ierr = Xf.SumIntoMyValue(f, 0, 0, XA_f[0][f]);
    ASSERT(!ierr);
    ierr = Xf.SumIntoMyValue(f, 1, 0, XB_f[0][f]);
    ASSERT(!ierr);
  }

  // Apply Schur Inverse,  Yf = Schur^-1 * Xf
  ierr = S_pc_->ApplyInverse(Xf, Yf);
  ASSERT(!ierr);

  // Backward Substitution, Yc = inv( A2c2c) [  (x_Ac,x_Bc)^T - A2c2f * Yf ]
  // Yc <-- A2c2f * Yf
  ierr = A2c2f_->Multiply(false, Yf, Yc);
  ASSERT(!ierr);


  // Yc -= (x_Ac,x_Bc)^T
  for (int c=0; c!=ncells; ++c){
    double a_c = -XA_c[0][c];
    double b_c = -XB_c[0][c];
    Yc.SumIntoMyValue(c, 0, 0, a_c);
    Yc.SumIntoMyValue(c, 1, 0, b_c);
  }

  // Yc <-- -Yc, now Yc =  (x_Ac,x_Bc)^T - A2c2f * Yf
  Yc.Scale(-1.0);

  // pull Y data
  Teuchos::RCP<CompositeVector> YA = Y.SubVector(0)->Data();
  Teuchos::RCP<CompositeVector> YB = Y.SubVector(1)->Data();
  Epetra_MultiVector& YA_c = *YA->ViewComponent("cell", false);
  Epetra_MultiVector& YB_c = *YB->ViewComponent("cell", false);
  Epetra_MultiVector& YA_f = *YA->ViewComponent("face", false);
  Epetra_MultiVector& YB_f = *YB->ViewComponent("face", false);

  // Yc <-- inv(A2c2c) * Yc
  // Copy data back into Y-vectors
  for (int c=0; c!=ncells; ++c){
    YA_c[0][c] = A2c2c_cells_Inv_[c](0,0)*Yc[0][2*c]
        + A2c2c_cells_Inv_[c](0,1)*Yc[0][2*c + 1];
    YB_c[0][c] = A2c2c_cells_Inv_[c](1,0)*Yc[0][2*c]
        + A2c2c_cells_Inv_[c](1,1)*Yc[0][2*c + 1];
  }

  for (int f=0; f!=nfaces; ++f){
    YA_f[0][f] = Yf[0][2*f];
    YB_f[0][f] = Yf[0][2*f + 1];
  }

  return ierr;
}


int MatrixMFD_Coupled::Apply(const TreeVector& X,
        TreeVector& Y) const {
  Teuchos::RCP<const CompositeVector> XA = X.SubVector(0)->Data();
  Teuchos::RCP<const CompositeVector> XB = X.SubVector(1)->Data();
  Teuchos::RCP<CompositeVector> YA = Y.SubVector(0)->Data();
  Teuchos::RCP<CompositeVector> YB = Y.SubVector(1)->Data();

  blockA_->Apply(*XA, *YA);
  blockB_->Apply(*XB, *YB);

  // add in the off-diagonals
  YA->ViewComponent("cell",false)->Multiply(scaling_, *Ccc_,
          *XB->ViewComponent("cell",false), 1.);
  YB->ViewComponent("cell",false)->Multiply(scaling_, *Dcc_,
          *XA->ViewComponent("cell",false), 1.);
  return 0;
}


void MatrixMFD_Coupled::ComputeSchurComplement(bool dump) {
  int ierr(0);

  const Epetra_BlockMap& cmap = mesh_->cell_map(false);
  const Epetra_BlockMap& fmap = mesh_->face_map(false);
  const Epetra_BlockMap& fmap_wghost = mesh_->face_map(true);

  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  ASSERT(Ccc_->MyLength() == ncells);
  ASSERT(Dcc_->MyLength() == ncells);

  // Get the assorted sub-blocks
  std::vector<Teuchos::SerialDenseMatrix<int, double> >& Aff = blockA_->Aff_cells();
  std::vector<double>& Acc = blockA_->Acc_cells();
  std::vector<Epetra_SerialDenseVector>& Afc = blockA_->Afc_cells();
  std::vector<Epetra_SerialDenseVector>& Acf = blockA_->Acf_cells();

  std::vector<Teuchos::SerialDenseMatrix<int, double> >& Bff = blockB_->Aff_cells();
  std::vector<double>& Bcc = blockB_->Acc_cells();
  std::vector<Epetra_SerialDenseVector>& Bfc = blockB_->Afc_cells();
  std::vector<Epetra_SerialDenseVector>& Bcf = blockB_->Acf_cells();

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

  if (is_matrix_constructed_) P2f2f_->PutScalar(0.0);

  // Assemble
  for (int c=0; c!=ncells; ++c){
    int cell_GID = cmap.GID(c);
    mesh_->cell_get_faces(c, &faces);
    int nfaces = faces.size();
    int nentries = nfaces; // not sure if this is required, but may be passed by ref

    // space for the local matrices
    Epetra_SerialDenseMatrix S2f2f(2*nfaces, 2*nfaces);
    Epetra_SerialDenseMatrix A2c2f(2, 2*nfaces);
    Epetra_SerialDenseMatrix A2f2c(2, 2*nfaces);

    // get IDs of faces
    for (int i=0; i!=nfaces; ++i) {
      faces_LID[i] = faces[i];
      faces_GID[i] = fmap_wghost.GID(faces_LID[i]);
    }

    // Invert the cell block
    double det_cell = Acc[c] * Bcc[c] - (*Ccc_)[0][c] * (*Dcc_)[0][c] * scaling_ * scaling_;
    if (std::abs(det_cell) > 1.e-30) {
      cell_inv(0, 0) = Bcc[c]/det_cell;
      cell_inv(1, 1) = Acc[c]/det_cell;
      cell_inv(0, 1) = -(*Ccc_)[0][c] * scaling_ / det_cell;
      cell_inv(1, 0) = -(*Dcc_)[0][c] * scaling_ / det_cell;
    } else {
      std::cout << "MatrixMFD_Coupled: Division by zero: determinant of the cell block is zero" << std::endl;
      //      ASSERT(0);
      Exceptions::amanzi_throw(Errors::CutTimeStep());
    }
    A2c2c_cells_Inv_[c] = cell_inv;

    /*
    if (c == 131) {
      std::cout << "MatrixMFD_Coupled c(131), f(595) terms:\n"
		<< "   scaling=" << scaling_ << "\n"
		<< "  Acc(131,131) = " << Acc[c] << "\n"
		<< "  Bcc(131,131) = " << Bcc[c] << "\n"
		<< "  Ccc(131,131) = " << (*Ccc_)[0][c] * scaling_ << "\n"
		<< "  Dcc(131,131) = " << (*Dcc_)[0][c] * scaling_ << "\n"
		<< "     (det = " << det_cell << ")\n"
		<< "  ------" << std::endl;
    }
    */


    // Make the cell-local Schur complement
    for (int i=0; i!=nfaces; ++i) {
      for (int j=0; j!=nfaces; ++j) {
        S2f2f(i, j) = Aff[c](i, j) - Afc[c](i)*cell_inv(0, 0)*Acf[c](j);
        if ((i == j) && std::abs(S2f2f(i,j)) < 1.e-40) {
	  std::cout << "MatrixMFD_Coupled: Schur complement pressure diagonal is zero" << std::endl;
          //          ASSERT(0);
          //          Exceptions::amanzi_throw(Errors::CutTimeStep());
        }

	//	if (c == 131 && faces_GID[i] == 595 && i == j)
	  //	  std::cout << "  ( Aff[595]= " << Aff[c](i,j) << " ) - ( Afc*c_inv*Acf = " << Afc[c](i)*cell_inv(0, 0)*Acf[c](j) << ")" << std::endl;

      }
    }

    for (int i=0; i!=nfaces; ++i) {
      for (int j=0; j!=nfaces; ++j) {
        S2f2f(nfaces + i, nfaces + j) = Bff[c](i, j) - Bfc[c](i)*cell_inv(1, 1)*Bcf[c](j);
        if ((i == j) && std::abs(S2f2f(nfaces+i,nfaces+j)) < 1.e-40) {
          std::cout << "MatrixMFD_Coupled: Schur complement temperature diagonal is zero" << std::endl;
          //          ASSERT(0);
          //          Exceptions::amanzi_throw(Errors::CutTimeStep());
        }
	//	if (c == 131 && faces_GID[i] == 595 && i == j)
	//	  std::cout << "  ( Bff[595]= " << Bff[c](i,j) << " ) - ( Bfc*c_inv*Bcf = " << Bfc[c](i)*cell_inv(1, 1)*Bcf[c](j) << ")" << std::endl;

      }
    }

    for (int i=0; i!=nfaces; ++i) {
      for (int j=0; j!=nfaces; ++j) {
        S2f2f(i, nfaces + j) = - Afc[c](i)*cell_inv(0, 1)*Bcf[c](j);
        if (std::abs(S2f2f(i,nfaces+j)) > 1.e+21) {
          std::cout << "BREAKING!" << std::endl;
        }

	//	if (c == 131 && faces_GID[i] == 595 && i == j)
	//	  std::cout << "   - ( Afc*c_inv*Bcf = " << Afc[c](i)*cell_inv(0, 1)*Bcf[c](j) << std::endl;

      }
    }

    for (int i=0; i!=nfaces; ++i) {
      for (int j=0; j!=nfaces; ++j) {
        S2f2f(nfaces + i, j) = - Bfc[c](i)*cell_inv(1, 0)*Acf[c](j);

	//	if (c == 131 && faces_GID[i] == 595 && i == j)
	//	  std::cout << "   - ( Bfc*c_inv*Acf = " << Bfc[c](i)*cell_inv(1,0)*Acf[c](j) << "\n  -------" << std::endl;

      }
    }

    // Make the local A2c2f, A2f2c
    for (int i=0; i!=nfaces; ++i) {
      A2c2f(0,i) =           Acf[c](i);
      A2c2f(0,i + nfaces)  = 0;
      A2c2f(1,i) =           0;
      A2c2f(1,i + nfaces)  = Bcf[c](i);

      A2f2c(0,i) =           Afc[c](i);
      A2f2c(0,i + nfaces)  = 0;
      A2f2c(1,i) =           0;
      A2f2c(1,i + nfaces)  = Bfc[c](i);
    }

    // -- Assemble Schur complement
    for (int i=0; i!=nfaces; ++i) {
      ierr = P2f2f_->BeginSumIntoGlobalValues(faces_GID[i], nentries, faces_GID);
      ASSERT(!ierr);

      for (int j=0; j!=nfaces; ++j){
        values(0,0) = S2f2f(i,j);
        values(0,1) = S2f2f(i,j + nfaces);
        values(1,0) = S2f2f(i + nfaces,j);
        values(1,1) = S2f2f(i+ nfaces,j+ nfaces);

	/*
	if (c == 131 && faces_GID[i] == 595 && i == j) {
	  std::cout << "total values =" << "\n"
		    << "      Sff_pp[65,65] = " << values(0,0) << "\n"
		    << "      Sff_TT[65,65] = " << values(1,1) << "\n"
		    << "      Sff_pT[65,65] = " << values(0,1) << "\n"
		    << "      Sff_Tp[65,65] = " << values(1,0) << "\n  ------" << std::endl;
	}
	*/

        //ierr = P2f2f_->SubmitBlockEntry(values);  // Bug in Trilinos 10.10 FeVbrMatrix
        ierr = P2f2f_->SubmitBlockEntry(values.A(), values.LDA(),
                values.M(), values.N());
        ASSERT(!ierr);
      }

      ierr = P2f2f_->EndSubmitEntries();
      ASSERT(!ierr);
    }

    // -- Assemble A2f2c, which is stored as transpose
    ierr = A2f2c_->BeginReplaceGlobalValues(cell_GID, nentries, faces_GID);
    ASSERT(!ierr);

    for (int i=0; i!=nfaces; ++i) {
      values(0,0) = A2f2c(0,i);
      values(0,1) = A2f2c(0,i + nfaces);
      values(1,0) = A2f2c(1,i);
      values(1,1) = A2f2c(1,i + nfaces);
      ierr = A2f2c_->SubmitBlockEntry(values);
      ASSERT(!ierr);
    }

    ierr = A2f2c_->EndSubmitEntries();
    ASSERT(!ierr);

    // -- Assemble A2c2f
    ierr = A2c2f_->BeginReplaceGlobalValues(cell_GID, nentries, faces_GID);
    ASSERT(!ierr);

    for (int i=0; i!=nfaces; ++i) {
      values(0,0) = A2c2f(0,i);
      values(0,1) = A2c2f(0,i + nfaces);
      values(1,0) = A2c2f(1,i);
      values(1,1) = A2c2f(1,i + nfaces);
      ierr = A2c2f_->SubmitBlockEntry(values);
      ASSERT(!ierr);
    }

    ierr = A2c2f_->EndSubmitEntries();
    ASSERT(!ierr);

  }

  // Finish assembly
  ierr = A2f2c_->FillComplete(*double_fmap_, *double_cmap_);
  ASSERT(!ierr);
  ierr = A2c2f_->FillComplete(*double_fmap_, *double_cmap_);
  ASSERT(!ierr);
  ierr = P2f2f_->GlobalAssemble();
  ASSERT(!ierr);
  is_matrix_constructed_ = true;

  // DEBUG dump
  if (dump || dump_schur_) {
    std::stringstream filename_s;
    filename_s << "schur_MatrixMFD_Coupled_" << 0 << ".txt";
    EpetraExt::RowMatrixToMatlabFile(filename_s.str().c_str(), *P2f2f_);
  }

}


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
  ASSERT(!ierr);
  ierr = ff_graph->GlobalAssemble();  // Symbolic graph is complete.
  ASSERT(!ierr);

  // Create the matrices
  A2f2c_ = Teuchos::rcp(new Epetra_VbrMatrix(Copy, *cf_graph)); // stored in transpose
  A2c2f_ = Teuchos::rcp(new Epetra_VbrMatrix(Copy, *cf_graph));
  P2f2f_ = Teuchos::rcp(new Epetra_FEVbrMatrix(Copy, *ff_graph, false));
  ierr = P2f2f_->GlobalAssemble();
  ASSERT(!ierr);
}


/* ******************************************************************
 * Initialization of the preconditioner
 ****************************************************************** */
void MatrixMFD_Coupled::InitPreconditioner() {}

/* ******************************************************************
 * Rebuild preconditioner.
 ****************************************************************** */
void MatrixMFD_Coupled::UpdatePreconditioner() {
  if (S_pc_ == Teuchos::null) {
    Errors::Message msg("MatrixMFD::UpdatePreconditioner called but no preconditioner sublist was provided");
    Exceptions::amanzi_throw(msg);
  }
  S_pc_->Destroy();
  S_pc_->Update(P2f2f_);
}


void MatrixMFD_Coupled::UpdateConsistentFaceCorrection(const TreeVector& u,
        const Teuchos::Ptr<TreeVector>& Pu) {
  blockA_->UpdateConsistentFaceCorrection(*u.SubVector(0)->Data(),
          Pu->SubVector(0)->Data().ptr());
  blockB_->UpdateConsistentFaceCorrection(*u.SubVector(1)->Data(),
          Pu->SubVector(1)->Data().ptr());
}

} // namespace
} // namespace
