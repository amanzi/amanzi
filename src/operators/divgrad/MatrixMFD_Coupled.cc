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
        const Teuchos::RCP<const AmanziMesh::Mesh> mesh) :
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
  prec_method_ = PREC_METHOD_NULL;
  if (plist_.isParameter("preconditioner")) {
    std::string precmethodstring = plist_.get<string>("preconditioner");
    if (precmethodstring == "ML") {
      prec_method_ = TRILINOS_ML;
    } else if (precmethodstring == "ILU" ) {
      prec_method_ = TRILINOS_ILU;
    } else if (precmethodstring == "Block ILU" ) {
      prec_method_ = TRILINOS_BLOCK_ILU;
#ifdef HAVE_HYPRE
    } else if (precmethodstring == "HYPRE AMG") {
      prec_method_ = HYPRE_AMG;
    } else if (precmethodstring == "HYPRE Euclid") {
      prec_method_ = HYPRE_EUCLID;
    } else if (precmethodstring == "HYPRE ParaSails") {
      prec_method_ = HYPRE_EUCLID;
#endif
    } else {
#ifdef HAVE_HYPRE
      Errors::Message msg("Matrix_MFD: The specified preconditioner "+precmethodstring+" is not supported, we only support ML, ILU, HYPRE AMG, HYPRE Euclid, and HYPRE ParaSails");
#else
      Errors::Message msg("Matrix_MFD: The specified preconditioner "+precmethodstring+" is not supported, we only support ML, and ILU");
#endif
      Exceptions::amanzi_throw(msg);
    }
  }
}


int MatrixMFD_Coupled::ApplyInverse(const TreeVector& X,
        TreeVector& Y) const {
  if (prec_method_ == PREC_METHOD_NULL) {
    Errors::Message msg("MatrixMFD::ApplyInverse requires a specified preconditioner method");
    Exceptions::amanzi_throw(msg);
  }

  // pull X data
  Teuchos::RCP<const CompositeVector> XA = X.SubVector(0)->data();
  Teuchos::RCP<const CompositeVector> XB = X.SubVector(1)->data();
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
  if (prec_method_ == TRILINOS_ML) {
    ierr = ml_prec_->ApplyInverse(Xf, Yf);
  } else if (prec_method_ == TRILINOS_ILU) {
    ierr = ilu_prec_->ApplyInverse(Xf, Yf);
  } else if (prec_method_ == TRILINOS_BLOCK_ILU) {
    ierr = ifp_prec_->ApplyInverse(Xf, Yf);
#ifdef HAVE_HYPRE
  } else if (prec_method_ == HYPRE_AMG || prec_method_ == HYPRE_EUCLID || prec_method_ == HYPRE_PARASAILS) {
    ierr = IfpHypre_Sff_->ApplyInverse(Xf, Yf);
#endif
  } else {
    ASSERT(0);
  }
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
  Teuchos::RCP<CompositeVector> YA = Y.SubVector(0)->data();
  Teuchos::RCP<CompositeVector> YB = Y.SubVector(1)->data();
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
  Teuchos::RCP<const CompositeVector> XA = X.SubVector(0)->data();
  Teuchos::RCP<const CompositeVector> XB = X.SubVector(1)->data();
  Teuchos::RCP<CompositeVector> YA = Y.SubVector(0)->data();
  Teuchos::RCP<CompositeVector> YB = Y.SubVector(1)->data();

  blockA_->Apply(*XA, *YA);
  blockB_->Apply(*XB, *YB);

  // add in the off-diagonals
  YA->ViewComponent("cell",false)->Multiply(1., *Ccc_,
          *XB->ViewComponent("cell",false), 1.);
  YB->ViewComponent("cell",false)->Multiply(1., *Dcc_,
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
  std::vector<int> dirs;
  const int MFD_MAX_FACES = 14;
  int faces_LID[MFD_MAX_FACES];  // Contigious memory is required.
  int faces_GID[MFD_MAX_FACES];

  // clear global space
  A2c2c_cells_Inv_.clear();
  if (is_matrix_constructed_) P2f2f_->PutScalar(0.0);

  // Assemble
  for (int c=0; c!=ncells; ++c){
    int cell_GID = cmap.GID(c);
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
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
    double det_cell = Acc[c] * Bcc[c] - (*Ccc_)[0][c] * (*Dcc_)[0][c];
    if (std::abs(det_cell) > 1.e-30) {
      cell_inv(0, 0) = Bcc[c]/det_cell;
      cell_inv(1, 1) = Acc[c]/det_cell;
      cell_inv(0, 1) = -(*Ccc_)[0][c]/det_cell;
      cell_inv(1, 0) = -(*Dcc_)[0][c]/det_cell;
    } else {
      std::cout << "MatrixMFD_Coupled: Division by zero: determinant of the cell block is zero" << std::endl;
      //      ASSERT(0);
      Exceptions::amanzi_throw(Errors::CutTimeStep());
    }
    A2c2c_cells_Inv_.push_back(cell_inv);

    // Make the cell-local Schur complement
    for (int i=0; i!=nfaces; ++i) {
      for (int j=0; j!=nfaces; ++j) {
        S2f2f(i, j) = Aff[c](i, j) - Afc[c](i)*cell_inv(0, 0)*Acf[c](j);
        if ((i == j) && std::abs(S2f2f(i,j)) < 1.e-40) {
          std::cout << "MatrixMFD_Coupled: Schur complement pressure diagonal is zero" << std::endl;
          //          ASSERT(0);
          //          Exceptions::amanzi_throw(Errors::CutTimeStep());
        }

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
      }
    }

    for (int i=0; i!=nfaces; ++i) {
      for (int j=0; j!=nfaces; ++j) {
        S2f2f(i, nfaces + j) = - Afc[c](i)*cell_inv(0, 1)*Bcf[c](j);
        if (std::abs(S2f2f(i,nfaces+j)) > 1.e+21) {
          std::cout << "BREAKING!" << std::endl;
        }
      }
    }

    for (int i=0; i!=nfaces; ++i) {
      for (int j=0; j!=nfaces; ++j) {
        S2f2f(nfaces + i, j) = - Bfc[c](i)*cell_inv(1, 0)*Acf[c](j);
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
  if (dump) {
    std::stringstream filename_s;
    filename_s << "schur_" << 0 << ".txt";
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
void MatrixMFD_Coupled::InitPreconditioner() {
  if (prec_method_ == TRILINOS_ML) {
    ml_plist_ =  plist_.sublist("ML Parameters");
    ml_prec_ = Teuchos::rcp(new ML_Epetra::MultiLevelPreconditioner(*P2f2f_, ml_plist_, false));
  } else if (prec_method_ == TRILINOS_ILU) {
    ilu_plist_ = plist_.sublist("ILU Parameters");
  } else if (prec_method_ == TRILINOS_BLOCK_ILU) {
    ifp_plist_ = plist_.sublist("Block ILU Parameters");
#ifdef HAVE_HYPRE
  } else if (prec_method_ == HYPRE_AMG) {
    // read some boomer amg parameters
    hypre_plist_ = plist_.sublist("HYPRE AMG Parameters");

    hypre_ncycles_ = hypre_plist_.get<int>("number of cycles",5);
    hypre_nsmooth_ = hypre_plist_.get<int>("number of smoothing iterations",3);
    hypre_tol_ = hypre_plist_.get<double>("tolerance",0.0);
    hypre_strong_threshold_ = hypre_plist_.get<double>("strong threshold",0.25);
    hypre_cycle_type_ = hypre_plist_.get<int>("cycle type",1);
    hypre_relax_type_ = hypre_plist_.get<int>("relax type",6);
    hypre_coarsen_type_ = hypre_plist_.get<int>("coarsen type",0);
    hypre_max_row_sum_ = hypre_plist_.get<double>("max row sum",0.9);
    hypre_max_levels_ = hypre_plist_.get<int>("max levels",25);
    hypre_relax_wt_ = hypre_plist_.get<double>("relax wt",1.0);
    hypre_interp_type_ = hypre_plist_.get<int>("interpolation type",0);
    hypre_agg_num_levels_ = hypre_plist_.get<int>("aggressive coarsening levels",0);
    hypre_agg_num_paths_ = hypre_plist_.get<int>("aggressive coarsening paths",1);
    hypre_print_level_ =  hypre_plist_.get<int>("print level",0);

  } else if (prec_method_ == HYPRE_EUCLID) {
    hypre_plist_ = plist_.sublist("HYPRE Euclid Parameters");
  } else if (prec_method_ == HYPRE_PARASAILS) {
    hypre_plist_ = plist_.sublist("HYPRE ParaSails Parameters");
#endif
  }
}

/* ******************************************************************
 * Rebuild preconditioner.
 ****************************************************************** */
void MatrixMFD_Coupled::UpdatePreconditioner() {
  ASSERT(is_matrix_constructed_);

  if (prec_method_ == TRILINOS_ML) {
    if (ml_prec_->IsPreconditionerComputed()) ml_prec_->DestroyPreconditioner();
    ml_prec_->SetParameterList(ml_plist_);
    ml_prec_->ComputePreconditioner();
  } else if (prec_method_ == TRILINOS_ILU) {
    ilu_prec_ = Teuchos::rcp(new Ifpack_ILU(&*P2f2f_));
    ilu_prec_->SetParameters(ilu_plist_);
    ilu_prec_->Initialize();
    ilu_prec_->Compute();
  } else if (prec_method_ == TRILINOS_BLOCK_ILU) {
    Ifpack factory;
    std::string prectype("ILU");
    int ovl = ifp_plist_.get<int>("overlap",0);
    ifp_plist_.set<std::string>("schwarz: combine mode","Add");
    ifp_prec_ = Teuchos::rcp(factory.Create(prectype, &*P2f2f_, ovl));
    ifp_prec_->SetParameters(ifp_plist_);
    ifp_prec_->Initialize();
    ifp_prec_->Compute();
#ifdef HAVE_HYPRE
  } else if (prec_method_ == HYPRE_AMG) {
    IfpHypre_Sff_ = Teuchos::rcp(new Ifpack_Hypre(&*P2f2f_));
    Teuchos::RCP<FunctionParameter> functs[14];
    
    functs[0] = Teuchos::rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetCoarsenType, hypre_coarsen_type_));
    functs[1] = Teuchos::rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetPrintLevel, hypre_print_level_));
    functs[2] = Teuchos::rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetNumSweeps, hypre_nsmooth_));
    functs[3] = Teuchos::rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetMaxIter, hypre_ncycles_));
    functs[4] = Teuchos::rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetRelaxType, hypre_relax_type_));
    functs[5] = Teuchos::rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetStrongThreshold, hypre_strong_threshold_));
    functs[6] = Teuchos::rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetTol, hypre_tol_));
    functs[7] = Teuchos::rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetCycleType, hypre_cycle_type_));
    functs[8] = Teuchos::rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetMaxRowSum, hypre_max_row_sum_));
    functs[9] = Teuchos::rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetMaxLevels, hypre_max_levels_));
    functs[10] = Teuchos::rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetInterpType, hypre_interp_type_));
    functs[11] = Teuchos::rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetAggNumLevels, hypre_agg_num_levels_));
    functs[12] = Teuchos::rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetNumPaths, hypre_agg_num_paths_));
    functs[13] = Teuchos::rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetRelaxWt, hypre_relax_wt_));
        

    Teuchos::ParameterList hypre_list;
    hypre_list.set("Preconditioner", BoomerAMG);
    hypre_list.set("SolveOrPrecondition", Preconditioner);
    hypre_list.set("SetPreconditioner", true);
    hypre_list.set("NumFunctions", 14);
    hypre_list.set<Teuchos::RCP<FunctionParameter>*>("Functions", functs);

    IfpHypre_Sff_->SetParameters(hypre_list);
    IfpHypre_Sff_->Initialize();
    IfpHypre_Sff_->Compute();
  } else if (prec_method_ == HYPRE_EUCLID) {
    IfpHypre_Sff_ = Teuchos::rcp(new Ifpack_Hypre(&*P2f2f_));

    Teuchos::ParameterList hypre_list;
    hypre_list.set("Preconditioner", Euclid);
    hypre_list.set("SolveOrPrecondition", Preconditioner);
    hypre_list.set("SetPreconditioner", true);
    hypre_list.set("NumFunctions", 0);

    IfpHypre_Sff_->SetParameters(hypre_list);
    IfpHypre_Sff_->Initialize();
    IfpHypre_Sff_->Compute();
  } else if (prec_method_ == HYPRE_PARASAILS) {
    IfpHypre_Sff_ = Teuchos::rcp(new Ifpack_Hypre(&*P2f2f_));

    Teuchos::ParameterList hypre_list;
    hypre_list.set("Preconditioner", ParaSails);
    hypre_list.set("SolveOrPrecondition", Preconditioner);
    hypre_list.set("SetPreconditioner", true);
    hypre_list.set("NumFunctions", 0);

    IfpHypre_Sff_->SetParameters(hypre_list);
    IfpHypre_Sff_->Initialize();
    IfpHypre_Sff_->Compute();
#endif
  }
}


void MatrixMFD_Coupled::UpdateConsistentFaceCorrection(const TreeVector& u,
        const Teuchos::Ptr<TreeVector>& Pu) {
  blockA_->UpdateConsistentFaceCorrection(*u.SubVector(0)->data(),
          Pu->SubVector(0)->data().ptr());
  blockB_->UpdateConsistentFaceCorrection(*u.SubVector(1)->data(),
          Pu->SubVector(1)->data().ptr());
}

} // namespace
} // namespace
