/*
  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)

  MatrixCoupledMFD takes two MatrixMFD objects, along with the cell coupling
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

#include "errors.hh"
#include "matrix_mfd.hh"

#include "matrix_coupled_mfd.hh"

namespace Amanzi {
namespace Operators {

MatrixCoupledMFD::MatrixCoupledMFD(Teuchos::ParameterList& plist,
        const Teuchos::RCP<const AmanziMesh::Mesh> mesh) :
    plist_(plist), mesh_(mesh) {
  InitializeFromPList_();
}


MatrixCoupledMFD::MatrixCoupledMFD(const MatrixCoupledMFD& other) :
    plist_(other.plist_), mesh_(other.mesh_) {
  InitializeFromPList_();
}


void MatrixCoupledMFD::InitializeFromPList_() {
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


void MatrixCoupledMFD::ApplyInverse(const TreeVector& X,
        const Teuchos::Ptr<TreeVector>& Y) const {
  // pull X data
  Teuchos::RCP<const CompositeVector> XA = X.SubVector(0)->data();
  Teuchos::RCP<const CompositeVector> XB = X.SubVector(1)->data();
  Epetra_MultiVector XA_c = *XA->ViewComponent("cell", false);
  Epetra_MultiVector XB_c = *XA->ViewComponent("cell", false);
  Epetra_MultiVector XA_f = *XA->ViewComponent("face", false);
  Epetra_MultiVector XB_f = *XA->ViewComponent("face", false);

  // Temporary cell and face vectors.
  Epetra_MultiVector Xc(*double_cmap_, 1);
  Epetra_MultiVector Xf(*double_fmap_, 1);

  Epetra_MultiVector Yc(*double_cmap_, 1);
  Epetra_MultiVector Yf(*double_fmap_, 1);

  Epetra_MultiVector test_res(*double_fmap_, 1);

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
  // Yc <-- A2f2c * Yf
  ierr = A2f2c_->Multiply(false, Yf, Yc);
  ASSERT(!ierr);

  // Yc -= (x_Ac,x_Bc)^T
  for (int c=0; c!=ncells; ++c){
    double a_c = -XA_c[0][c];
    double b_c = -XB_c[0][c];
    Yc.SumIntoMyValue(c, 0, 0, a_c);
    Yc.SumIntoMyValue(c, 1, 0, b_c);
  }

  // Yc <-- -Yc
  Yc.Scale(-1.0);

  // pull Y data
  Teuchos::RCP<const CompositeVector> YA = Y->SubVector(0)->data();
  Teuchos::RCP<const CompositeVector> YB = Y->SubVector(1)->data();
  Epetra_MultiVector YA_c = *YA->ViewComponent("cell", false);
  Epetra_MultiVector YB_c = *YA->ViewComponent("cell", false);
  Epetra_MultiVector YA_f = *YA->ViewComponent("face", false);
  Epetra_MultiVector YB_f = *YA->ViewComponent("face", false);

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
}


void MatrixCoupledMFD::Apply(const TreeVector& X,
        const Teuchos::Ptr<TreeVector>& Y) const {
  ASSERT(0);
}


void MatrixCoupledMFD::ComputeSchurComplement(const Epetra_MultiVector& Ccc,
                                              const Epetra_MultiVector& Dcc) {
  int ierr(0);

  const Epetra_BlockMap& cmap = mesh_->cell_map(false);
  const Epetra_BlockMap& fmap = mesh_->face_map(false);
  const Epetra_BlockMap& fmap_wghost = mesh_->face_map(true);

  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  ASSERT(Ccc.MyLength() == ncells);
  ASSERT(Dcc.MyLength() == ncells);

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
  int nfaces = 6;
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;
  const int MFD_MAX_FACES = 14;
  int faces_LID[MFD_MAX_FACES];  // Contigious memory is required.
  int faces_GID[MFD_MAX_FACES];

  // space for the assembled matrices
  Epetra_SerialDenseMatrix S2f2f(2*nfaces, 2*nfaces);
  Epetra_SerialDenseMatrix A2c2f(2, 2*nfaces);
  A2c2c_cells_Inv_.clear();
  if (is_matrix_constructed_) P2f2f_->PutScalar(0.0);

  // Assemble
  for (int c=0; c!=ncells; ++c){
    int cell_GID = cmap.GID(c);
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);

    // Invert the cell block
    double det_cell = Acc[c] * Bcc[c] - Ccc[0][c] * Dcc[0][c];
    if (det_cell != 0.) {
      cell_inv(0, 0) = Bcc[c]/det_cell;
      cell_inv(1, 1) = Acc[c]/det_cell;
      cell_inv(0, 1) = -Ccc[0][c]/det_cell;
      cell_inv(1, 0) = -Dcc[0][c]/det_cell;
    } else {
      Errors::Message m("Division by zero: determinant of the cell block is zero");
      Exceptions::amanzi_throw(m);
    }
    A2c2c_cells_Inv_.push_back(cell_inv);

    // Make the cell-local Schur complement
    for (int i=0; i!=nfaces; ++i) {
      for (int j=0; j!=nfaces; ++j) {
        S2f2f(i, j) = Aff[c](i, j) - Afc[c](i)*cell_inv(0, 0)*Acf[c](j);
      }
    }

    for (int i=0; i!=nfaces; ++i) {
      for (int j=0; j!=nfaces; ++j) {
        S2f2f(nfaces + i, nfaces + j) = Bff[c](i, j) - Bfc[c](i)*cell_inv(1, 1)*Bcf[c](j);
      }
    }

    for (int i=0; i!=nfaces; ++i) {
      for (int j=0; j!=nfaces; ++j) {
        S2f2f(i, nfaces + j) = - Afc[c](i)*cell_inv(0, 1)*Bcf[c](j);
      }
    }

    for (int i=0; i!=nfaces; ++i) {
      for (int j=0; j!=nfaces; ++j) {
        S2f2f(nfaces + i, j) = - Bfc[c](i)*cell_inv(1, 0)*Acf[c](j);
      }
    }

    // Make the local A2c2f
    for (int i=0; i!=nfaces; ++i) {
      A2c2f(0,i) =           Acf[c](i);
      A2c2f(0,i + nfaces)  = 0;
      A2c2f(1,i) =           0;
      A2c2f(1,i + nfaces)  = Bfc[c](i);
    }

    // Assemble
    int nentries = nfaces;
    for (int i=0; i!=nfaces; ++i) {
      faces_LID[i] = faces[i];
      faces_GID[i] = fmap_wghost.GID(faces_LID[i]);
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

    // -- Assemble A2f2c, which is stored as transpose?
    ierr = A2f2c_->BeginReplaceGlobalValues(cell_GID, nentries, faces_GID);
    ASSERT(!ierr);

    for (int i=0; i!=nfaces; ++i) {
      values(0,0) = A2c2f(0,i);
      values(0,1) = A2c2f(0,i + nfaces);
      values(1,0) = A2c2f(1,i);
      values(1,1) = A2c2f(1,i + nfaces);
      ierr = A2f2c_->SubmitBlockEntry(values);
      ASSERT(!ierr);
    }

    ierr = A2f2c_->EndSubmitEntries();
    ASSERT(!ierr);

  }

  // Finish assembly
  ierr = A2f2c_->FillComplete(*double_fmap_, *double_cmap_);
  ASSERT(!ierr);
  ierr = P2f2f_->GlobalAssemble();
  ASSERT(!ierr);
  is_matrix_constructed_ = true;
}


void MatrixCoupledMFD::SymbolicAssembleGlobalMatrices() {
  int ierr(0);
  const Epetra_BlockMap& cmap = mesh_->cell_map(false);
  const Epetra_BlockMap& fmap = mesh_->face_map(false);
  const Epetra_BlockMap& fmap_wghost = mesh_->face_map(true);

  // Make the maps
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
  Epetra_CrsGraph cf_graph(Copy, *double_cmap_, *double_fmap_wghost_,
                           avg_entries_row, false);
  Epetra_FECrsGraph ff_graph(Copy, *double_fmap_, 2*avg_entries_row - 1, false);

  // Workspace
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;
  const int MFD_MAX_FACES = 14;
  int faces_LID[MFD_MAX_FACES];  // Contigious memory is required.
  int faces_GID[MFD_MAX_FACES];

  // Fill the graphs
  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c=0; c!=ncells; ++c) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    for (int n=0; n!=nfaces; ++n) {
      faces_LID[n] = faces[n];
      faces_GID[n] = double_fmap_wghost_->GID(faces_LID[n]);
    }

    cf_graph.InsertMyIndices(c, nfaces, faces_LID);
    ierr = ff_graph.InsertGlobalIndices(nfaces, faces_GID, nfaces, faces_GID);
    ASSERT(!ierr);
  }

  // Assemble the graphs
  cf_graph.FillComplete(*double_fmap_, *double_cmap_);
  ierr = ff_graph.GlobalAssemble();  // Symbolic graph is complete.
  ASSERT(!ierr);

  // Create the matrices
  A2f2c_ = Teuchos::rcp(new Epetra_VbrMatrix(Copy, cf_graph));
  P2f2f_ = Teuchos::rcp(new Epetra_FEVbrMatrix(Copy, ff_graph, false));
  ierr = P2f2f_->GlobalAssemble();
  ASSERT(!ierr);
}


/* ******************************************************************
 * Initialization of the preconditioner
 ****************************************************************** */
void MatrixCoupledMFD::InitPreconditioner() {
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
void MatrixCoupledMFD::UpdatePreconditioner() {
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
    Teuchos::RCP<FunctionParameter> functs[8];
    functs[0] = Teuchos::rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetCoarsenType, 0));
    functs[1] = Teuchos::rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetPrintLevel, 0));
    functs[2] = Teuchos::rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetNumSweeps, hypre_nsmooth_));
    functs[3] = Teuchos::rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetMaxIter, hypre_ncycles_));
    functs[4] = Teuchos::rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetRelaxType, 6));
    functs[5] = Teuchos::rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetStrongThreshold, hypre_strong_threshold_));
    functs[6] = Teuchos::rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetTol, hypre_tol_));
    functs[7] = Teuchos::rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetCycleType, 1));

    Teuchos::ParameterList hypre_list;
    hypre_list.set("Preconditioner", BoomerAMG);
    hypre_list.set("SolveOrPrecondition", Preconditioner);
    hypre_list.set("SetPreconditioner", true);
    hypre_list.set("NumFunctions", 8);
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


} // namespace
} // namespace
