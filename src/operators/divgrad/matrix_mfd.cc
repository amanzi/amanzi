/*
  This is the flow component of the Amanzi code.
  License: BSD
  Authors: Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
*/

#include "errors.hh"
#include "Epetra_FECrsGraph.h"
#include "matrix_mfd.hh"

namespace Amanzi {
namespace Operators {

MatrixMFD::MatrixMFD(Teuchos::ParameterList& plist,
                     const Teuchos::RCP<const AmanziMesh::Mesh> mesh) :
    plist_(plist), mesh_(mesh) {
  std::string methodstring = plist.get<string>("MFD method");

  if (methodstring == "polyhedra") {
    method_ = MFD_POLYHEDRA;
  } else if (methodstring == "polyhedra monotone") {
    method_ = MFD_POLYHEDRA_MONOTONE;
  } else if (methodstring == "hexahedra monotone") {
    method_ = MFD_HEXAHEDRA_MONOTONE;
  }

  std::string precmethodstring = plist.get<string>("preconditioner", "Trilinos ML");
  
  if (precmethodstring == "Trilinos ML") {
    prec_method_ = TRILINOS_ML;
  } else if (precmethodstring == "HYPRE AMG") {
    prec_method_ = HYPRE_AMG;
  }
}

// main computational methods
/* ******************************************************************
 * Calculate elemental inverse mass matrices.
 * WARNING: The original Aff_ matrices are destroyed.
 ****************************************************************** */
void MatrixMFD::CreateMFDmassMatrices(const Teuchos::Ptr<std::vector<WhetStone::Tensor> >& K) {
  int dim = mesh_->space_dimension();

  // TODO: fix -- MFD3D should NOT need a non-const mesh
  Teuchos::RCP<AmanziMesh::Mesh> nc_mesh = Teuchos::rcp_const_cast<AmanziMesh::Mesh>(mesh_);
  WhetStone::MFD3D mfd(nc_mesh);
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  Mff_cells_.clear();

  int ok;
  nokay_ = npassed_ = 0;

  WhetStone::Tensor Kc;
  if (K == Teuchos::null) {
    Kc.init(mesh_->space_dimension(), 1);
    Kc(0,0) = 1.0;
  }

  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c=0; c!=ncells; ++c) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    Teuchos::SerialDenseMatrix<int, double> Mff(nfaces, nfaces);

    if (K != Teuchos::null) {
      Kc = (*K)[c];
    }

    if (method_ == MFD_HEXAHEDRA_MONOTONE) {
      if ((nfaces == 6 && dim == 3) || (nfaces == 4 && dim == 2)) {
        ok = mfd.darcy_mass_inverse_hex(c, Kc, Mff);
      } else {
        ok = mfd.darcy_mass_inverse(c, Kc, Mff);
      }
    } else if (method_ == MFD_TWO_POINT_FLUX) {
      ok = mfd.darcy_mass_inverse_diagonal(c, Kc, Mff);
    } else if (method_ == MFD_SUPPORT_OPERATOR) {
      ok = mfd.darcy_mass_inverse_SO(c, Kc, Mff);
    } else if (method_ == MFD_OPTIMIZED) {
      ok = mfd.darcy_mass_inverse_optimized(c, Kc, Mff);
    } else {
      ok = mfd.darcy_mass_inverse(c, Kc, Mff);
    }

    Mff_cells_.push_back(Mff);

    if (ok == WhetStone::WHETSTONE_ELEMENTAL_MATRIX_FAILED) {
      Errors::Message msg("Matrix_MFD: unexpected failure of LAPACK in WhetStone.");
      Exceptions::amanzi_throw(msg);
    }

    if (ok == WhetStone::WHETSTONE_ELEMENTAL_MATRIX_OK) nokay_++;
    if (ok == WhetStone::WHETSTONE_ELEMENTAL_MATRIX_PASSED) npassed_++;
  }

  // sum up the numbers across processors
  int nokay_tmp = nokay_, npassed_tmp = npassed_;
  mesh_->get_comm()->SumAll(&nokay_tmp, &nokay_, 1);
  mesh_->get_comm()->SumAll(&npassed_tmp, &npassed_, 1);
}


/* ******************************************************************
 * Calculate elemental stiffness matrices.
 ****************************************************************** */
void MatrixMFD::CreateMFDstiffnessMatrices(const CompositeVector& Krel) {
  int dim = mesh_->space_dimension();

  WhetStone::MFD3D mfd(mesh_);
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  Aff_cells_.clear();
  Afc_cells_.clear();
  Acf_cells_.clear();
  Acc_cells_.clear();

  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c=0; c!=ncells; ++c) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    Teuchos::SerialDenseMatrix<int, double>& Mff = Mff_cells_[c];
    Teuchos::SerialDenseMatrix<int, double> Bff(nfaces,nfaces);
    Epetra_SerialDenseVector Bcf(nfaces), Bfc(nfaces);

    if (!Krel.has_component("face")) {
      for (int n=0; n!=nfaces; ++n) {
        for (int m=0; m!=nfaces; ++m) {
          Bff(m, n) = Mff(m,n) * Krel("cell",0,c);
        }
      }
    } else if (!Krel.has_component("cell")) {
      for (int n=0; n!=nfaces; ++n) {
        for (int m=0; m!=nfaces; ++m) {
          Bff(m, n) = Mff(m,n) * Krel("face",0,faces[m]);
        }
      }
    } else {
      for (int n=0; n!=nfaces; ++n) {
        for (int m=0; m!=nfaces; ++m) {
          Bff(m, n) = Mff(m,n) * Krel("cell",0,c) * Krel("face",0,faces[m]);
        }
      }
    }

    double matsum = 0.0;  // elimination of mass matrix
    for (int n=0; n!=nfaces; ++n) {
      double rowsum = 0.0, colsum = 0.0;
      for (int m=0; m!=nfaces; ++m) {
        colsum += Bff(m, n);
        rowsum += Bff(n, m);
      }
      Bcf(n) = -colsum;
      Bfc(n) = -rowsum;
      matsum += colsum;
    }

    Aff_cells_.push_back(Bff);  // This the only place where memory can be allocated.
    Afc_cells_.push_back(Bfc);
    Acf_cells_.push_back(Bcf);
    Acc_cells_.push_back(matsum);
  }
}


/* ******************************************************************
 * May be used in the future.
 ****************************************************************** */
void MatrixMFD::RescaleMFDstiffnessMatrices(const Epetra_Vector& old_scale,
        const Epetra_Vector& new_scale) {

  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c=0; c!=ncells; ++c) {
    Teuchos::SerialDenseMatrix<int, double>& Bff = Aff_cells_[c];
    Epetra_SerialDenseVector& Bcf = Acf_cells_[c];

    int n = Bff.numRows();
    double scale = old_scale[c] / new_scale[c];

    for (int i=0; i<n; i++) {
      for (int j=0; j<n; j++) Bff(i, j) *= scale;
      Bcf(i) *= scale;
    }
    Acc_cells_[c] *= scale;
  }
}


/* ******************************************************************
 * Simply allocates memory.
 ****************************************************************** */
void MatrixMFD::CreateMFDrhsVectors() {
  Ff_cells_.clear();
  Fc_cells_.clear();

  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  for (int c=0; c!=ncells; ++c) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    Epetra_SerialDenseVector Ff(nfaces);  // Entries are initilaized to 0.0.
    double Fc = 0.0;

    Ff_cells_.push_back(Ff);
    Fc_cells_.push_back(Fc);
  }
}

/* ******************************************************************
 *  Create work vectors for Apply-ing the operator from/to Epetra_Vectors
 * (Supervectors) instead of CompositeVectors -- for use with AztecOO.
 * This should likely only be called with a non-ghosted sample vector.
 ****************************************************************** */
// void MatrixMFD::InitializeSuperVecs(const CompositeVector& sample) {
//   vector_x_ = Teuchos::rcp(new CompositeVector(sample));
//   vector_y_ = Teuchos::rcp(new CompositeVector(sample));
//   supermap_ = vector_x_->supermap();
// }

/* ******************************************************************
 * Applies boundary conditions to elemental stiffness matrices and
 * creates elemental rigth-hand-sides.
 ****************************************************************** */
void MatrixMFD::ApplyBoundaryConditions(const std::vector<Matrix_bc>& bc_markers,
        const std::vector<double>& bc_values) {
  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  for (int c=0; c!=ncells; ++c) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    Teuchos::SerialDenseMatrix<int, double>& Bff = Aff_cells_[c];  // B means elemental.
    Epetra_SerialDenseVector& Bfc = Afc_cells_[c];
    Epetra_SerialDenseVector& Bcf = Acf_cells_[c];

    Epetra_SerialDenseVector& Ff = Ff_cells_[c];
    double& Fc = Fc_cells_[c];

    for (int n=0; n!=nfaces; ++n) {
      int f=faces[n];
      if (bc_markers[f] == MFD_BC_DIRICHLET) {
        for (int m=0; m!=nfaces; ++m) {
          Ff[m] -= Bff(m, n) * bc_values[f];
          Bff(n, m) = Bff(m, n) = 0.0;
        }
        Fc -= Bcf(n) * bc_values[f];
        Bcf(n) = Bfc(n) = 0.0;

        Bff(n, n) = 1.0;
        Ff[n] = bc_values[f];
      } else if (bc_markers[f] == MFD_BC_FLUX) {
        Ff[n] -= bc_values[f] * mesh_->face_area(f);
      }
    }
  }
}


/* ******************************************************************
 * Initialize Trilinos matrices. It must be called only once.
 * If matrix is non-symmetric, we generate transpose of the matrix
 * block Afc_ to reuse cf_graph; otherwise, pointer Afc_ = Acf_.
 ****************************************************************** */
void MatrixMFD::SymbolicAssembleGlobalMatrices() {
  const Epetra_Map& cmap = mesh_->cell_map(false);
  const Epetra_Map& fmap = mesh_->face_map(false);
  const Epetra_Map& fmap_wghost = mesh_->face_map(true);

  int avg_entries_row = (mesh_->space_dimension() == 2) ? MFD_QUAD_FACES : MFD_HEX_FACES;
  Epetra_CrsGraph cf_graph(Copy, cmap, fmap_wghost, avg_entries_row, false);  // FIX (lipnikov@lanl.gov)
  Epetra_FECrsGraph ff_graph(Copy, fmap, 2*avg_entries_row);

  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;
  int faces_LID[MFD_MAX_FACES];  // Contigious memory is required.
  int faces_GID[MFD_MAX_FACES];

  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c=0; c!=ncells; ++c) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    for (int n=0; n!=nfaces; ++n) {
      faces_LID[n] = faces[n];
      faces_GID[n] = fmap_wghost.GID(faces_LID[n]);
    }
    cf_graph.InsertMyIndices(c, nfaces, faces_LID);
    ff_graph.InsertGlobalIndices(nfaces, faces_GID, nfaces, faces_GID);
  }
  cf_graph.FillComplete(fmap, cmap);
  ff_graph.GlobalAssemble();  // Symbolic graph is complete.

  // create global matrices
  Acc_ = Teuchos::rcp(new Epetra_Vector(cmap));
  Acf_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, cf_graph));
  Aff_ = Teuchos::rcp(new Epetra_FECrsMatrix(Copy, ff_graph));
  Sff_ = Teuchos::rcp(new Epetra_FECrsMatrix(Copy, ff_graph));
  Aff_->GlobalAssemble();
  Sff_->GlobalAssemble();

  if (flag_symmetry_) {
    Afc_ = Acf_;
  } else {
    Afc_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, cf_graph));
  }

  std::vector<std::string> names(2);
  names[0] = "cell"; names[1] = "face";

  std::vector<AmanziMesh::Entity_kind> locations(2);
  locations[0] = AmanziMesh::CELL; locations[1] = AmanziMesh::FACE;

  std::vector<int> num_dofs(2,1);
  rhs_ = Teuchos::rcp(new CompositeVector(mesh_, names, locations, num_dofs, true));
  rhs_->CreateData();
}


/* ******************************************************************
 * Convert elemental mass matrices into stiffness matrices and
 * assemble them into four global matrices.
 * We need an auxiliary GHOST-based vector to assemble the RHS.
 ****************************************************************** */
void MatrixMFD::AssembleGlobalMatrices() {
  Aff_->PutScalar(0.0);

  const Epetra_Map& fmap_wghost = mesh_->face_map(true);
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;
  int faces_LID[MFD_MAX_FACES];
  int faces_GID[MFD_MAX_FACES];

  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c=0; c!=ncells; ++c) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    for (int n=0; n!=nfaces; ++n) {
      faces_LID[n] = faces[n];
      faces_GID[n] = fmap_wghost.GID(faces_LID[n]);
    }
    (*Acc_)[c] = Acc_cells_[c];
    (*Acf_).ReplaceMyValues(c, nfaces, Acf_cells_[c].Values(), faces_LID);
    (*Aff_).SumIntoGlobalValues(nfaces, faces_GID, Aff_cells_[c].values());

    if (!flag_symmetry_)
      (*Afc_).ReplaceMyValues(c, nfaces, Afc_cells_[c].Values(), faces_LID);
  }
  (*Aff_).GlobalAssemble();

  // We repeat some of the loops for code clarity.
  rhs_->ViewComponent("face", true)->PutScalar(0.0);
  for (int c=0; c!=ncells; ++c) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    (*rhs_)("cell",c) = Fc_cells_[c];

    for (int n=0; n!=nfaces; ++n) {
      int f = faces[n];
      (*rhs_)("face",f) += Ff_cells_[c][n];
    }
  }
  rhs_->GatherGhostedToMaster("face");
}


/* ******************************************************************
 * Compute the face Schur complement of 2x2 block matrix.
 ****************************************************************** */
void MatrixMFD::ComputeSchurComplement(const std::vector<Matrix_bc>& bc_markers,
        const std::vector<double>& bc_values) {
  Sff_->PutScalar(0.0);

  AmanziMesh::Entity_ID_List faces_LID;
  std::vector<int> dirs;
  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);

  for (int c=0; c!=ncells; ++c) {
    mesh_->cell_get_faces_and_dirs(c, &faces_LID, &dirs);
    int nfaces = faces_LID.size();
    Epetra_SerialDenseMatrix Schur(nfaces, nfaces);

    Epetra_SerialDenseVector& Bcf = Acf_cells_[c];
    Epetra_SerialDenseVector& Bfc = Afc_cells_[c];

    for (int n=0; n!=nfaces; ++n) {
      for (int m=0; m!=nfaces; ++m) {
        Schur(n, m) = Aff_cells_[c](n, m) - Bfc[n] * Bcf[m] / (*Acc_)[c];
      }
    }

    for (int n=0; n!=nfaces; ++n) {  // Symbolic boundary conditions
      int f=faces_LID[n];
      if (bc_markers[f] == MFD_BC_DIRICHLET) {
        for (int m=0; m!=nfaces; ++m) Schur(n, m) = Schur(m, n) = 0.0;
        Schur(n, n) = 1.0;
      }
    }

    Epetra_IntSerialDenseVector faces_GID(nfaces);
    for (int n=0; n!=nfaces; ++n) faces_GID[n] = (*Acf_).ColMap().GID(faces_LID[n]);
    (*Sff_).SumIntoGlobalValues(faces_GID, Schur);
  }
  (*Sff_).GlobalAssemble();
}

/* ******************************************************************
 * Parallel matvec product A * X.
 ****************************************************************** */
// int MatrixMFD::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {
//   vector_x_->DataFromSuperVector(*X(0));
//   vector_y_->DataFromSuperVector(*Y(0));
//   Apply(*vector_x_, vector_y_);
// }

// int MatrixMFD::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {
//   vector_x_->DataFromSuperVector(*X(0));
//   vector_y_->DataFromSuperVector(*Y(0));
//   ApplyInverse(*vector_x_, vector_y_);
// }

/* ******************************************************************
 * Parallel matvec product A * X.
 ****************************************************************** */
void MatrixMFD::Apply(const CompositeVector& X,
                     const Teuchos::RCP<CompositeVector>& Y) const {
  int ierr;

  // Face unknowns:  Yf = Aff_ * Xf + Afc_ * Xc
  ierr = (*Aff_).Multiply(false, *X.ViewComponent("face",false),
                          *Y->ViewComponent("face", false));

  Epetra_MultiVector Tf(*Y->ViewComponent("face", false));
  ierr |= (*Afc_).Multiply(true, *X.ViewComponent("cell",false), Tf);  // Afc_ is kept in transpose form
  Y->ViewComponent("face",false)->Update(1.0, Tf, 1.0);

  // Cell unknowns:  Yc = Acf_ * Xf + Acc_ * Xc
  ierr |= (*Acf_).Multiply(false, *X.ViewComponent("face", false),
                           *Y->ViewComponent("cell", false));  // It performs the required parallel communications.
  ierr |= Y->ViewComponent("cell", false)->Multiply(1.0, *Acc_, *X.ViewComponent("cell", false), 1.0);

  if (ierr) {
    Errors::Message msg("MatrixMFD::Apply has failed to calculate Y = inv(A) * X.");
    Exceptions::amanzi_throw(msg);
  }
}


/* ******************************************************************
 * The OWNED cell-based and face-based d.o.f. are packed together into
 * the X and Y Epetra vectors, with the cell-based in the first part.
 *
 * WARNING: When invoked by AztecOO the arguments X and Y may be
 * aliased: possibly the same object or different views of the same
 * underlying data. Thus, we do not assign to Y until the end.
 *
 * NOTE however that this is broken for AztecOO since we use CompositeVectors,
 * not Epetra_MultiVectors.
 *
 ****************************************************************** */
void MatrixMFD::ApplyInverse(const CompositeVector& X,
                     const Teuchos::RCP<CompositeVector>& Y) const {
  // Temporary cell and face vectors.
  Epetra_MultiVector Tc(*Y->ViewComponent("cell", false));
  Epetra_MultiVector Tf(*Y->ViewComponent("face", false));

  // FORWARD ELIMINATION:  Tf = Xf - Afc_ inv(Acc_) Xc
  int ierr;

  ierr  = Tc.ReciprocalMultiply(1.0, *Acc_, *X.ViewComponent("cell", false), 0.0);
  ierr |= (*Afc_).Multiply(true, Tc, Tf);  // Afc_ is kept in transpose form
  Tf.Update(1.0, *X.ViewComponent("face", false), -1.0);

  // Solve the Schur complement system Sff_ * Yf = Tf.
  if (prec_method_ == TRILINOS_ML) {
    ierr |= ml_prec_->ApplyInverse(Tf, *Y->ViewComponent("face", false));
  } else if (prec_method_ == HYPRE_AMG) {
    IfpHypre_Sff_->ApplyInverse(Tf, *Y->ViewComponent("face", false));
  }

  // BACKWARD SUBSTITUTION:  Yc = inv(Acc_) (Xc - Acf_ Yf)
  ierr |= (*Acf_).Multiply(false, *Y->ViewComponent("face", false), Tc);  // It performs the required parallel communications.
  Tc.Update(1.0, *X.ViewComponent("cell", false), -1.0);
  ierr |= Y->ViewComponent("cell", false)->ReciprocalMultiply(1.0, *Acc_, Tc, 0.0);

  if (ierr) {
    Errors::Message msg("MatrixMFD::ApplyInverse has failed in calculating y = A*x.");
    Exceptions::amanzi_throw(msg);
  }
}


/* ******************************************************************
 * Linear algebra operations with matrices: r = f - A * x
 ****************************************************************** */
void MatrixMFD::ComputeResidual(const CompositeVector& solution,
          const Teuchos::RCP<CompositeVector>& residual) const {
  Apply(solution, residual);
  residual->Update(1.0, *rhs_, -1.0);
}


/* ******************************************************************
 * Linear algebra operations with matrices: r = A * x - f
 ****************************************************************** */
void MatrixMFD::ComputeNegativeResidual(const CompositeVector& solution,
        const Teuchos::RCP<CompositeVector>& residual) const {
  Apply(solution, residual);
  residual->Update(-1.0, *rhs_, 1.0);
}


/* ******************************************************************
 * Initialization of the preconditioner
 ****************************************************************** */
void MatrixMFD::InitPreconditioner(Teuchos::ParameterList& prec_plist) {
  if (prec_method_ == TRILINOS_ML) {
    ml_plist_ =  prec_plist.sublist("ML Parameters");
    ml_prec_ = Teuchos::rcp(new ML_Epetra::MultiLevelPreconditioner(*Sff_, ml_plist_, false));
  } else if (prec_method_ == HYPRE_AMG) {
    // read some boomer amg parameters 
    hypre_plist_ = prec_plist.sublist("HYPRE AMG Parameters");
    hypre_ncycles_ = hypre_plist_.get<int>("number of cycles",5);
    hypre_nsmooth_ = hypre_plist_.get<int>("number of smoothing iterations",3);
    hypre_tol_ = hypre_plist_.get<double>("tolerance",0.0);
    hypre_strong_threshold_ = hypre_plist_.get<double>("strong threshold",0.25);
  }
}

/* ******************************************************************
 * Rebuild preconditioner.
 ****************************************************************** */
void MatrixMFD::UpdatePreconditioner() {
  if (prec_method_ == TRILINOS_ML) {
    if (ml_prec_->IsPreconditionerComputed()) ml_prec_->DestroyPreconditioner();
    ml_prec_->SetParameterList(ml_plist_);
    ml_prec_->ComputePreconditioner();
  } else if (prec_method_ == HYPRE_AMG) {
    IfpHypre_Sff_ = Teuchos::rcp(new Ifpack_Hypre(&*Sff_));
    Teuchos::RCP<FunctionParameter> functs[10];
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
  }
}



/* ******************************************************************
 * WARNING: Routines requires original mass matrices (Aff_cells_), i.e.
 * before boundary conditions were imposed.
 *
 * WARNING: Since diffusive flux is not continuous, we derive it only
 * once (using flag) and in exactly the same manner as in routine
 * Flow_PK::addGravityFluxes_DarcyFlux.
 ****************************************************************** */
void MatrixMFD::DeriveFlux(const CompositeVector& solution,
                           const Teuchos::RCP<CompositeVector>& flux) const {

  AmanziMesh::Entity_ID_List faces;
  std::vector<double> dp;
  std::vector<int> dirs;

  solution.ScatterMasterToGhosted("face");
  flux->PutScalar(0.);

  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  int nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  std::vector<int> flag(solution.ViewComponent("face", true)->MyLength(), 0);

  for (int c=0; c!=ncells; ++c) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    dp.resize(nfaces);
    for (int n=0; n!=nfaces; ++n) {
      int f = faces[n];
      dp[n] = solution("cell", c) - solution("face", f);
    }

    for (int n=0; n!=nfaces; ++n) {
      int f = faces[n];
      if (f < nfaces_owned && !flag[f]) {
        double s = 0.0;
        for (int m=0; m!=nfaces; ++m) s += Aff_cells_[c](n, m) * dp[m];
        (*flux)("face", 0, f) = s * dirs[n];
        flag[f] = 1;
      }
    }
  }
}


/* ******************************************************************
 * Derive Darcy velocity in cells.
 * WARNING: It cannot be consistent with the Darcy flux.
 ****************************************************************** */
void MatrixMFD::DeriveCellVelocity(const CompositeVector& flux,
        const Teuchos::RCP<CompositeVector>& velocity) const {

  Teuchos::LAPACK<int, double> lapack;

  int dim = mesh_->space_dimension();
  Teuchos::SerialDenseMatrix<int, double> matrix(dim, dim);
  double rhs_cell[dim];

  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c=0; c!=ncells_owned; ++c) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    for (int i=0; i!=dim; ++i) rhs_cell[i] = 0.0;
    matrix.putScalar(0.0);

    for (int n=0; n!=nfaces; ++n) {  // populate least-square matrix
      int f = faces[n];
      const AmanziGeometry::Point& normal = mesh_->face_normal(f);
      double area = mesh_->face_area(f);

      for (int i=0; i!=dim; ++i) {
        rhs_cell[i] += normal[i] * flux("face",0,f);
        matrix(i,i) += normal[i] * normal[i];
        for (int j=i+1; j!=dim; ++j) {
          matrix(j,i) = matrix(i,j) += normal[i] * normal[j];
        }
      }
    }

    int info;
    lapack.POSV('U', dim, 1, matrix.values(), dim, rhs_cell, dim, &info);

    for (int i=0; i!=dim; ++i) (*velocity)("cell",i,c) = rhs_cell[i];
  }
}

}  // namespace AmanziFlow
}  // namespace Amanzi
