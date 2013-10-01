/*
  This is the flow component of the Amanzi code.
  License: BSD
  Authors: Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
*/

#include "errors.hh"
#include "Epetra_FECrsGraph.h"
#include "EpetraExt_RowMatrixOut.h"
#include "MatrixMFD.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
 * Constructor
 ****************************************************************** */
MatrixMFD::MatrixMFD(Teuchos::ParameterList& plist,
                     const Teuchos::RCP<const AmanziMesh::Mesh> mesh) :
    plist_(plist),
    mesh_(mesh),
    flag_symmetry_(false),
    assembled_operator_(false),
    assembled_schur_(false),
    method_(MFD3D_NULL),
    hypre_ncycles_(5),
    hypre_nsmooth_(3),
    hypre_tol_(0.0),
    hypre_strong_threshold_(0.25) {
  InitializeFromPList_();
}


/* ******************************************************************
 * Copy constructor
 ****************************************************************** */
MatrixMFD::MatrixMFD(const MatrixMFD& other) :
    plist_(other.plist_),
    mesh_(other.mesh_),
    flag_symmetry_(other.flag_symmetry_),
    assembled_operator_(false),
    assembled_schur_(false),
    method_(other.method_),
    hypre_ncycles_(5),
    hypre_nsmooth_(3),
    hypre_tol_(0.0),
    hypre_strong_threshold_(0.25) {
  InitializeFromPList_();
}


/* ******************************************************************
 * Initialization of method, solver, etc.
 ****************************************************************** */
void MatrixMFD::InitializeFromPList_() {
  std::string methodstring = plist_.get<string>("MFD method");
  method_ = MFD3D_NULL;

  // standard MFD
  if (methodstring == "polyhedra") {
    method_ = MFD3D_POLYHEDRA;
  } else if (methodstring == "polyhedra scaled") {
    method_ = MFD3D_POLYHEDRA_SCALED;
  } else if (methodstring == "optimized") {
    method_ = MFD3D_OPTIMIZED;
  } else if (methodstring == "optimized scaled") {
    method_ = MFD3D_OPTIMIZED_SCALED;
  } else if (methodstring == "hexahedra monotone") {
    method_ = MFD3D_HEXAHEDRA_MONOTONE;
  } else if (methodstring == "two point flux") {
    method_ = MFD3D_TWO_POINT_FLUX;
  } else if (methodstring == "support operator") {
    method_ = MFD3D_SUPPORT_OPERATOR;
  } else {
	Errors::Message msg("MatrixMFD: unexpected discretization methods");
	Exceptions::amanzi_throw(msg);
  }

  // method for inversion
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


// main computational methods
/* ******************************************************************
 * Calculate elemental inverse mass matrices.
 ****************************************************************** */
void MatrixMFD::CreateMFDmassMatrices(
    const Teuchos::Ptr<std::vector<WhetStone::Tensor> >& K) {
  // tag global matrices as invalid
  assembled_schur_ = false;
  assembled_operator_ = false;

  int dim = mesh_->space_dimension();
  WhetStone::MFD3D_Diffusion mfd(mesh_);
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

    WhetStone::DenseMatrix Mff(nfaces, nfaces);

    if (K != Teuchos::null) {
      Kc = (*K)[c];
    }

    if (method_ == MFD3D_POLYHEDRA_SCALED) {
      ok = mfd.MassMatrixInverseScaled(c, Kc, Mff);
    } else if (method_ == MFD3D_POLYHEDRA) {
      ok = mfd.MassMatrixInverse(c, Kc, Mff);
    } else if (method_ == MFD3D_OPTIMIZED_SCALED) {
      ok = mfd.MassMatrixInverseOptimizedScaled(c, Kc, Mff);
    } else if (method_ == MFD3D_OPTIMIZED) {
      ok = mfd.MassMatrixInverseOptimized(c, Kc, Mff);
    } else if (method_ == MFD3D_HEXAHEDRA_MONOTONE) {
      if ((nfaces == 6 && dim == 3) || (nfaces == 4 && dim == 2))
        ok = mfd.MassMatrixInverseHex(c, Kc, Mff);
      else
        ok = mfd.MassMatrixInverse(c, Kc, Mff);
    } else if (method_ == MFD3D_TWO_POINT_FLUX) {
      ok = mfd.MassMatrixInverseDiagonal(c, Kc, Mff);
    } else if (method_ == MFD3D_SUPPORT_OPERATOR) {
      ok = mfd.MassMatrixInverseSO(c, Kc, Mff);
    } else {
      Errors::Message msg("Flow PK: unexpected discretization methods (contact lipnikov@lanl.gov).");
      Exceptions::amanzi_throw(msg);
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
void MatrixMFD::CreateMFDstiffnessMatrices(
    const Teuchos::Ptr<const CompositeVector>& Krel) {
  // tag global matrices as invalid
  assembled_schur_ = false;
  assembled_operator_ = false;

  int dim = mesh_->space_dimension();
  WhetStone::MFD3D_Diffusion mfd(mesh_);
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

    WhetStone::DenseMatrix& Mff = Mff_cells_[c];
    Teuchos::SerialDenseMatrix<int, double> Bff(nfaces,nfaces);
    Epetra_SerialDenseVector Bcf(nfaces), Bfc(nfaces);

    if (Krel == Teuchos::null ||
        (!Krel->HasComponent("cell") && !Krel->HasComponent("face"))) {
      for (int n=0; n!=nfaces; ++n) {
        for (int m=0; m!=nfaces; ++m) {
          Bff(m, n) = Mff(m,n);
        }
      }
    } else if (Krel->HasComponent("cell") && !Krel->HasComponent("face")) {
      const Epetra_MultiVector& Krel_c = *Krel->ViewComponent("cell",false);

      for (int n=0; n!=nfaces; ++n) {
        for (int m=0; m!=nfaces; ++m) {
          Bff(m, n) = Mff(m,n) * Krel_c[0][c];
        }
      }
    } else if (!Krel->HasComponent("cell") && Krel->HasComponent("face")) {
      Krel->ScatterMasterToGhosted("face");
      const Epetra_MultiVector& Krel_f = *Krel->ViewComponent("face",true);

      for (int n=0; n!=nfaces; ++n) {
        for (int m=0; m!=nfaces; ++m) {
          Bff(m, n) = Mff(m,n) * Krel_f[0][faces[m]];
        }
      }
    } else if (Krel->HasComponent("cell") && Krel->HasComponent("face")) {
      Krel->ScatterMasterToGhosted("face");
      const Epetra_MultiVector& Krel_f = *Krel->ViewComponent("face",true);
      const Epetra_MultiVector& Krel_c = *Krel->ViewComponent("cell",false);

      for (int n=0; n!=nfaces; ++n) {
        for (int m=0; m!=nfaces; ++m) {
          Bff(m, n) = Mff(m,n) * Krel_c[0][c] * Krel_f[0][faces[m]];
        }
      }
    }

	double matsum = 0.0;
	for (int n=0; n!=nfaces; ++n) {
	  double rowsum = 0.0;
      double colsum = 0.0;

	  for (int m=0; m!=nfaces; ++m) {
        colsum += Bff(m, n);
		rowsum += Bff(n, m);
	  }

	  Bcf(n) = -colsum;
	  Bfc(n) = -rowsum;
	  matsum += colsum;
	}

    Aff_cells_.push_back(Bff);
    Afc_cells_.push_back(Bfc);
    Acf_cells_.push_back(Bcf);
    Acc_cells_.push_back(matsum);
  }
}


/* ******************************************************************
 * Create elemental rhs vectors.
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
 * Applies boundary conditions to elemental stiffness matrices and
 * adds contributions to elemental rigth-hand-sides.
 ****************************************************************** */
void MatrixMFD::ApplyBoundaryConditions(const std::vector<MatrixBC>& bc_markers,
        const std::vector<double>& bc_values) {
  // tag global matrices as invalid
  assembled_schur_ = false;
  assembled_operator_ = false;

  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  int nfaces = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  AmanziMesh::Entity_ID_List faces;
  AmanziMesh::Entity_ID_List cells;
  std::vector<int> dirs;

  for (int c=0; c!=ncells; ++c) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    Teuchos::SerialDenseMatrix<int, double>& Bff = Aff_cells_[c];
    Epetra_SerialDenseVector& Bfc = Afc_cells_[c];
    Epetra_SerialDenseVector& Bcf = Acf_cells_[c];

    Epetra_SerialDenseVector& Ff = Ff_cells_[c];
    double& Fc = Fc_cells_[c];

    for (int n=0; n!=nfaces; ++n) {
      int f=faces[n];
      if (bc_markers[f] == MATRIX_BC_DIRICHLET) {
        for (int m=0; m!=nfaces; ++m) {
          Ff[m] -= Bff(m, n) * bc_values[f];
          Bff(n, m) = Bff(m, n) = 0.0;
        }
        Fc -= Bcf(n) * bc_values[f];
        Bcf(n) = Bfc(n) = 0.0;

        Bff(n, n) = 1.0;
        Ff[n] = bc_values[f];
      } else if (bc_markers[f] == MATRIX_BC_FLUX) {
        Ff[n] -= bc_values[f] * mesh_->face_area(f);
      }
    }
  }
}


/* ******************************************************************
 * Initialize global matrices.
 *
 * This likely should only be called once.
 * If matrix is non-symmetric, we generate transpose of the matrix
 * block Afc_ to reuse cf_graph; otherwise, pointer Afc_ = Acf_.
 ****************************************************************** */
void MatrixMFD::SymbolicAssembleGlobalMatrices() {
  const Epetra_Map& cmap = mesh_->cell_map(false);
  const Epetra_Map& fmap = mesh_->face_map(false);
  const Epetra_Map& fmap_wghost = mesh_->face_map(true);

  int avg_entries_row = (mesh_->space_dimension() == 2) ? MFD_QUAD_FACES : MFD_HEX_FACES;

  // allocate the graphs
  Teuchos::RCP<Epetra_CrsGraph> cf_graph =
      Teuchos::rcp(new Epetra_CrsGraph(Copy, cmap, fmap_wghost, avg_entries_row, false));
  Teuchos::RCP<Epetra_FECrsGraph> ff_graph =
      Teuchos::rcp(new Epetra_FECrsGraph(Copy, fmap, 2*avg_entries_row));

  // fill the graphs
  FillMatrixGraphs_(cf_graph.ptr(), ff_graph.ptr());

  // assemble the graphs
  int ierr = cf_graph->FillComplete(fmap, cmap);
  ASSERT(!ierr);
  ierr = ff_graph->GlobalAssemble();  // Symbolic graph is complete.
  ASSERT(!ierr);

  // allocate the matrices
  CreateMatrices_(*cf_graph, *ff_graph);
}


/* ******************************************************************
 * Fill sparsity structure graphs of global matrices (only done once)
 ****************************************************************** */
void MatrixMFD::FillMatrixGraphs_(const Teuchos::Ptr<Epetra_CrsGraph> cf_graph,
          const Teuchos::Ptr<Epetra_FECrsGraph> ff_graph) {
  const Epetra_Map& cmap = mesh_->cell_map(false);
  const Epetra_Map& fmap = mesh_->face_map(false);
  const Epetra_Map& fmap_wghost = mesh_->face_map(true);

  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;
  int faces_LID[MFD_MAX_FACES];  // Contigious memory is required.
  int faces_GID[MFD_MAX_FACES];

  // fill the graphs
  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c=0; c!=ncells; ++c) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    for (int n=0; n!=nfaces; ++n) {
      faces_LID[n] = faces[n];
      faces_GID[n] = fmap_wghost.GID(faces_LID[n]);
    }
    cf_graph->InsertMyIndices(c, nfaces, faces_LID);
    ff_graph->InsertGlobalIndices(nfaces, faces_GID, nfaces, faces_GID);
  }
}


/* ******************************************************************
 * Allocate global matrices
 ****************************************************************** */
void MatrixMFD::CreateMatrices_(const Epetra_CrsGraph& cf_graph,
        const Epetra_FECrsGraph& ff_graph) {
  // create global matrices
  const Epetra_Map& cmap = mesh_->cell_map(false);
  Acc_ = Teuchos::rcp(new Epetra_Vector(cmap));

  Acf_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, cf_graph));
  Aff_ = Teuchos::rcp(new Epetra_FECrsMatrix(Copy, ff_graph));
  Sff_ = Teuchos::rcp(new Epetra_FECrsMatrix(Copy, ff_graph));
  Aff_->GlobalAssemble();
  Sff_->GlobalAssemble();

  if (symmetric()) {
    Afc_ = Acf_;
  } else {
    Afc_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, cf_graph));
  }

  // create the RHS
  std::vector<std::string> names(2);
  names[0] = "cell";
  names[1] = "face";

  std::vector<AmanziMesh::Entity_kind> locations(2);
  locations[0] = AmanziMesh::CELL;
  locations[1] = AmanziMesh::FACE;

  std::vector<int> num_dofs(2,1);
  CompositeVectorSpace space;
  space.SetMesh(mesh_)->SetGhosted()->SetComponents(names,locations,num_dofs);
  rhs_ = Teuchos::rcp(new CompositeVector(space));
}


/* ******************************************************************
 * Assemble global matrices from local matrices.
 *
 * Convert elemental mass matrices into stiffness matrices and
 * assemble them into four global matrices.
 * We need an auxiliary GHOST-based vector to assemble the RHS.
 *
 * Precondition: SymbolicAssembleGlobalMatrices() has been called.
 ****************************************************************** */
void MatrixMFD::AssembleGlobalMatrices() {
  ASSERT(Aff_.get()); // precondition: matrices have been created

  // reinitialize to zero if adding
  Aff_->PutScalar(0.0);

  rhs_->ViewComponent("face", true)->PutScalar(0.0);
  Epetra_MultiVector& rhs_c = *rhs_->ViewComponent("cell", false);
  Epetra_MultiVector& rhs_f = *rhs_->ViewComponent("face", true);

  // loop over cells and fill
  const Epetra_Map& fmap_wghost = mesh_->face_map(true);
  int faces_LID[MFD_MAX_FACES];
  int faces_GID[MFD_MAX_FACES];

  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c=0; c!=ncells; ++c) {
    AmanziMesh::Entity_ID_List faces;
    std::vector<int> dirs;
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    for (int n=0; n!=nfaces; ++n) {
      faces_LID[n] = faces[n];
      faces_GID[n] = fmap_wghost.GID(faces_LID[n]);
    }

    // assemble matrices
    (*Acc_)[c] = Acc_cells_[c];
    (*Acf_).ReplaceMyValues(c, nfaces, Acf_cells_[c].Values(), faces_LID);
    if (!symmetric()) {
      (*Afc_).ReplaceMyValues(c, nfaces, Afc_cells_[c].Values(), faces_LID);
    }
    (*Aff_).SumIntoGlobalValues(nfaces, faces_GID, Aff_cells_[c].values());

    // assemble rhs
    rhs_c[0][c] = Fc_cells_[c];
    for (int n=0; n!=nfaces; ++n) {
      rhs_f[0][faces[n]] += Ff_cells_[c][n];
    }
  }

  // communicate
  (*Aff_).GlobalAssemble();
  rhs_->GatherGhostedToMaster("face");

  // tag matrices as assembled
  assembled_operator_ = true;
}


/* ******************************************************************
 * Compute the face Schur complement of 2x2 block matrix.
 ****************************************************************** */
void MatrixMFD::ComputeSchurComplement(const std::vector<MatrixBC>& bc_markers,
        const std::vector<double>& bc_values) {
  // initialize to zero
  Sff_->PutScalar(0.0);

  // loop over cells and assemble
  AmanziMesh::Entity_ID_List faces_LID;
  std::vector<int> dirs;
  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);

  for (int c=0; c!=ncells; ++c) {
    mesh_->cell_get_faces_and_dirs(c, &faces_LID, &dirs);
    int nfaces = faces_LID.size();
    Epetra_SerialDenseMatrix Tff(nfaces, nfaces); // T implies local S
    Epetra_SerialDenseVector& Bcf = Acf_cells_[c];
    Epetra_SerialDenseVector& Bfc = Afc_cells_[c];

    for (int n=0; n!=nfaces; ++n) {
      for (int m=0; m!=nfaces; ++m) {
        Tff(n, m) = Aff_cells_[c](n, m) - Bfc[n] * Bcf[m] / (*Acc_)[c];
      }
    }

    for (int n=0; n!=nfaces; ++n) {  // boundary conditions
      int f=faces_LID[n];
      if (bc_markers[f] == MATRIX_BC_DIRICHLET) {
        for (int m=0; m!=nfaces; ++m) Tff(n, m) = Tff(m, n) = 0.0;
        Tff(n, n) = 1.0;
      }
    }

    Epetra_IntSerialDenseVector faces_GID(nfaces);
    for (int n=0; n!=nfaces; ++n) faces_GID[n] = (*Acf_).ColMap().GID(faces_LID[n]);
    (*Sff_).SumIntoGlobalValues(faces_GID, Tff);
  }
  (*Sff_).GlobalAssemble();

  // tag matrices as assembled
  assembled_schur_ = true;

}


/* ******************************************************************
 * Parallel matvec product Y <-- A * X.
 ****************************************************************** */
int MatrixMFD::Apply(const CompositeVector& X,
                      CompositeVector& Y) const {
  AssertAssembledOperator_or_die_();

  int ierr;

  // Face unknowns:  Yf = Aff_ * Xf + Afc_ * Xc
  ierr = (*Aff_).Multiply(false, *X.ViewComponent("face",false),
                          *Y.ViewComponent("face", false));

  Epetra_MultiVector Tf(*Y.ViewComponent("face", false));
  ierr |= (*Afc_).Multiply(true, *X.ViewComponent("cell",false), Tf);
  Y.ViewComponent("face",false)->Update(1.0, Tf, 1.0);

  // Cell unknowns:  Yc = Acf_ * Xf + Acc_ * Xc
  ierr |= (*Acf_).Multiply(false, *X.ViewComponent("face", false),
                           *Y.ViewComponent("cell", false));
  ierr |= Y.ViewComponent("cell", false)->Multiply(1.0, *Acc_,
          *X.ViewComponent("cell", false), 1.0);

  if (ierr) {
    Errors::Message msg("MatrixMFD::Apply has failed to calculate Y = inv(A) * X.");
    Exceptions::amanzi_throw(msg);
  }
  return ierr;
}


/* ******************************************************************
 * Parallel solve, Y <-- A^-1 X
 ****************************************************************** */
int MatrixMFD::ApplyInverse(const CompositeVector& X,
                             CompositeVector& Y) const {
  AssertAssembledOperator_or_die_();
  AssertAssembledSchur_or_die_();

  if (prec_method_ == PREC_METHOD_NULL) {
    Errors::Message msg("MatrixMFD::ApplyInverse requires a specified preconditioner method");
    Exceptions::amanzi_throw(msg);
  }

  // Temporary cell and face vectors.
  Epetra_MultiVector Tc(*Y.ViewComponent("cell", false));
  Epetra_MultiVector Tf(*Y.ViewComponent("face", false));

  // FORWARD ELIMINATION:  Tf = Xf - Afc_ inv(Acc_) Xc
  int ierr;
  ierr  = Tc.ReciprocalMultiply(1.0, *Acc_, *X.ViewComponent("cell", false), 0.0);
  ASSERT(!ierr);

  ierr |= (*Afc_).Multiply(true, Tc, Tf);  // Afc_ is kept in transpose form
  ASSERT(!ierr);

  Tf.Update(1.0, *X.ViewComponent("face", false), -1.0);

  // Solve the Schur complement system Sff_ * Yf = Tf.
  if (prec_method_ == TRILINOS_ML) {
    ierr |= ml_prec_->ApplyInverse(Tf, *Y.ViewComponent("face", false));
  } else if (prec_method_ == TRILINOS_ILU) {
    ierr |= ilu_prec_->ApplyInverse(Tf, *Y.ViewComponent("face", false));
  } else if (prec_method_ == TRILINOS_BLOCK_ILU) {
    ierr |= ifp_prec_->ApplyInverse(Tf, *Y.ViewComponent("face", false));
#ifdef HAVE_HYPRE
  } else if (prec_method_ == HYPRE_AMG || prec_method_ == HYPRE_EUCLID) {
    ierr |= IfpHypre_Sff_->ApplyInverse(Tf, *Y.ViewComponent("face", false));
#endif
  } else {
    ASSERT(0);
  }
  ASSERT(!ierr);

  // BACKWARD SUBSTITUTION:  Yc = inv(Acc_) (Xc - Acf_ Yf)
  ierr |= (*Acf_).Multiply(false, *Y.ViewComponent("face", false), Tc);  // It performs the required parallel communications.
  ASSERT(!ierr);

  Tc.Update(1.0, *X.ViewComponent("cell", false), -1.0);
  ierr |= Y.ViewComponent("cell", false)->ReciprocalMultiply(1.0, *Acc_, Tc, 0.0);

  if (ierr) {
    Errors::Message msg("MatrixMFD::ApplyInverse has failed in calculating y = A*x.");
    Exceptions::amanzi_throw(msg);
  }
  return ierr;
}


/* ******************************************************************
 * Linear algebra operations with matrices: r = f - A * x
 ****************************************************************** */
void MatrixMFD::ComputeResidual(const CompositeVector& solution,
        const Teuchos::Ptr<CompositeVector>& residual) const {
  Apply(solution, *residual);
  residual->Update(1.0, *rhs_, -1.0);
}


/* ******************************************************************
 * Linear algebra operations with matrices: r = A * x - f
 ****************************************************************** */
void MatrixMFD::ComputeNegativeResidual(const CompositeVector& solution,
        const Teuchos::Ptr<CompositeVector>& residual) const {
  Apply(solution, *residual);
  residual->Update(-1.0, *rhs_, 1.0);
}


/* ******************************************************************
 * Initialization of the preconditioner
 ****************************************************************** */
void MatrixMFD::InitPreconditioner() {
  if (prec_method_ == TRILINOS_ML) {
    ml_plist_ =  plist_.sublist("ML Parameters");
    ml_prec_ = Teuchos::rcp(new ML_Epetra::MultiLevelPreconditioner(*Sff_, ml_plist_, false));
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
    hypre_print_level_ = hypre_plist_.get<int>("print level",0);
    hypre_max_row_sum_ = hypre_plist_.get<double>("max row sum",0.9);
    hypre_max_levels_ = hypre_plist_.get<int>("max levels",25);
    hypre_relax_wt_ = hypre_plist_.get<double>("relax wt",1.0);
    hypre_interp_type_ = hypre_plist_.get<int>("interpolation type",0);
    hypre_agg_num_levels_ = hypre_plist_.get<int>("aggressive coarsening levels",0);
    hypre_agg_num_paths_ = hypre_plist_.get<int>("aggressive coarsening paths",1);

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
void MatrixMFD::UpdatePreconditioner() {
  if (prec_method_ == TRILINOS_ML) {
    if (ml_prec_->IsPreconditionerComputed()) ml_prec_->DestroyPreconditioner();
    ml_prec_->SetParameterList(ml_plist_);
    ml_prec_->ComputePreconditioner();
  } else if (prec_method_ == TRILINOS_ILU) {
    ilu_prec_ = Teuchos::rcp(new Ifpack_ILU(&*Sff_));
    ilu_prec_->SetParameters(ilu_plist_);
    ilu_prec_->Initialize();
    ilu_prec_->Compute();
  } else if (prec_method_ == TRILINOS_BLOCK_ILU) {
    Ifpack factory;
    std::string prectype("ILU");
    int ovl = ifp_plist_.get<int>("overlap",0);
    ifp_plist_.set<std::string>("schwarz: combine mode","Add");
    ifp_prec_ = Teuchos::rcp(factory.Create(prectype, &*Sff_, ovl));
    ifp_prec_->SetParameters(ifp_plist_);
    ifp_prec_->Initialize();
    ifp_prec_->Compute();
#ifdef HAVE_HYPRE
  } else if (prec_method_ == HYPRE_AMG) {
    IfpHypre_Sff_ = Teuchos::rcp(new Ifpack_Hypre(&*Sff_));
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
    IfpHypre_Sff_ = Teuchos::rcp(new Ifpack_Hypre(&*Sff_));

    Teuchos::ParameterList hypre_list;
    hypre_list.set("Preconditioner", Euclid);
    hypre_list.set("SolveOrPrecondition", Preconditioner);
    hypre_list.set("SetPreconditioner", true);
    hypre_list.set("NumFunctions", 0);

    IfpHypre_Sff_->SetParameters(hypre_list);
    IfpHypre_Sff_->Initialize();
    IfpHypre_Sff_->Compute();
  } else if (prec_method_ == HYPRE_PARASAILS) {
    IfpHypre_Sff_ = Teuchos::rcp(new Ifpack_Hypre(&*Sff_));

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


/* ******************************************************************
 * WARNING: Routines requires original mass matrices (Aff_cells_), i.e.
 * before boundary conditions were imposed.
 *
 * WARNING: Since diffusive flux is not continuous, we derive it only
 * once (using flag) and in exactly the same manner as in routine
 * Flow_PK::addGravityFluxes_DarcyFlux.
 *
 * WARNING: THIS ASSUMES solution has previously be communicated to update
 * ghost faces.
 ****************************************************************** */
void MatrixMFD::DeriveFlux(const CompositeVector& solution,
                           const Teuchos::Ptr<CompositeVector>& flux) const {

  AmanziMesh::Entity_ID_List faces;
  std::vector<double> dp;
  std::vector<int> dirs;

  flux->PutScalar(0.);

  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  int nfaces_owned = flux->size("face",false);
  solution.ScatterMasterToGhosted("face");

  std::vector<bool> done(nfaces_owned, false);
  const Epetra_MultiVector& soln_cells = *solution.ViewComponent("cell",false);
  const Epetra_MultiVector& soln_faces = *solution.ViewComponent("face",true);
  Epetra_MultiVector& flux_v = *flux->ViewComponent("face",false);

  for (int c=0; c!=ncells_owned; ++c) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    dp.resize(nfaces);
    for (int n=0; n!=nfaces; ++n) {
      int f = faces[n];
      dp[n] = soln_cells[0][c] - soln_faces[0][f];
    }

    for (int n=0; n!=nfaces; ++n) {
      int f = faces[n];
      if (f < nfaces_owned && !done[f]) {
        double s = 0.0;
        for (int m=0; m!=nfaces; ++m) {
          s += Aff_cells_[c](n, m) * dp[m];
        }

        flux_v[0][f] = s * dirs[n];
        done[f] = true;
      }
    }
  }

  // ensure post-condition - we got them all
  for (int f=0; f!=nfaces_owned; ++f) {
    ASSERT(done[f]);
  }

}


/* ******************************************************************
 * Derive Darcy velocity in cells.
 * WARNING: It cannot be consistent with the Darcy flux.
 * WARNING: It is assumed that flux faces have been communicated.
 ****************************************************************** */
void MatrixMFD::DeriveCellVelocity(const CompositeVector& flux,
        const Teuchos::Ptr<CompositeVector>& velocity) const {

  flux.ScatterMasterToGhosted("face");
  const Epetra_MultiVector& flux_f = *flux.ViewComponent("face",true);
  Epetra_MultiVector& velocity_c = *velocity->ViewComponent("cell",false);

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
        rhs_cell[i] += normal[i] * flux_f[0][f];
        matrix(i,i) += normal[i] * normal[i];
        for (int j=i+1; j!=dim; ++j) {
          matrix(j,i) = matrix(i,j) += normal[i] * normal[j];
        }
      }
    }

    int info;
    lapack.POSV('U', dim, 1, matrix.values(), dim, rhs_cell, dim, &info);

    for (int i=0; i!=dim; ++i) velocity_c[i][c] = rhs_cell[i];
  }
}


/* ******************************************************************
 * Solve the bottom row of the block system for lambda, given p.
 ****************************************************************** */
void MatrixMFD::UpdateConsistentFaceConstraints(const Teuchos::Ptr<CompositeVector>& u) {
  AssertAssembledOperator_or_die_();

  Teuchos::RCP<Epetra_MultiVector> uc = u->ViewComponent("cell", false);
  Teuchos::RCP<Epetra_MultiVector> rhs_f = rhs_->ViewComponent("face", false);
  Teuchos::RCP<Epetra_MultiVector> update_f =
      Teuchos::rcp(new Epetra_MultiVector(*rhs_f));

  Afc_->Multiply(true, *uc, *update_f);  // Afc is kept in the transpose form.
  rhs_f->Update(-1.0, *update_f, 1.0);

  // Replace the schur complement so it can be used as a face-only system
  Teuchos::RCP<Epetra_FECrsMatrix> Sff_tmp = Sff_;
  Sff_ = Aff_;

  // Update the preconditioner with a solver
  UpdatePreconditioner();

  // Use this entry to get appropriate faces.
  int ierr;
  if (prec_method_ == TRILINOS_ML) {
    ierr = ml_prec_->ApplyInverse(*rhs_f, *u->ViewComponent("face",false));
  } else if (prec_method_ == TRILINOS_ILU) {
    ierr = ilu_prec_->ApplyInverse(*rhs_f, *u->ViewComponent("face",false));
  } else if (prec_method_ == TRILINOS_BLOCK_ILU) {
    ierr = ifp_prec_->ApplyInverse(*rhs_f, *u->ViewComponent("face",false));
#ifdef HAVE_HYPRE
  } else if (prec_method_ == HYPRE_AMG || prec_method_ == HYPRE_EUCLID) {
    ierr = IfpHypre_Sff_->ApplyInverse(*rhs_f, *u->ViewComponent("face",false));
#endif
  } else {
    ASSERT(0);
  }

  ASSERT(!ierr);

  // revert schur complement
  Sff_ = Sff_tmp;
}


/* ******************************************************************
 * Solve the bottom row of the block system for lambda, given p.
 ****************************************************************** */
void MatrixMFD::UpdateConsistentFaceCorrection(const CompositeVector& u,
        const Teuchos::Ptr<CompositeVector>& Pu) {
  AssertAssembledOperator_or_die_();

  Teuchos::RCP<const Epetra_MultiVector> Pu_c = Pu->ViewComponent("cell", false);
  Epetra_MultiVector& Pu_f = *Pu->ViewComponent("face", false);
  const Epetra_MultiVector& u_f = *u.ViewComponent("face", false);
  Epetra_MultiVector update_f(u_f);

  Afc_->Multiply(true, *Pu_c, update_f);  // Afc is kept in the transpose form.
  update_f.Update(1., u_f, -1.);

  // Replace the schur complement so it can be used as a face-only system
  Teuchos::RCP<Epetra_FECrsMatrix> Sff_tmp = Sff_;
  Sff_ = Aff_;

  // Update the preconditioner with a solver
  UpdatePreconditioner();

  // Use this entry to get appropriate faces.
  int ierr;
  if (prec_method_ == TRILINOS_ML) {
    ierr = ml_prec_->ApplyInverse(update_f, Pu_f);
  } else if (prec_method_ == TRILINOS_ILU) {
    ierr = ilu_prec_->ApplyInverse(update_f, Pu_f);
  } else if (prec_method_ == TRILINOS_BLOCK_ILU) {
    ierr = ifp_prec_->ApplyInverse(update_f, Pu_f);
#ifdef HAVE_HYPRE
  } else if (prec_method_ == HYPRE_AMG || prec_method_ == HYPRE_EUCLID) {
    ierr = IfpHypre_Sff_->ApplyInverse(update_f, Pu_f);
#endif
  }
  ASSERT(!ierr);

  // revert schur complement
  Sff_ = Sff_tmp;
}


/* ******************************************************************
 * Solve the top row of the block system for p, given lambda.
 ****************************************************************** */
void MatrixMFD::UpdateConsistentCellCorrection(const CompositeVector& u,
        const Teuchos::Ptr<CompositeVector>& Pu) {
  AssertAssembledOperator_or_die_();

  int ierr(0);
  Epetra_MultiVector Tc(*Pu->ViewComponent("cell", false));

  // BACKWARD SUBSTITUTION:  Yc = inv(Acc_) (Xc - Acf_ Yf)
  ierr |= (*Acf_).Multiply(false, *Pu->ViewComponent("face", false), Tc);
  Tc.Update(1.0, *u.ViewComponent("cell", false), -1.0);

  ierr |= Pu->ViewComponent("cell", false)->ReciprocalMultiply(1.0, *Acc_, Tc, 0.0);

  if (ierr) {
    Errors::Message msg("MatrixMFD::ApplyInverse has failed in calculating y = A*x.");
    Exceptions::amanzi_throw(msg);
  }
}


/* ******************************************************************
 * Checks that the global matrices have been assembled.
 ****************************************************************** */
void MatrixMFD::AssertAssembledOperator_or_die_() const {
  if (!assembled_operator_) {
    Errors::Message msg("MatrixMFD: Operator has not been assembled.");
    Exceptions::amanzi_throw(msg);
  }
}


/* ******************************************************************
 * Checks that the Schur complement have been assembled.
 ****************************************************************** */
void MatrixMFD::AssertAssembledSchur_or_die_() const {
  if (!assembled_schur_) {
    Errors::Message msg("MatrixMFD: Schur complement has not been assembled.");
    Exceptions::amanzi_throw(msg);
  }
}

}  // namespace Operators
}  // namespace Amanzi
