/*
  This is the flow component of the Amanzi code.
  License: BSD
  Authors: Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
*/

#include "Teuchos_LAPACK.hpp"
#include "Epetra_FECrsGraph.h"
#include "EpetraExt_RowMatrixOut.h"

#include "errors.hh"
#include "MatrixMFD.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
 * Constructor
 ****************************************************************** */
MatrixMFD::MatrixMFD(Teuchos::ParameterList& plist,
                     const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) :
    plist_(plist),
    mesh_(mesh),
    flag_symmetry_(false),
    assembled_operator_(false),
    assembled_schur_(false),
    method_(MFD3D_NULL) 
{
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
    method_(other.method_)
{
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

  // vector space
  std::vector<std::string> names;
  names.push_back("cell"); names.push_back("face");
  std::vector<AmanziMesh::Entity_kind> locations;
  locations.push_back(AmanziMesh::CELL); locations.push_back(AmanziMesh::FACE);
  std::vector<int> ndofs(2,1);
  space_ = Teuchos::rcp(new CompositeVectorSpace());
  space_->SetMesh(mesh_)->SetComponents(names,locations,ndofs);

  // preconditioner
  if (plist_.isSublist("preconditioner")) {
    Teuchos::ParameterList pc_list = plist_.sublist("preconditioner");
    AmanziPreconditioners::PreconditionerFactory pc_fac;
    S_pc_ = pc_fac.Create(pc_list);
  }

  // verbose object
  vo_ = Teuchos::rcp(new VerboseObject("MatrixMFD", plist_));

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
  //  AmanziMesh::Entity_ID_List faces;

  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);

  if (Mff_cells_.size() != ncells) {
   Mff_cells_.resize(static_cast<size_t>(ncells));
  }

  int ok;
  nokay_ = npassed_ = 0;

  WhetStone::Tensor Kc;
  if (K == Teuchos::null) {
    Kc.init(mesh_->space_dimension(), 1);
    Kc(0,0) = 1.0;
  }

  //int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c=0; c!=ncells; ++c) {
    int nfaces = mesh_->cell_get_num_faces(c);

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
    
    Mff_cells_[c] = Mff;
    
    if (ok == WhetStone::WHETSTONE_ELEMENTAL_MATRIX_FAILED) {
      Errors::Message msg("Matrix_MFD: unexpected failure of LAPACK in WhetStone.");
      Exceptions::amanzi_throw(msg);
    }

    if (ok == WhetStone::WHETSTONE_ELEMENTAL_MATRIX_OK) nokay_++;
    if (ok == WhetStone::WHETSTONE_ELEMENTAL_MATRIX_PASSED) npassed_++;
  }

  // sum up the numbers across processors
  int nokay_tmp = nokay_, npassed_tmp = npassed_;
  Comm().SumAll(&nokay_tmp, &nokay_, 1);
  Comm().SumAll(&npassed_tmp, &npassed_, 1);
}


/* ******************************************************************
 * Calculate elemental stiffness matrices.
 ****************************************************************** */
void MatrixMFD::CreateMFDstiffnessMatrices(
    const Teuchos::Ptr<const CompositeVector>& Krel) {
  // tag global matrices as invalid
  assembled_schur_ = false;
  assembled_operator_ = false;

  // communicate as necessary
  if (Krel.get() && Krel->HasComponent("face"))
    Krel->ScatterMasterToGhosted("face");

  int dim = mesh_->space_dimension();
  WhetStone::MFD3D_Diffusion mfd(mesh_);
  AmanziMesh::Entity_ID_List faces;

  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  
  if (Aff_cells_.size() != ncells) {
    Aff_cells_.resize(static_cast<size_t>(ncells));
  }
  if (Afc_cells_.size() != ncells) {
    Afc_cells_.resize(static_cast<size_t>(ncells));
  }
  if (Acf_cells_.size() != ncells) {
    Acf_cells_.resize(static_cast<size_t>(ncells));
  }
  if (Acc_cells_.size() != ncells) {
    Acc_cells_.resize(static_cast<size_t>(ncells));
  }

  for (int c=0; c!=ncells; ++c) {
    int nfaces = mesh_->cell_get_num_faces(c);

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
      const Epetra_MultiVector& Krel_f = *Krel->ViewComponent("face",true);

      mesh_->cell_get_faces(c, &faces);

      for (int m=0; m!=nfaces; ++m) {
        AmanziMesh::Entity_ID f = faces[m];
        for (int n=0; n!=nfaces; ++n) {
          Bff(m, n) = Mff(m,n) * Krel_f[0][f];
        }
      }
    } else if (Krel->HasComponent("cell") && Krel->HasComponent("face")) {
      const Epetra_MultiVector& Krel_f = *Krel->ViewComponent("face",true);
      const Epetra_MultiVector& Krel_c = *Krel->ViewComponent("cell",false);

      mesh_->cell_get_faces(c, &faces);

      for (int m=0; m!=nfaces; ++m) {
        AmanziMesh::Entity_ID f = faces[m];
        for (int n=0; n!=nfaces; ++n) {
          Bff(m, n) = Mff(m,n) * Krel_c[0][c] * Krel_f[0][f];
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
    
    Aff_cells_[c] = Bff;
    Afc_cells_[c] = Bfc;
    Acf_cells_[c] = Bcf;
    Acc_cells_[c] = matsum;
  }
}


/* ******************************************************************
 * Create elemental rhs vectors.
 ****************************************************************** */
void MatrixMFD::CreateMFDrhsVectors() {
  
  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  
  if (Ff_cells_.size() != ncells) {
    Ff_cells_.resize(static_cast<size_t>(ncells));
  }
  if (Fc_cells_.size() != ncells) {
    Fc_cells_.resize(static_cast<size_t>(ncells));
  }

  for (int c=0; c!=ncells; ++c) {
    int nfaces = mesh_->cell_get_num_faces(c);

    Epetra_SerialDenseVector Ff(nfaces);  // Entries are initilaized to 0.0.
    double Fc = 0.0;

    Ff_cells_[c] = Ff;
    Fc_cells_[c] = Fc;
  }
}


/* ******************************************************************
 * Applies boundary conditions to elemental stiffness matrices and
 * adds contributions to elemental rigth-hand-sides.
 ****************************************************************** */
void MatrixMFD::ApplyBoundaryConditions(const std::vector<MatrixBC>& bc_markers,
					const std::vector<double>& bc_values, bool ADD_BC_FLUX) {
  // tag global matrices as invalid
  assembled_schur_ = false;
  assembled_operator_ = false;

  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  int nfaces = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  AmanziMesh::Entity_ID_List faces;
  AmanziMesh::Entity_ID_List cells;

  for (int c=0; c!=ncells; ++c) {
    mesh_->cell_get_faces(c, &faces);
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
      } else if ((bc_markers[f] == MATRIX_BC_FLUX)&&(ADD_BC_FLUX)) {
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
  int faces_LID[MFD_MAX_FACES];  // Contigious memory is required.
  int faces_GID[MFD_MAX_FACES];

  // fill the graphs
  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c=0; c!=ncells; ++c) {
    mesh_->cell_get_faces(c, &faces);
    int nfaces = faces.size();

    for (int n=0; n!=nfaces; ++n) {
      //     faces_LID[n] = faces[n];
      //      faces_GID[n] = fmap_wghost.GID(faces_LID[n]);
      faces_GID[n] = fmap_wghost.GID(faces[n]);
    }
    //    cf_graph->InsertMyIndices(c, nfaces, faces_LID);
    cf_graph->InsertMyIndices(c, nfaces, &(faces[0]));
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
    mesh_->cell_get_faces(c, &faces);
    int nfaces = faces.size();

    // assemble rhs (and simultaneously get GIDs of faces
    rhs_c[0][c] = Fc_cells_[c];
    for (int n=0; n!=nfaces; ++n) {
      AmanziMesh::Entity_ID f = faces[n];

      rhs_f[0][f] += Ff_cells_[c][n];

      faces_GID[n] = fmap_wghost.GID(f);
    }

    //    for (int n=0; n!=nfaces; ++n) {
      //      faces_LID[n] = faces[n];
      //      faces_GID[n] = fmap_wghost.GID(faces_LID[n]);
    //    }

    // assemble matrices
    (*Acc_)[c] = Acc_cells_[c];
    (*Acf_).ReplaceMyValues(c, nfaces, Acf_cells_[c].Values(), &(faces[0]));
    if (!symmetric()) {
      (*Afc_).ReplaceMyValues(c, nfaces, Afc_cells_[c].Values(), &(faces[0]));
    }
    (*Aff_).SumIntoGlobalValues(nfaces, faces_GID, Aff_cells_[c].values());

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
  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);

  for (int c=0; c!=ncells; ++c) {
    mesh_->cell_get_faces(c, &faces_LID);
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

  if (S_pc_ == Teuchos::null) {
    Errors::Message msg("MatrixMFD::ApplyInverse called but no preconditioner sublist was provided");
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
  ierr = S_pc_->ApplyInverse(Tf, *Y.ViewComponent("face",false));
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
void MatrixMFD::InitPreconditioner() {}

/* ******************************************************************
 * Rebuild preconditioner.
 ****************************************************************** */
void MatrixMFD::UpdatePreconditioner() {
  if (S_pc_ == Teuchos::null) {
    Errors::Message msg("MatrixMFD::UpdatePreconditioner called but no preconditioner sublist was provided");
    Exceptions::amanzi_throw(msg);
  }
  S_pc_->Destroy();
  S_pc_->Update(Sff_);
}


/* ******************************************************************
 * WARNING: Routines requires original mass matrices (Aff_cells_), i.e.
 * before boundary conditions were imposed.
 *
 * WARNING: Since diffusive flux is not continuous, we derive it only
 * once (using flag) and in exactly the same manner as in routine
 * Flow_PK::addGravityFluxes_DarcyFlux.
 *
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

  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c=0; c!=ncells_owned; ++c) {
    mesh_->cell_get_faces(c, &faces);
    int nfaces = faces.size();

    for (int i=0; i!=dim; ++i) rhs_cell[i] = 0.0;
    matrix.putScalar(0.0);

    for (int n=0; n!=nfaces; ++n) {  // populate least-square matrix
      int f = faces[n];
      const AmanziGeometry::Point& normal = mesh_->face_normal(f);
      //      double area = mesh_->face_area(f);

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
  
  // Aff solutions
  if (Aff_solver_ == Teuchos::null) {
    if (plist_.isSublist("consistent face solver")) {
      Teuchos::ParameterList Aff_plist = plist_.sublist("consistent face solver");
      Aff_op_ = Teuchos::rcp(new EpetraMatrixDefault<Epetra_FECrsMatrix>(Aff_plist));
      Aff_op_->Update(Aff_);

      if (Aff_plist.isParameter("iterative method")) {
        AmanziSolvers::LinearOperatorFactory<EpetraMatrix,Epetra_Vector,Epetra_BlockMap> op_fac;
        Aff_solver_ = op_fac.Create(Aff_plist, Aff_op_);
      } else {
        Aff_solver_ = Aff_op_;
      }
    } else {
      Errors::Message msg("MatrixMFD::UpdateConsistentFaceConstraints was called, but no consistent face solver sublist was provided.");
      Exceptions::amanzi_throw(msg);
    }
  }

  Teuchos::RCP<Epetra_MultiVector> uc = u->ViewComponent("cell", false);
  Teuchos::RCP<Epetra_MultiVector> rhs_f = rhs_->ViewComponent("face", false);
  Teuchos::RCP<Epetra_MultiVector> update_f =
      Teuchos::rcp(new Epetra_MultiVector(*rhs_f));

  Afc_->Multiply(true, *uc, *update_f);  // Afc is kept in the transpose form.
  update_f->Update(1.0, *rhs_f, -1.0);

  Aff_op_->Destroy();
  Aff_op_->Update(Aff_);
  int ierr = Aff_solver_->ApplyInverse(*(*update_f)(0), *(*u->ViewComponent("face",false))(0));
  ASSERT(!ierr);
}


/* ******************************************************************
 * Solve the bottom row of the block system for lambda, given p.
 ****************************************************************** */
void MatrixMFD::UpdateConsistentFaceCorrection(const CompositeVector& u,
        const Teuchos::Ptr<CompositeVector>& Pu) {
  AssertAssembledOperator_or_die_();

  // Aff solutions
  if (Aff_solver_ == Teuchos::null) {
    if (plist_.isSublist("consistent face solver")) {
      Teuchos::ParameterList Aff_plist = plist_.sublist("consistent face solver");
      Aff_op_ = Teuchos::rcp(new EpetraMatrixDefault<Epetra_FECrsMatrix>(Aff_plist));
      Aff_op_->Update(Aff_);

      if (Aff_plist.isParameter("iterative method")) {
        AmanziSolvers::LinearOperatorFactory<EpetraMatrix,Epetra_Vector,Epetra_BlockMap> op_fac;
        Aff_solver_ = op_fac.Create(Aff_plist, Aff_op_);
      } else {
        Aff_solver_ = Aff_op_;
      }
    } else {
      Errors::Message msg("MatrixMFD::UpdateConsistentFaceConstraints was called, but no consistent face solver sublist was provided.");
      Exceptions::amanzi_throw(msg);
    }
  }

  Teuchos::RCP<const Epetra_MultiVector> Pu_c = Pu->ViewComponent("cell", false);
  Epetra_MultiVector& Pu_f = *Pu->ViewComponent("face", false);
  const Epetra_MultiVector& u_f = *u.ViewComponent("face", false);
  Epetra_MultiVector update_f(u_f);

  Afc_->Multiply(true, *Pu_c, update_f);  // Afc is kept in the transpose form.
  update_f.Update(1., u_f, -1.);

  Aff_op_->Destroy();
  Aff_op_->Update(Aff_);
  int ierr = Aff_solver_->ApplyInverse(*update_f(0), *Pu_f(0));
  ASSERT(!ierr);
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

void MatrixMFD::Add2MFDstiffnessMatrices(std::vector<double> *Acc_ptr,
			      std::vector<Teuchos::SerialDenseMatrix<int, double> > *Aff_ptr,
			      std::vector<Epetra_SerialDenseVector> *Acf_ptr,
			      std::vector<Epetra_SerialDenseVector> *Afc_ptr){

  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);

  if (Acc_ptr){
    for (int c=0; c!=ncells; ++c) {
      //cout<<Acc_cells_[c]<<" "<<(*Acc_ptr)[c]<<endl;
      Acc_cells_[c] += (*Acc_ptr)[c];
    }
  }
  if (Afc_ptr){
    for (int c=0; c!=ncells; ++c) {
      Afc_cells_[c] += Afc_ptr->at(c);
    }
  }
  if (Acf_ptr){
    for (int c=0; c!=ncells; ++c) {
      Acf_cells_[c] += Acf_ptr->at(c);
    }
  }
  if (Aff_ptr){
    for (int c=0; c!=ncells; ++c) {
      Aff_cells_[c] += Aff_ptr->at(c);
    }
  }

}

}  // namespace Operators
}  // namespace Amanzi
