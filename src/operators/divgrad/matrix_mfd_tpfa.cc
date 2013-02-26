/*
  This is the flow component of the Amanzi code.

  Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
  Daniil Svyatskiy (dasvyat@lanl.gov)
  Ethan Coon (ATS version) (ecoon@lanl.gov)
*/

#include <vector>

#include "Epetra_FECrsGraph.h"

#include "errors.hh"
#include "matrix_mfd_tpfa.hh"


namespace Amanzi {
namespace Operators {

template<class T>
int FindPosition(const std::vector<T>& v, const T& value) {
  for (int i = 0; i < v.size(); i++)
    if (v[i] == value) return i;
  return -1;
}



/* ******************************************************************
 * Calculate elemental stiffness matrices.
 ****************************************************************** */
void MatrixMFD_TPFA::CreateMFDstiffnessMatrices(
    const Teuchos::Ptr<const CompositeVector>& Krel) {

  int dim = mesh_->space_dimension();
  WhetStone::MFD3D mfd(mesh_);
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  Aff_cells_.clear();
  Afc_cells_.clear();
  Acf_cells_.clear();
  Acc_cells_.clear();

  Teuchos::RCP<const Epetra_MultiVector> Krel_cells;
  Teuchos::RCP<const Epetra_MultiVector> Krel_faces;
  if (Krel != Teuchos::null) {
    if (Krel->has_component("cell")) {
      Krel_cells = Krel->ViewComponent("cell",false);
    }
    if (Krel->has_component("face")) {
      Krel_faces = Krel->ViewComponent("face",true);
    }
  }

  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c = 0; c < ncells; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    Teuchos::SerialDenseMatrix<int, double>& Mff = Mff_cells_[c];
    Teuchos::SerialDenseMatrix<int, double> Bff(nfaces, nfaces);
    Epetra_SerialDenseVector Bcf(nfaces), Bfc(nfaces);

    if (Krel == Teuchos::null) {
      for (int n = 0; n < nfaces; n++) Bff(n, n) = Mff(n, n);
    } else if (Krel_faces == Teuchos::null) {
      for (int n = 0; n < nfaces; n++) Bff(n, n) = Mff(n, n) * (*Krel_cells)[0][c];
    } else if (Krel_cells == Teuchos::null) {
      for (int n = 0; n < nfaces; n++) Bff(n, n) = Mff(n, n) * (*Krel_faces)[0][faces[n]];
    } else {
      for (int n = 0; n < nfaces; n++) Bff(n, n) = Mff(n, n)
                                           * (*Krel_faces)[0][faces[n]] * (*Krel_cells)[0][c];
    }

    double matsum = 0.0;  // elimination of mass matrix
    for (int n = 0; n < nfaces; n++) {
      double rowsum = Bff(n, n), colsum = Bff(n, n);
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
 * Initialize Trilinos matrices. It must be called only once.
 * If matrix is non-symmetric, we generate transpose of the matrix
 * block Afc to reuse cf_graph; otherwise, pointer Afc = Acf.
 ****************************************************************** */
void MatrixMFD_TPFA::SymbolicAssembleGlobalMatrices() {
  MatrixMFD::SymbolicAssembleGlobalMatrices();

  const Epetra_Map& cmap = mesh_->cell_map(false);
  const Epetra_Map& cmap_wghost = mesh_->cell_map(true);
  const Epetra_Map& fmap_wghost = mesh_->face_map(true);

  int avg_entries_row = (mesh_->space_dimension() == 2) ? MFD_QUAD_FACES : MFD_HEX_FACES;
  Epetra_FECrsGraph pp_graph(Copy, cmap, avg_entries_row + 1);

  AmanziMesh::Entity_ID_List cells;
  int cells_GID[2];

  int nfaces = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  for (int f = 0; f < nfaces; f++) {
    mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
    int ncells = cells.size();

    for (int n = 0; n < ncells; n++)
      cells_GID[n] = cmap_wghost.GID(cells[n]);

    pp_graph.InsertGlobalIndices(ncells, cells_GID, ncells, cells_GID);
  }
  pp_graph.GlobalAssemble();  // Symbolic graph is complete.

  // create global matrices
  std::vector<std::string> names(1,"face");
  std::vector<AmanziMesh::Entity_kind> locations(1,AmanziMesh::FACE);
  std::vector<int> ndofs(1,1);
  Dff_ = Teuchos::rcp(new CompositeVector(mesh_, names, locations, ndofs, true));
  Dff_->CreateData();

  Spp_ = Teuchos::rcp(new Epetra_FECrsMatrix(Copy, pp_graph));
  Spp_->GlobalAssemble();
}


/* ******************************************************************
 * Convert elemental mass matrices into stiffness matrices and
 * assemble them into four global matrices.
 * We need an auxiliary GHOST-based vector to assemble the RHS.
 ****************************************************************** */
void MatrixMFD_TPFA::AssembleGlobalMatrices() {
  std::vector<std::string> names_c(1,"cell");
  std::vector<AmanziMesh::Entity_kind> locations_c(1,AmanziMesh::CELL);
  std::vector<std::string> names_f(1,"face");
  std::vector<AmanziMesh::Entity_kind> locations_f(1,AmanziMesh::FACE);
  std::vector<int> ndofs(1,1);

  MatrixMFD::AssembleGlobalMatrices();

  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;
  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);

  // Dff
  Dff_->PutScalar(0.0);
  Epetra_MultiVector& Dff_f = *Dff_->ViewComponent("face",true);
  for (int c = 0; c < ncells_owned; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int mfaces = faces.size();

    for (int n = 0; n < mfaces; n++) {
      int f = faces[n];
      Dff_f[0][f] += Aff_cells_[c](n, n);
    }
  }
  Dff_->GatherGhostedToMaster();

  // convert right-hand side to a cell-based vector
  const Epetra_Map& cmap = mesh_->cell_map(false);
  const Epetra_Map& cmap_wghost = mesh_->cell_map(true);
  Epetra_Vector Tc(cmap);
  int nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);

  Teuchos::RCP<Epetra_MultiVector> rhs_faces = rhs_->ViewComponent("face",false);
  Teuchos::RCP<Epetra_MultiVector> rhs_cells = rhs_->ViewComponent("cell",false);
  for (int f=0; f!=nfaces_owned; ++f) (*rhs_faces)[0][f] /= Dff_f[0][f];
  (*Acf_).Multiply(false, *rhs_faces, Tc);
  for (int c=0; c!=ncells_owned; ++c) (*rhs_cells)[0][c] -= Tc[c];
  for (int f=0; f!=nfaces_owned; ++f) (*rhs_faces)[0][f] *= Dff_f[0][f];

  // create a with-ghost copy of Acc
  CompositeVector Dcc(mesh_,names_c,locations_c,ndofs,true);
  Dcc.CreateData();
  Epetra_MultiVector& Dcc_c = *Dcc.ViewComponent("cell",true);

  for (int c=0; c!=ncells_owned; ++c) Dcc_c[0][c] = (*Acc_)[c];
  Dcc.ScatterMasterToGhosted();

  AmanziMesh::Entity_ID_List cells;
  int cells_GID[2];
  double Acf_copy[2];

  // create auxiliaty with-ghost copy of Acf_cells
  CompositeVector Acf_parallel(mesh_,names_f,locations_f,ndofs,true);
  Acf_parallel.CreateData();
  Acf_parallel.PutScalar(0.);
  Epetra_MultiVector& Acf_parallel_f = *Acf_parallel.ViewComponent("face",true);

  for (int c=0; c!=ncells_owned; ++c) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int mfaces = faces.size();

    for (int i=0; i!=mfaces; ++i) {
      int f = faces[i];
      if (f >= nfaces_owned) Acf_parallel_f[0][f] = Acf_cells_[c][i];
    }
  }
  Acf_parallel.GatherGhostedToMaster();

  // populate the global matrix
  Spp_->PutScalar(0.0);
  for (AmanziMesh::Entity_ID f=0; f!=nfaces_owned; ++f) {
    mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
    int mcells = cells.size();

    // populate face-based matrix.
    Teuchos::SerialDenseMatrix<int, double> Bpp(mcells, mcells);
    for (int n = 0; n < mcells; n++) {
      int c = cells[n];
      cells_GID[n] = cmap_wghost.GID(c);

      mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
      int i = FindPosition<AmanziMesh::Entity_ID>(faces, f);
      Bpp(n, n) = Dcc_c[0][c] / faces.size();
      if (c < ncells_owned) {
        int i = FindPosition<AmanziMesh::Entity_ID>(faces, f);
        Acf_copy[n] = Acf_cells_[c][i];
      } else {
        Acf_copy[n] = Acf_parallel_f[0][f];
      }
    }

    for (int n = 0; n < mcells; n++) {
      for (int m = n; m < mcells; m++) {
        Bpp(n, m) -= Acf_copy[n] * Acf_copy[m] / Dff_f[0][f];
        Bpp(m, n) = Bpp(n, m);
      }
    }

    (*Spp_).SumIntoGlobalValues(mcells, cells_GID, Bpp.values());
  }
  (*Spp_).GlobalAssemble();
}



/* ******************************************************************
 * Initialization of the preconditioner
 ****************************************************************** */
void MatrixMFD_TPFA::InitPreconditioner(Teuchos::ParameterList& prec_plist) {
  if (prec_method_ == TRILINOS_ML) {
    ml_plist_ =  prec_plist.sublist("ML Parameters");
    ml_prec_ = Teuchos::rcp(new ML_Epetra::MultiLevelPreconditioner(*Spp_, ml_plist_, false));
  } else if (prec_method_ == TRILINOS_ILU) {
    ilu_plist_ = prec_plist.sublist("ILU Parameters");
  } else if (prec_method_ == TRILINOS_BLOCK_ILU) {
    ifp_plist_ = prec_plist.sublist("Block ILU Parameters");
#ifdef HAVE_HYPRE
  } else if (prec_method_ == HYPRE_AMG) {
    // read some boomer amg parameters
    hypre_plist_ = prec_plist.sublist("HYPRE AMG Parameters");
    hypre_ncycles_ = hypre_plist_.get<int>("number of cycles",5);
    hypre_nsmooth_ = hypre_plist_.get<int>("number of smoothing iterations",3);
    hypre_tol_ = hypre_plist_.get<double>("tolerance",0.0);
    hypre_strong_threshold_ = hypre_plist_.get<double>("strong threshold",0.25);
  } else if (prec_method_ == HYPRE_EUCLID) {
    hypre_plist_ = prec_plist.sublist("HYPRE Euclid Parameters");
  } else if (prec_method_ == HYPRE_PARASAILS) {
    hypre_plist_ = prec_plist.sublist("HYPRE ParaSails Parameters");
#endif
  }
}


/* ******************************************************************
 * Rebuild the preconditioner.
 ****************************************************************** */
void MatrixMFD_TPFA::UpdatePreconditioner() {
  if (prec_method_ == TRILINOS_ML) {
    if (ml_prec_->IsPreconditionerComputed()) ml_prec_->DestroyPreconditioner();
    ml_prec_->SetParameterList(ml_plist_);
    ml_prec_->ComputePreconditioner();
  } else if (prec_method_ == TRILINOS_ILU) {
    ilu_prec_ = Teuchos::rcp(new Ifpack_ILU(&*Spp_));
    ilu_prec_->SetParameters(ilu_plist_);
    ilu_prec_->Initialize();
    ilu_prec_->Compute();
  } else if (prec_method_ == TRILINOS_BLOCK_ILU) {
    Ifpack factory;
    std::string prectype("ILU");
    int ovl = ifp_plist_.get<int>("overlap",0);
    ifp_plist_.set<std::string>("schwarz: combine mode","Add");
    ifp_prec_ = Teuchos::rcp(factory.Create(prectype, &*Spp_, ovl));
    ifp_prec_->SetParameters(ifp_plist_);
    ifp_prec_->Initialize();
    ifp_prec_->Compute();
#ifdef HAVE_HYPRE
  } else if (prec_method_ == HYPRE_AMG) {
    IfpHypre_Spp_ = Teuchos::rcp(new Ifpack_Hypre(&*Spp_));
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

    IfpHypre_Spp_->SetParameters(hypre_list);
    IfpHypre_Spp_->Initialize();
    IfpHypre_Spp_->Compute();
  } else if (prec_method_ == HYPRE_EUCLID) {
    IfpHypre_Spp_ = Teuchos::rcp(new Ifpack_Hypre(&*Spp_));

    Teuchos::ParameterList hypre_list;
    hypre_list.set("Preconditioner", Euclid);
    hypre_list.set("SolveOrPrecondition", Preconditioner);
    hypre_list.set("SetPreconditioner", true);
    hypre_list.set("NumFunctions", 0);

    IfpHypre_Spp_->SetParameters(hypre_list);
    IfpHypre_Spp_->Initialize();
    IfpHypre_Spp_->Compute();
  } else if (prec_method_ == HYPRE_PARASAILS) {
    IfpHypre_Spp_ = Teuchos::rcp(new Ifpack_Hypre(&*Spp_));

    Teuchos::ParameterList hypre_list;
    hypre_list.set("Preconditioner", ParaSails);
    hypre_list.set("SolveOrPrecondition", Preconditioner);
    hypre_list.set("SetPreconditioner", true);
    hypre_list.set("NumFunctions", 0);

    IfpHypre_Spp_->SetParameters(hypre_list);
    IfpHypre_Spp_->Initialize();
    IfpHypre_Spp_->Compute();
#endif
  }
}



/* ******************************************************************
 * Parallel matvec product Spp * Xc.
 ****************************************************************** */
void MatrixMFD_TPFA::Apply(const CompositeVector& X,
                      const Teuchos::Ptr<CompositeVector>& Y) const {

  int ierr = (*Spp_).Multiply(false, *X.ViewComponent("cell",false),
          *Y->ViewComponent("cell",false));

  if (ierr) {
    Errors::Message msg("MatrixMFD_TPFA::Apply has failed to calculate y = A*x.");
    Exceptions::amanzi_throw(msg);
  }

  // Yf = Xf;
  Y->ViewComponent("face",false)->PutScalar(0.);
}



void MatrixMFD_TPFA::ApplyInverse(const CompositeVector& X,
        const Teuchos::Ptr<CompositeVector>& Y) const {
  // Solve the Schur complement system Spp * Yc = Xc. Since AztecOO may
  // use the same memory for X and Y, we introduce auxiliaty vector Tc.
  int ierr = 0;
  const Epetra_MultiVector& Xc = *X.ViewComponent("cell",false);
  Epetra_MultiVector Tc(Xc);

  // Solve the pp system
  if (prec_method_ == TRILINOS_ML) {
    ierr |= ml_prec_->ApplyInverse(Xc, Tc);
  } else if (prec_method_ == TRILINOS_ILU) {
    ierr |= ilu_prec_->ApplyInverse(Xc, Tc);
  } else if (prec_method_ == TRILINOS_BLOCK_ILU) {
    ierr |= ifp_prec_->ApplyInverse(Xc, Tc);
#ifdef HAVE_HYPRE
  } else if (prec_method_ == HYPRE_AMG || prec_method_ == HYPRE_EUCLID) {
    ierr != IfpHypre_Sff_->ApplyInverse(Xc, Tc);
#endif
  }
  *Y->ViewComponent("cell",false) = Tc;

  if (ierr) {
    Errors::Message msg("MatrixMFD_TPFA::ApplyInverse has failed in calculating y = inv(A)*x.");
    Exceptions::amanzi_throw(msg);
  }

  // copy over face values
  *Y->ViewComponent("face",false) = *X.ViewComponent("face",false);

}


void MatrixMFD_TPFA::UpdateConsistentFaceConstraints(
    const Teuchos::Ptr<CompositeVector>& u) {

  Teuchos::RCP<const Epetra_MultiVector> uc = u->ViewComponent("cell", false);
  Epetra_MultiVector& uf = *u->ViewComponent("face", false);

  Epetra_MultiVector& Dff_f = *Dff_->ViewComponent("face",false);
  Epetra_MultiVector& rhs_f = *rhs_->ViewComponent("face", false);
  Epetra_MultiVector update_f(rhs_f);


  Afc_->Multiply(true,*uc, update_f);  // Afc is kept in the transpose form.
  rhs_f.Update(-1.0, update_f, 1.0);

  int nfaces = rhs_f.MyLength();
  for (int f=0; f!=nfaces; ++f) {
    uf[0][f] = rhs_f[0][f] / Dff_f[0][f];
  }
}

void MatrixMFD_TPFA::UpdateConsistentFaceCorrection(const CompositeVector& u,
    const Teuchos::Ptr<CompositeVector>& Pu) {

  Teuchos::RCP<const Epetra_MultiVector> Pu_c = Pu->ViewComponent("cell", false);
  Epetra_MultiVector& Pu_f = *Pu->ViewComponent("face", false);
  const Epetra_MultiVector& u_f = *u.ViewComponent("face", false);

  Epetra_MultiVector& Dff_f = *Dff_->ViewComponent("face",false);

  Afc_->Multiply(true, *Pu_c, Pu_f);  // Afc is kept in the transpose form.
  Pu_f.Update(1., u_f, -1.);

  int nfaces = Pu_f.MyLength();
  for (int f=0; f!=nfaces; ++f) {
    Pu_f[0][f] /= Dff_f[0][f];
  }
}

}  // namespace AmanziFlow
}  // namespace Amanzi
