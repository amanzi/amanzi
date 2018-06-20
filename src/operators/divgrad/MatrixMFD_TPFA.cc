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
#include "EpetraExt_RowMatrixOut.h"

#include "errors.hh"
#include "MatrixMFD_TPFA.hh"


namespace Amanzi {
namespace Operators {

template<class T>
int FindPosition(const std::vector<T>& v, const T& value) {
  for (unsigned int i=0; i!=v.size(); ++i)
    if (v[i] == value) return i;
  return -1;
}


/* ******************************************************************
 * Calculate elemental stiffness matrices.
 ****************************************************************** */
void MatrixMFD_TPFA::CreateMFDstiffnessMatrices(
    const Teuchos::Ptr<const CompositeVector>& Krel) {
  // tag global matrices as invalid
  MarkLocalMatricesAsChanged_();

  // communicate as necessary
  if (Krel.get() && Krel->HasComponent("face"))
    Krel->ScatterMasterToGhosted("face");

  int dim = mesh_->space_dimension();
  WhetStone::MFD3D_Diffusion mfd(mesh_);
  AmanziMesh::Entity_ID_List faces;

  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  
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
    Acc_ = Teuchos::rcp(new Epetra_Vector(View,mesh_->cell_map(false),&Acc_cells_[0]));
  }  

  for (int c=0; c!=ncells; ++c) {
    mesh_->cell_get_faces(c, &faces);
    int nfaces = faces.size();

    WhetStone::DenseMatrix& Mff = Mff_cells_[c];
    Teuchos::SerialDenseMatrix<int, double> Bff(nfaces,nfaces);
    Epetra_SerialDenseVector Bcf(nfaces), Bfc(nfaces);

    if (Krel == Teuchos::null ||
        (!Krel->HasComponent("cell") && !Krel->HasComponent("face"))) {
      for (int n=0; n!=nfaces; ++n) {
        Bff(n, n) = Mff(n, n);
      }
    } else if (Krel->HasComponent("cell") && !Krel->HasComponent("face")) {
      const Epetra_MultiVector& Krel_c = *Krel->ViewComponent("cell",false);

      for (int n=0; n!=nfaces; ++n) {
        Bff(n, n) = Mff(n,n) * Krel_c[0][c];
      }
    } else if (!Krel->HasComponent("cell") && Krel->HasComponent("face")) {
      const Epetra_MultiVector& Krel_f = *Krel->ViewComponent("face",true);

      for (int n=0; n!=nfaces; ++n) {
        Bff(n, n) = Mff(n,n) * Krel_f[0][faces[n]];
      }
    } else if (Krel->HasComponent("cell") && Krel->HasComponent("face")) {
      const Epetra_MultiVector& Krel_c = *Krel->ViewComponent("cell",false);
      const Epetra_MultiVector& Krel_f = *Krel->ViewComponent("face",true);

      for (int n=0; n!=nfaces; ++n) {
        Bff(n, n) = Mff(n,n) * Krel_c[0][c] * Krel_f[0][faces[n]];
      }
    }

    double matsum = 0.0;  // elimination of mass matrix
    for (int n=0; n!=nfaces; ++n) {
      double rowsum = Bff(n, n);
      double colsum = Bff(n, n);

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
 * Initialize Trilinos matrices. It must be called only once.
 * If matrix is non-symmetric, we generate transpose of the matrix
 * block Afc to reuse cf_graph; otherwise, pointer Afc = Acf.
 ****************************************************************** */
void MatrixMFD_TPFA::SymbolicAssembleGlobalMatrices() {
  // get the standard matrices
  MatrixMFD::SymbolicAssembleGlobalMatrices();

  // also create a cell-cell matrix for the TPFA
  const Epetra_Map& cmap = mesh_->cell_map(false);
  const Epetra_Map& cmap_wghost = mesh_->cell_map(true);
  const Epetra_Map& fmap_wghost = mesh_->face_map(true);

  int avg_entries_row = (mesh_->space_dimension() == 2) ? MFD_QUAD_FACES : MFD_HEX_FACES;
  Epetra_FECrsGraph pp_graph(Copy, cmap, avg_entries_row + 1);

  AmanziMesh::Entity_ID_List cells;
  int cells_GID[2];

  int nfaces = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  for (int f = 0; f < nfaces; f++) {
    mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
    int ncells = cells.size();

    for (int n = 0; n < ncells; n++) cells_GID[n] = cmap_wghost.GID(cells[n]);

    pp_graph.InsertGlobalIndices(ncells, cells_GID, ncells, cells_GID);
  }
  pp_graph.GlobalAssemble();  // Symbolic graph is complete.

  // create global matrices
  std::vector<std::string> names(1,"face");
  std::vector<AmanziMesh::Entity_kind> locations(1,AmanziMesh::FACE);
  std::vector<int> ndofs(1,1);
  
  CompositeVectorSpace space;
  space.SetMesh(mesh_)->SetGhosted()->SetComponents(names,locations,ndofs);
  Dff_ = Teuchos::rcp(new CompositeVector(space));
  Dff_->PutScalar(0.);

  App_ = Teuchos::rcp(new Epetra_FECrsMatrix(Copy, pp_graph));
  App_->GlobalAssemble();
}

/* ******************************************************************
 * Assemble Dff, the diagonal Aff
 ****************************************************************** */
void MatrixMFD_TPFA::AssembleDff_() const {
  // Dff
  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  AmanziMesh::Entity_ID_List faces;

  int c=0;
  // std::cout<<"Aff\n"<<Aff_cells_[c]<<"\n";
  // std::cout<<"Acc\n"<<Acc_cells_[c]<<"\n";
  // std::cout<<"Afc\n"<<Acf_cells_[c]<<"\n";
  // Teuchos::RCP<Epetra_MultiVector> rhs_cells_tmp = rhs_->ViewComponent("cell",false);
  // std::cout<<"rhs_cell\n"<<*rhs_cells_tmp<<"\n";
  // Teuchos::RCP<Epetra_MultiVector> rhs_faces_tmp = rhs_->ViewComponent("face",false);

  //exit(0);

  Epetra_MultiVector& Dff_f = *Dff_->ViewComponent("face",true);
  Dff_f.PutScalar(0.);
  for (int c=0; c!=ncells_owned; ++c) {
    mesh_->cell_get_faces(c, &faces);

    int mfaces = faces.size();
    for (int n=0; n!=mfaces; ++n) {
      int f = faces[n];
      Dff_f[0][f] += Aff_cells_[c](n, n);
    }
  }
  Dff_->GatherGhostedToMaster();

  assembled_dff_ = true;
}




/* ******************************************************************
 * Assemble elemental rhs matrices into global RHS
 ****************************************************************** */
void MatrixMFD_TPFA::AssembleRHS_() const {
  MatrixMFD::AssembleRHS_();
  if (!assembled_dff_) AssembleDff_();
  const Epetra_MultiVector& Dff_f = *Dff_->ViewComponent("face",true);

  int nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  Teuchos::RCP<Epetra_MultiVector> rhs_faces = rhs_->ViewComponent("face",false);
  Teuchos::RCP<Epetra_MultiVector> rhs_cells = rhs_->ViewComponent("cell",false);

  // Schur complement RHS
  for (int f=0; f!=nfaces_owned; ++f) {
    if (std::abs(Dff_f[0][f]) > 0) {
      (*rhs_faces)[0][f] /= Dff_f[0][f];
    } else {
      AMANZI_ASSERT( (*rhs_faces)[0][f] == 0. );
    }
  }
  ApplyAcf(*rhs_, *rhs_, -1.);
  rhs_cells->Scale(-1.);

  // std::cout<<"rhs_faces\n"<<*rhs_faces<<"\n";
  // std::cout<<"Tc\n"<<Tc<<"\n";
  // std::cout<<"rhs_cell\n"<<*rhs_cells_tmp<<"\n";
  // exit(0);


  // unscale the rhs
  for (int f=0; f!=nfaces_owned; ++f) {
    (*rhs_faces)[0][f] *= Dff_f[0][f];
  }

  assembled_rhs_ = true;
}


/* ******************************************************************
 * App is also Schur complement
 ****************************************************************** */
void MatrixMFD_TPFA::AssembleSchur_() const {
  if (!assembled_app_) AssembleApp_();
  Sff_ = App_;
  assembled_schur_ = true;
}


/* ******************************************************************
 * Convert elemental mass matrices into stiffness matrices and
 * assemble them into global matrices.
 ****************************************************************** */
void MatrixMFD_TPFA::AssembleApp_() const {
  if (!assembled_dff_) AssembleDff_();

  std::vector<std::string> names_c(1,"cell");
  std::vector<AmanziMesh::Entity_kind> locations_c(1,AmanziMesh::CELL);
  std::vector<std::string> names_f(1,"face");
  std::vector<AmanziMesh::Entity_kind> locations_f(1,AmanziMesh::FACE);
  std::vector<int> ndofs(1,1);

  AmanziMesh::Entity_ID_List faces;
  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

  // convert right-hand side to a cell-based vector
  const Epetra_Map& fmap_wghost = mesh_->face_map(true);
  const Epetra_Map& cmap = mesh_->cell_map(false);
  const Epetra_Map& cmap_wghost = mesh_->cell_map(true);
  Epetra_Vector Tc(cmap);
  int nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);

  // create a with-ghost copy of Acc
  CompositeVectorSpace space_c;
  space_c.SetMesh(mesh_)->SetGhosted()->SetComponents(names_c,locations_c,ndofs);
  CompositeVector Dcc(space_c);

  { // context for non-const access
    Epetra_MultiVector& Dcc_c = *Dcc.ViewComponent("cell",true);
    for (int c=0; c!=ncells_owned; ++c) Dcc_c[0][c] = (*Acc_)[c];
  }
  Dcc.ScatterMasterToGhosted();
  Epetra_MultiVector& Dcc_c = *Dcc.ViewComponent("cell",true);

  AmanziMesh::Entity_ID_List cells;
  int cells_GID[2];
  double Acf_copy[2];
  double Afc_copy[2];

  // create auxiliaty with-ghost copy of Acf_cells
  CompositeVectorSpace space_f;
  space_f.SetMesh(mesh_)->SetGhosted()->SetComponents(names_f,locations_f,ndofs);
  CompositeVector Acf_parallel(space_f);
  CompositeVector Afc_parallel(space_f);

  Epetra_MultiVector& Acf_parallel_f = *Acf_parallel.ViewComponent("face",true);
  Epetra_MultiVector& Afc_parallel_f = *Afc_parallel.ViewComponent("face",true);
  // note the PutScalar is called on MultiVector to ensure ghost entries are
  // initialized
  Acf_parallel_f.PutScalar(0.);
  Afc_parallel_f.PutScalar(0.);

  for (int c=0; c!=ncells_owned; ++c) {
    mesh_->cell_get_faces(c, &faces);

    int mfaces = faces.size();
    for (int i=0; i!=mfaces; ++i) {
      int f = faces[i];
      if (f >= nfaces_owned) Acf_parallel_f[0][f] = Acf_cells_[c][i];
      if (f >= nfaces_owned) Afc_parallel_f[0][f] = Afc_cells_[c][i];
    }
  }

  Acf_parallel.GatherGhostedToMaster();
  Afc_parallel.GatherGhostedToMaster();

  // populate the global matrix
  const Epetra_MultiVector& Dff_f = *Dff_->ViewComponent("face",true);

  App_->PutScalar(0.0);
  for (AmanziMesh::Entity_ID f=0; f!=nfaces_owned; ++f) {
    mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
    int mcells = cells.size();

    // populate face-based matrix.
    Teuchos::SerialDenseMatrix<int, double> Bpp(mcells, mcells);
    for (int n=0; n!=mcells; ++n) {
      int c = cells[n];
      cells_GID[n] = cmap_wghost.GID(c);

      mesh_->cell_get_faces(c, &faces);
      Bpp(n, n) = Dcc_c[0][c] / faces.size();
      //      AMANZI_ASSERT(std::abs(Dcc_c[0][c] / faces.size()) < 1.e40);

      if (c < ncells_owned) {
        int i = FindPosition<AmanziMesh::Entity_ID>(faces, f);
        AMANZI_ASSERT(i>=0);
        Acf_copy[n] = Acf_cells_[c][i];
        Afc_copy[n] = Afc_cells_[c][i];
      } else {
        Acf_copy[n] = Acf_parallel_f[0][f];
        Afc_copy[n] = Afc_parallel_f[0][f];
      }
    }

    for (int n = 0; n < mcells; n++) {
      for (int m = 0; m < mcells; m++) {
        if (std::abs(Dff_f[0][f] * Afc_copy[m]) == 0.) {
          AMANZI_ASSERT(std::abs(Acf_copy[n]) == 0.);
        } else {
          Bpp(n, m) -= Acf_copy[n] / Dff_f[0][f] * Afc_copy[m];
        }
        //AMANZI_ASSERT(std::abs(Acf_copy[n] / Dff_f[0][f] * Afc_copy[m]) < 1.e40);
      }
    }

    (*App_).SumIntoGlobalValues(mcells, cells_GID, Bpp.values());
  }
  (*App_).GlobalAssemble();
  assembled_app_ = true;





  // Check min
#ifdef ENABLE_DBC
  Epetra_Vector App_diag(mesh_->cell_map(false));
  App_->ExtractDiagonalCopy(App_diag);
  double minval;
  double maxval;
  App_diag.MinValue(&minval);
  App_diag.MaxValue(&maxval);
  AMANZI_ASSERT(std::abs(minval) < 1.e40);
  AMANZI_ASSERT(std::abs(maxval) < 1.e40);
#endif

  // std::stringstream filename_s;
  // filename_s << "App_" << 0 << ".txt";
  // EpetraExt::RowMatrixToMatlabFile(filename_s.str().c_str(), *App_);

}


/* ******************************************************************
 * Parallel matvec product App * Xc.
 ****************************************************************** */
int MatrixMFD_TPFA::Apply(const CompositeVector& X,
			   CompositeVector& Y) const {
  if (!assembled_app_) AssembleApp_();

  int ierr = App_->Multiply(false, *X.ViewComponent("cell",false),
                            *Y.ViewComponent("cell",false));

  if (Y.HasComponent("face")) {
    ApplyAfc(X, Y, 0.);
    ierr |= Y.ViewComponent("face",false)->Multiply(1.,
	 *Dff_->ViewComponent("face",false), *X.ViewComponent("face",false), 1.);
  }

  if (ierr) {
    Errors::Message msg("MatrixMFD_TPFA::Apply has failed to calculate y = A*x.");
    Exceptions::amanzi_throw(msg);
  }
  return ierr;
}

// Apply, for the cell-only system
int MatrixMFD_TPFA::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {
  if (!assembled_app_) AssembleApp_();

  int ierr = App_->Multiply(false, X, Y);
  if (ierr) {
    Errors::Message msg("MatrixMFD_TPFA::Apply has failed to calculate y = A*x.");
    Exceptions::amanzi_throw(msg);
  }
  return ierr;
}




int MatrixMFD_TPFA::ApplyInverse(const CompositeVector& X,
				 CompositeVector& Y) const {
  if (!assembled_schur_) {
    AssembleSchur_();
    UpdatePreconditioner_();
  }

  if (S_pc_ == Teuchos::null) {
    Errors::Message msg("MatrixMFD_TPFA::ApplyInverse called but no preconditioner sublist was provided");
    Exceptions::amanzi_throw(msg);
  }

  // Solve the Schur complement system App * Yc = Xc.
  int ierr = 0;
  const Epetra_MultiVector& Xc = *X.ViewComponent("cell",false);
  Epetra_MultiVector Tc(Xc);

  // Solve the pp system
  ierr = S_pc_->ApplyInverse(Xc, Tc);
  AMANZI_ASSERT(!ierr);

  *Y.ViewComponent("cell",false) = Tc;

  if (ierr) {
    Errors::Message msg("MatrixMFD_TPFA::ApplyInverse has failed in calculating y = inv(A)*x.");
    Exceptions::amanzi_throw(msg);
  }

  // Update the faces
  if (Y.HasComponent("face")) {
    const Epetra_MultiVector& Xf = *X.ViewComponent("face", false);
    const Epetra_MultiVector& Dff_f = *Dff_->ViewComponent("face",false);

    ierr |= ApplyAfc_(Tc, Y, 0.);  // Afc is kept in the transpose form.
    AMANZI_ASSERT(!ierr);
    Epetra_MultiVector& Yf = *Y.ViewComponent("face", false);

    Yf.Update(1., Xf, -1.);

    int nfaces = Yf.MyLength();
    for (int f=0; f!=nfaces; ++f) {
      Yf[0][f] /= Dff_f[0][f];
    }
  }
  // std::cout<<*X.ViewComponent("cell",false)<<*Y.ViewComponent("cell",false);
  // std::cout<<*X.ViewComponent("face",false)<<*Y.ViewComponent("face",false);


  // CompositeVector Z(X);

  // Z.PutScalar(0.);
  // Apply(Y,Z);

  // std::cout<<*Z.ViewComponent("cell",false)<<*Z.ViewComponent("boundary_face",false);




  //  exit(0);


  return ierr;
}


void MatrixMFD_TPFA::UpdateConsistentFaceConstraints(
    const Teuchos::Ptr<CompositeVector>& u) {
  if (!assembled_rhs_) AssembleRHS_();
  if (!assembled_dff_) AssembleDff_();

  Teuchos::Ptr<const CompositeVector> rhs = rhs_.ptr();
  const Epetra_MultiVector& rhs_f = *rhs->ViewComponent("face", false);

  Teuchos::RCP<CompositeVector> work =
      Teuchos::rcp(new CompositeVector(*rhs_));

  ApplyAfc(*u, *work, 0.);  // Afc is kept in the transpose form.
  Epetra_MultiVector& work_f = *work->ViewComponent("face",false);
  work_f.Update(1.0, rhs_f, -1.0);

  Teuchos::Ptr<const CompositeVector> Dff = Dff_.ptr();
  const Epetra_MultiVector& Dff_f = *Dff->ViewComponent("face",false);

  Epetra_MultiVector& uf = *u->ViewComponent("face",false);
  for (int f=0; f!=uf.MyLength(); ++f) {
    uf[0][f] = work_f[0][f] / Dff_f[0][f];
  }
}

void MatrixMFD_TPFA::UpdateConsistentFaceCorrection(const CompositeVector& u,
    const Teuchos::Ptr<CompositeVector>& Pu) {
  if (!assembled_rhs_) AssembleRHS_();
  if (!assembled_dff_) AssembleDff_();

  ApplyAfc(*Pu, *Pu, 0.);
  Epetra_MultiVector& Pu_f = *Pu->ViewComponent("face",false);
  Pu_f.Update(1., *u.ViewComponent("face",false), -1);

  Teuchos::Ptr<const CompositeVector> Dff = Dff_.ptr();
  const Epetra_MultiVector& Dff_f = *Dff->ViewComponent("face",false);

  for (int f=0; f!=Pu_f.MyLength(); ++f) {
    Pu_f[0][f] /= Dff_f[0][f];
  }
}


// MANY ASSUMPTIONS!  THIS IS NOT GENERAL!
// -- K_abs = 1
// -- upwinding is potential difference, with height overlap
void MatrixMFD_TPFA::AnalyticJacobian(const Upwinding& upwinding,
				      const Teuchos::Ptr<State>& S,
				      std::string potential_key,
				      const CompositeVector& dconductivity,
				      const std::vector<MatrixBC>& bc_markers,
				      const std::vector<double>& bc_values) {
  if (!assembled_app_) AssembleApp_();

  // maps and counts
  AmanziMesh::Entity_ID_List faces;
  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  const Epetra_Map& cmap_wghost = mesh_->cell_map(true);

  // local work arrays
  AmanziMesh::Entity_ID_List cells;
  int cells_GID[2];
  int ierr = 0;

  // Get the derivatives
  std::vector<Teuchos::RCP<Teuchos::SerialDenseMatrix<int, double> > > Jpp_faces;
  upwinding.UpdateDerivatives(S, potential_key, dconductivity, bc_markers, bc_values, &Jpp_faces);
  AMANZI_ASSERT(Jpp_faces.size() == nfaces_owned);

  // Assemble into App
  for (unsigned int f=0; f!=nfaces_owned; ++f) {
    mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);

    int mcells = cells.size();
    for (int n=0; n!=mcells; ++n) {
      cells_GID[n] = cmap_wghost.GID(cells[n]);
    }
    ierr = (*App_).SumIntoGlobalValues(mcells, cells_GID, Jpp_faces[f]->values());
    AMANZI_ASSERT(!ierr);
  }

  // finish assembly
  ierr = App_->GlobalAssemble();
  AMANZI_ASSERT(!ierr);
}

}  // namespace AmanziFlow
}  // namespace Amanzi
