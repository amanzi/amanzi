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
  assembled_schur_ = false;
  assembled_operator_ = false;

  // communicate as necessary
  if (Krel.get() && Krel->HasComponent("face"))
    Krel->ScatterMasterToGhosted("face");

  int dim = mesh_->space_dimension();
  WhetStone::MFD3D_Diffusion mfd(mesh_);
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

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
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
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

  int nfaces = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  for (int f = 0; f < nfaces; f++) {
    mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
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

  Spp_ = Teuchos::rcp(new Epetra_FECrsMatrix(Copy, pp_graph));
  Spp_->GlobalAssemble();
}


/* ******************************************************************
 * Convert elemental mass matrices into stiffness matrices and
 * assemble them into four global matrices.
 * We need an auxiliary GHOST-based vector to assemble the RHS.
 ****************************************************************** */
void MatrixMFD_TPFA::AssembleGlobalMatrices() {
  MatrixMFD::AssembleGlobalMatrices();

  std::vector<std::string> names_c(1,"cell");
  std::vector<AmanziMesh::Entity_kind> locations_c(1,AmanziMesh::CELL);
  std::vector<std::string> names_f(1,"face");
  std::vector<AmanziMesh::Entity_kind> locations_f(1,AmanziMesh::FACE);
  std::vector<int> ndofs(1,1);

  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;
  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);

  // Dff
  Epetra_MultiVector& Dff_f = *Dff_->ViewComponent("face",true);
  Dff_f.PutScalar(0.);
  for (int c=0; c!=ncells_owned; ++c) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);

    int mfaces = faces.size();
    for (int n=0; n!=mfaces; ++n) {
      int f = faces[n];
      Dff_f[0][f] += Aff_cells_[c](n, n);
    }
  }
  Dff_->GatherGhostedToMaster();

  // convert right-hand side to a cell-based vector
  const Epetra_Map& fmap_wghost = mesh_->face_map(true);
  const Epetra_Map& cmap = mesh_->cell_map(false);
  const Epetra_Map& cmap_wghost = mesh_->cell_map(true);
  Epetra_Vector Tc(cmap);
  int nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);

  Teuchos::RCP<Epetra_MultiVector> rhs_faces = rhs_->ViewComponent("face",false);
  Teuchos::RCP<Epetra_MultiVector> rhs_cells = rhs_->ViewComponent("cell",false);

  // Schur complement RHS
  for (int f=0; f!=nfaces_owned; ++f) {
    (*rhs_faces)[0][f] /= Dff_f[0][f];
  }
  (*Acf_).Multiply(false, *rhs_faces, Tc);
  for (int c=0; c!=ncells_owned; ++c) (*rhs_cells)[0][c] -= Tc[c];

  // Zero the face RHS out, as the application of matrix results in zero in
  // the face system.
  rhs_faces->PutScalar(0.);

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
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);

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
  Spp_->PutScalar(0.0);
  for (AmanziMesh::Entity_ID f=0; f!=nfaces_owned; ++f) {
    mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
    int mcells = cells.size();

    // populate face-based matrix.
    Teuchos::SerialDenseMatrix<int, double> Bpp(mcells, mcells);
    for (int n=0; n!=mcells; ++n) {
      int c = cells[n];
      cells_GID[n] = cmap_wghost.GID(c);

      mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
      Bpp(n, n) = Dcc_c[0][c] / faces.size();
      ASSERT(std::abs(Dcc_c[0][c] / faces.size()) < 1.e40);

      if (c < ncells_owned) {
        int i = FindPosition<AmanziMesh::Entity_ID>(faces, f);
        ASSERT(i>=0);
        Acf_copy[n] = Acf_cells_[c][i];
        Afc_copy[n] = Afc_cells_[c][i];
      } else {
        Acf_copy[n] = Acf_parallel_f[0][f];
        Afc_copy[n] = Afc_parallel_f[0][f];
      }
    }

    for (int n = 0; n < mcells; n++) {
      for (int m = 0; m < mcells; m++) {
        Bpp(n, m) -= Acf_copy[n] / Dff_f[0][f] * Afc_copy[m];
	ASSERT(std::abs(Acf_copy[n] / Dff_f[0][f] * Afc_copy[m]) < 1.e40);
      }
    }

    (*Spp_).SumIntoGlobalValues(mcells, cells_GID, Bpp.values());
  }
  (*Spp_).GlobalAssemble();

  // Check min
#ifdef ENABLE_DBC
  Epetra_Vector Spp_diag(mesh_->cell_map(false));
  Spp_->ExtractDiagonalCopy(Spp_diag);
  double minval;
  double maxval;
  Spp_diag.MinValue(&minval);
  Spp_diag.MaxValue(&maxval);
  ASSERT(std::abs(minval) < 1.e40);
  ASSERT(std::abs(maxval) < 1.e40);
#endif
}


void
MatrixMFD_TPFA::ComputeSchurComplement(const std::vector<MatrixBC>& bc_markers,
				       const std::vector<double>& bc_values) {
  AssertAssembledOperator_or_die_();
  assembled_schur_ = true;
  Sff_ = Spp_;
}


/* ******************************************************************
 * Parallel matvec product Spp * Xc.
 ****************************************************************** */
int MatrixMFD_TPFA::Apply(const CompositeVector& X,
			   CompositeVector& Y) const {
  AssertAssembledOperator_or_die_();

  int ierr = Spp_->Multiply(false, *X.ViewComponent("cell",false),
                            *Y.ViewComponent("cell",false));

  if (Y.HasComponent("face")) {
    Epetra_MultiVector& Yf = *Y.ViewComponent("face",false);
    ierr |= Afc_->Multiply(true, *X.ViewComponent("cell",false), Yf);
    ierr |= Yf.Multiply(1., *Dff_->ViewComponent("face",false), *X.ViewComponent("face",false), 1.);
  }

  if (ierr) {
    Errors::Message msg("MatrixMFD_TPFA::Apply has failed to calculate y = A*x.");
    Exceptions::amanzi_throw(msg);
  }
  return ierr;
}

// Apply, for the cell-only system
int MatrixMFD_TPFA::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {
  AssertAssembledOperator_or_die_();

  int ierr = Spp_->Multiply(false, X, Y);
  if (ierr) {
    Errors::Message msg("MatrixMFD_TPFA::Apply has failed to calculate y = A*x.");
    Exceptions::amanzi_throw(msg);
  }
  return ierr;
}




int MatrixMFD_TPFA::ApplyInverse(const CompositeVector& X,
				 CompositeVector& Y) const {
  AssertAssembledOperator_or_die_();
  AssertAssembledSchur_or_die_();

  // Solve the Schur complement system Spp * Yc = Xc.
  int ierr = 0;
  const Epetra_MultiVector& Xc = *X.ViewComponent("cell",false);
  Epetra_MultiVector Tc(Xc);

  // Solve the pp system
  ierr = S_pc_->ApplyInverse(Xc, Tc);
  ASSERT(!ierr);

  *Y.ViewComponent("cell",false) = Tc;

  if (ierr) {
    Errors::Message msg("MatrixMFD_TPFA::ApplyInverse has failed in calculating y = inv(A)*x.");
    Exceptions::amanzi_throw(msg);
  }

  // Update the faces
  if (Y.HasComponent("face")) {
    Epetra_MultiVector& Yf = *Y.ViewComponent("face", false);
    const Epetra_MultiVector& Xf = *X.ViewComponent("face", false);
    const Epetra_MultiVector& Dff_f = *Dff_->ViewComponent("face",false);

    Afc_->Multiply(true, Tc, Yf);  // Afc is kept in the transpose form.
    Yf.Update(1., Xf, -1.);

    int nfaces = Yf.MyLength();
    for (int f=0; f!=nfaces; ++f) {
      Yf[0][f] /= Dff_f[0][f];
    }
  }
  return ierr;
}


void MatrixMFD_TPFA::UpdateConsistentFaceConstraints(
    const Teuchos::Ptr<CompositeVector>& u) {
  AssertAssembledOperator_or_die_();

  Teuchos::RCP<const Epetra_MultiVector> uc = u->ViewComponent("cell", false);
  Epetra_MultiVector& uf = *u->ViewComponent("face", false);

  Epetra_MultiVector& Dff_f = *Dff_->ViewComponent("face",false);
  Epetra_MultiVector& rhs_f = *rhs_->ViewComponent("face", false);
  Epetra_MultiVector update_f(rhs_f);

  Afc_->Multiply(true,*uc, update_f);  // Afc is kept in the transpose form.
  update_f.Update(1.0, rhs_f, -1.0);

  int nfaces = rhs_f.MyLength();
  for (int f=0; f!=nfaces; ++f) {
    uf[0][f] = update_f[0][f] / Dff_f[0][f];
  }
}

void MatrixMFD_TPFA::UpdateConsistentFaceCorrection(const CompositeVector& u,
    const Teuchos::Ptr<CompositeVector>& Pu) {
  AssertAssembledOperator_or_die_();

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


// MANY ASSUMPTIONS!  THIS IS NOT GENERAL!
// -- K_abs = 1
// -- upwinding is potential difference, with height overlap
void MatrixMFD_TPFA::AnalyticJacobian(const Upwinding& upwinding,
				      const Teuchos::Ptr<State>& S,
				      std::string potential_key,
				      const CompositeVector& dconductivity,
				      const std::vector<MatrixBC>& bc_markers,
				      const std::vector<double>& bc_values) {

  AssertAssembledOperator_or_die_();

  // maps and counts
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;
  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  int nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  const Epetra_Map& cmap_wghost = mesh_->cell_map(true);

  // local work arrays
  AmanziMesh::Entity_ID_List cells;
  int cells_GID[2];
  int ierr = 0;

  // Get the derivatives
  std::vector<Teuchos::RCP<Teuchos::SerialDenseMatrix<int, double> > > Jpp_faces;
  upwinding.UpdateDerivatives(S, potential_key, dconductivity, bc_markers, bc_values, &Jpp_faces);
  ASSERT(Jpp_faces.size() == nfaces_owned);

  // Assemble into Spp
  for (unsigned int f=0; f!=nfaces_owned; ++f) {
    mesh_->face_get_cells(f, AmanziMesh::USED, &cells);

    int mcells = cells.size();
    for (int n=0; n!=mcells; ++n) {
      cells_GID[n] = cmap_wghost.GID(cells[n]);
    }
    ierr = (*Spp_).SumIntoGlobalValues(mcells, cells_GID, Jpp_faces[f]->values());
    ASSERT(!ierr);
  }

  // finish assembly
  ierr = Spp_->GlobalAssemble();
  ASSERT(!ierr);
}

}  // namespace AmanziFlow
}  // namespace Amanzi
