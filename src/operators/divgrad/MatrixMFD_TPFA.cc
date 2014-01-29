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

    Aff_cells_.push_back(Bff);
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
      }
    }

    (*Spp_).SumIntoGlobalValues(mcells, cells_GID, Bpp.values());
  }
  (*Spp_).GlobalAssemble();
}



/* ******************************************************************
 * Parallel matvec product Spp * Xc.
 ****************************************************************** */
int MatrixMFD_TPFA::Apply(const CompositeVector& X,
			   CompositeVector& Y) const {

  int ierr = (*Spp_).Multiply(false, *X.ViewComponent("cell",false),
          *Y.ViewComponent("cell",false));

  if (Y.HasComponent("face")) {
    Epetra_MultiVector& Yf = *Y.ViewComponent("face",false);
    ierr |= Afc_->Multiply(true, *X.ViewComponent("cell",false), Yf);
    std::cout << std::setprecision(15) << "FUCK RESIDUAL LHS Afc*Xc(237): LHS = " << Yf[0][237] << std::endl;
    ierr |= Yf.Multiply(1., *Dff_->ViewComponent("face",false), *X.ViewComponent("face",false), 1.);
    std::cout << std::setprecision(15) << "FUCK RESIDUAL LHS Dff*Xf(237): LHS = " << (*Dff_->ViewComponent("face",false))[0][237]*(*X.ViewComponent("face",false))[0][237] << std::endl;
  }


  std::cout << std::setprecision(15) << "FUCK RESIDUAL LHS F(237): LHS = " << (*Y.ViewComponent("face",false))[0][237] << std::endl;
  std::cout << std::setprecision(15) << "FUCK RESIDUAL RHS F(237): RHS = " << (*rhs_->ViewComponent("face",false))[0][237] << std::endl;
  

  if (ierr) {
    Errors::Message msg("MatrixMFD_TPFA::Apply has failed to calculate y = A*x.");
    Exceptions::amanzi_throw(msg);
  }
  return ierr;
}



int MatrixMFD_TPFA::ApplyInverse(const CompositeVector& X,
				 CompositeVector& Y) const {
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
void MatrixMFD_TPFA::AnalyticJacobian(const CompositeVector& height,
        const CompositeVector& potential,
        const CompositeVector& Krel,
        const CompositeVector& dKrel_dp,
        const CompositeVector& Krel_cell,
        const CompositeVector& dKrel_cell_dp) {

  // maps and counts
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;
  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  int nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  const Epetra_Map& cmap_wghost = mesh_->cell_map(true);

  // local work arrays
  AmanziMesh::Entity_ID_List cells;
  int cells_LID[2], cells_GID[2];
  double perm_abs_vert[2];
  double perm_abs_horz[2];
  double k_rel[2];
  double dkrel_dp[2];
  double k_rel_cell[2];
  double dkrel_cell_dp[2];
  double height_l[2];
  double potential_l[2];
  AmanziGeometry::Point cntr_cell[2];
  double dist = 0;

  // local pointers
  height.ScatterMasterToGhosted();
  potential.ScatterMasterToGhosted();
  Krel.ScatterMasterToGhosted();
  dKrel_dp.ScatterMasterToGhosted();
  Krel_cell.ScatterMasterToGhosted();
  dKrel_cell_dp.ScatterMasterToGhosted();
  const Epetra_MultiVector& height_c = *height.ViewComponent("cell",true);
  const Epetra_MultiVector& potential_c = *potential.ViewComponent("cell",true);
  const Epetra_MultiVector& Krel_c = *Krel.ViewComponent("cell",true);
  const Epetra_MultiVector& dKrel_dp_c = *dKrel_dp.ViewComponent("cell",true);
  const Epetra_MultiVector& Krel_cell_c = *Krel_cell.ViewComponent("cell",true);
  const Epetra_MultiVector& dKrel_cell_dp_c = *dKrel_cell_dp.ViewComponent("cell",true);

  for (int f=0; f!=nfaces_owned; ++f){
    mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
    int mcells = cells.size();
    Teuchos::SerialDenseMatrix<int, double> Jpp(mcells, mcells);
    AmanziGeometry::Point face_cntr = mesh_->face_centroid(f);

    for (int n=0; n!=mcells; ++n) {
      cells_LID[n] = cells[n];
      cells_GID[n] = cmap_wghost.GID(cells_LID[n]);

      height_l[n] = height_c[0][cells_LID[n]];
      potential_l[n] = potential_c[0][cells_LID[n]];
      k_rel[n] = Krel_c[0][cells_LID[n]];
      dkrel_dp[n] = dKrel_dp_c[0][cells_LID[n]];
      k_rel_cell[n] = Krel_cell_c[0][cells_LID[n]];
      dkrel_cell_dp[n] = dKrel_cell_dp_c[0][cells_LID[n]];
      cntr_cell[n] = mesh_->cell_centroid(cells_LID[n]);
    }

    if (mcells == 2) {
      dist = norm(cntr_cell[0] - cntr_cell[1]);
    } else if (mcells == 1) {
      dist = norm(cntr_cell[0] - face_cntr);
      height_l[1] = height("face",f);
      potential_l[1] = potential("face",f);
      k_rel[1] = Krel_cell("face",f);
      dkrel_dp[1] = 0.;
      k_rel_cell[1] = 0.;
      dkrel_cell_dp[1] = 0.;

    } else {
      ASSERT(0);
    }

    AmanziGeometry::Point normal = mesh_->face_normal(f, false, cells_LID[0]);
    normal *= 1./ mesh_->face_area(f);

    ComputeJacobianLocal(mcells, f, dist, height_l, potential_l,
                         k_rel, dkrel_dp, k_rel_cell, dkrel_cell_dp, Jpp);

    (*Spp_).SumIntoGlobalValues(mcells, cells_GID, Jpp.values());
  }

  Spp_->GlobalAssemble();
}


/* ******************************************************************
* Computation of a local submatrix of 
* Analytical Jacobian (nonlinear part) on a particular face.
// MANY ASSUMPTIONS!  THIS IS NOT GENERAL!
// -- K_abs = 1
// -- upwinding is potential difference, with height overlap
****************************************************************** */
void MatrixMFD_TPFA::ComputeJacobianLocal(int mcells,
        int face_id,
        double dist,
        double *height,
        double *potential,
        double *k_rel,
        double *dk_rel_dp,
        double *k_rel_cell,
        double *dk_rel_cell_dp,
        Teuchos::SerialDenseMatrix<int, double>& Jpp) {

  // Determine Kface (via upwinding) and dKface_dp for each cell.
  double dKface_dp[2];
  double Kface;
  if (mcells == 1) {
    if (potential[0] >= potential[1]) {
      Kface = k_rel[0];
      dKface_dp[0] = dk_rel_dp[0];
    } else {
      Kface = k_rel[1];
      dKface_dp[0] = 0.;
    }
  } else {
    // determine the overlap tolerance
    double ol0 = std::max(0., height[0]);
    double ol1 = std::max(0., height[1]);
    double flow_eps = 0.;
    if ((ol0 > 0) || (ol1 > 0)) {
      flow_eps = (ol0 * ol1) / (ol0 + ol1);
    }
    flow_eps = std::max(flow_eps, 1.e-16);

    // upwind
    if (potential[0] - potential[1] > flow_eps) {
      // 0 is up
      Kface = k_rel[0];
      dKface_dp[0] = dk_rel_dp[0];
      dKface_dp[1] = 0.;
    } else if (potential[1] - potential[0] > flow_eps) {
      // 1 is up
      Kface = k_rel[1];
      dKface_dp[0] = 0.;
      dKface_dp[1] = dk_rel_dp[1];
    } else {
      // Parameterization of a linear scaling between 0 and 1
      double param;
      if (flow_eps < 2*1.e-16) {
        param = 0.5;
      } else {
        param = (potential[1] - potential[0]) / (2*flow_eps) + 0.5;
      }
      ASSERT(param >= 0.0);
      ASSERT(param <= 1.0);

      Kface = param * k_rel[1] + (1.-param) * k_rel[0];
      dKface_dp[0] = (1-param) * dk_rel_dp[0];
      dKface_dp[1] = param * dk_rel_dp[1];
    }
  }

  // Calculate J
  double dphi = (potential[0] - potential[1]) / dist;
  if (mcells == 1) {
    Jpp(0,0) = dphi * (Kface * dk_rel_cell_dp[0] + dKface_dp[0] * k_rel_cell[0])
        * mesh_->face_area(face_id);
  } else {
    Jpp(0, 0) = dphi * (Kface * dk_rel_cell_dp[0] + dKface_dp[0] * k_rel_cell[0])
        * mesh_->face_area(face_id);
    Jpp(0, 1) = dphi * (Kface * dk_rel_cell_dp[1] + dKface_dp[1] * k_rel_cell[1])
        * mesh_->face_area(face_id);

    Jpp(1, 0) = -Jpp(0, 0);
    Jpp(1, 1) = -Jpp(0, 1);
  }
}

}  // namespace AmanziFlow
}  // namespace Amanzi
