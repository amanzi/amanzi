/*
  This is the flow component of the Amanzi code.
  License: BSD
  Authors: Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
*/

#include "errors.hh"
#include "Epetra_FECrsGraph.h"
#include "EpetraExt_RowMatrixOut.h"
#include "matrix_mfd.hh"
#include "matrix_mfd_scaled_constraint.hh"

namespace Amanzi {
namespace Operators {

bool MatrixMFD_ScaledConstraint::set_symmetric(bool flag_symmetry) {
    if (flag_symmetry) {
      Errors::Message msg("MatrixMFD_ScaledConstraint: ScaledConstraint cannot be symmetric");
      Exceptions::amanzi_throw(msg);
    }
    MatrixMFD::set_symmetric(false);
  }


/* ******************************************************************
 * Allocate global matrices
 ****************************************************************** */
void MatrixMFD_ScaledConstraint::CreateMatrices_(const Epetra_CrsGraph& cf_graph,
        const Epetra_FECrsGraph& ff_graph) {
  MatrixMFD::CreateMatrices_(cf_graph, ff_graph);
  Krel_ = Teuchos::rcp(new Epetra_Vector(*(*rhs_->ViewComponent("face",true))(0)));
}


/* ******************************************************************
 * Calculate elemental stiffness matrices.
 ****************************************************************** */
void MatrixMFD_ScaledConstraint::CreateMFDstiffnessMatrices(
    const Teuchos::Ptr<const CompositeVector>& Krel) {

  if (Krel == Teuchos::null || !Krel->has_component("face")) {
    MatrixMFD::CreateMFDstiffnessMatrices(Krel);

  } else {
    // tag global matrices as invalid
    assembled_schur_ = false;
    assembled_operator_ = false;

    // store a copy of Krel on faces
    *Krel_ = *(*Krel->ViewComponent("face", true))(0);

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

      if (Krel->has_component("cell")) {
        const Epetra_MultiVector& Krel_c = *Krel->ViewComponent("cell",false);

        for (int n=0; n!=nfaces; ++n) {
          for (int m=0; m!=nfaces; ++m) {
            Bff(m, n) = Mff(m,n) * Krel_c[0][c];
          }
        }
      } else {
        for (int n=0; n!=nfaces; ++n) {
          for (int m=0; m!=nfaces; ++m) {
            Bff(m, n) = Mff(m,n);
          }
        }
      }

      const Epetra_MultiVector& Krel_f = *Krel->ViewComponent("face",true);
      double matsum = 0.0;
      for (int n=0; n!=nfaces; ++n) {
        double rowsum = 0.0, colsum = 0.0;
        for (int m=0; m!=nfaces; ++m) {
          colsum += Bff(m,n) * Krel_f[0][faces[m]];
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
}


/* ******************************************************************
 * Applies boundary conditions to elemental stiffness matrices and
 * creates elemental rigth-hand-sides.
 ****************************************************************** */
void MatrixMFD_ScaledConstraint::ApplyBoundaryConditions(
    const std::vector<MatrixBC>& bc_markers,
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

    Teuchos::SerialDenseMatrix<int, double>& Bff = Aff_cells_[c];  // B means elemental.
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
        if (std::abs(bc_values[f]) > 0.) {
          Ff[n] -= bc_values[f] * mesh_->face_area(f) / (*Krel_)[f];
        }
      }
    }
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
void MatrixMFD_ScaledConstraint::DeriveFlux(const CompositeVector& solution,
                           const Teuchos::Ptr<CompositeVector>& flux) const {

  AmanziMesh::Entity_ID_List faces;
  std::vector<double> dp;
  std::vector<int> dirs;

  flux->PutScalar(0.);

  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  int nfaces_owned = flux->size("face",false);

  std::vector<bool> done(nfaces_owned, false);
  const Epetra_MultiVector& soln_cells = *solution.ViewComponent("cell",false);
  const Epetra_MultiVector& soln_faces = *solution.ViewComponent("face",true);
  Epetra_MultiVector& flux_v = *flux->ViewComponent("face",false);

  for (int c=0; c!=ncells; ++c) {
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

        flux_v[0][f] = s * dirs[n] * (*Krel_)[f];
        done[f] = true;
      }
    }
  }

  // ensure post-condition - we got them all
  for (int f=0; f!=nfaces_owned; ++f) {
    ASSERT(done[f]);
  }

}


}  // namespace AmanziFlow
}  // namespace Amanzi
