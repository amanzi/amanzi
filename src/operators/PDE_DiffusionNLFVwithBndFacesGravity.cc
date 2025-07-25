/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatskiy (dasvyat@lanl.gov)
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Operators

*/

#include <vector>

// Amanzi
#include "DenseMatrix.hh"

// Operators
#include "Op_Face_Cell.hh"
#include "PDE_DiffusionNLFVwithBndFacesGravity.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Populate face-based matrices.
****************************************************************** */
void
PDE_DiffusionNLFVwithBndFacesGravity::UpdateMatrices(
  const Teuchos::Ptr<const CompositeVector>& flux,
  const Teuchos::Ptr<const CompositeVector>& u)
{
  // affine map of u. It is equivalent to calculating hydraulic head.
  auto hh = Teuchos::rcp(new CompositeVector(*u));
  Epetra_MultiVector& hh_c = *hh->ViewComponent("cell");
  Epetra_MultiVector& hh_bnd = *hh->ViewComponent("boundary_face");
  const Epetra_MultiVector& u_c = *u->ViewComponent("cell");
  const Epetra_MultiVector& u_bnd = *u->ViewComponent("boundary_face");

  for (int c = 0; c < ncells_owned; ++c) {
    double rho_g = GetDensity(c) * fabs(g_[dim_ - 1]);
    double zc = (mesh_->getCellCentroid(c))[dim_ - 1];
    hh_c[0][c] = u_c[0][c] + rho_g * zc;
    auto faces = mesh_->getCellFaces(c);
    for (auto f : faces) {
      int bf = mesh_->getMap(AmanziMesh::Entity_kind::BOUNDARY_FACE, false)
                 .LID(mesh_->getMap(AmanziMesh::Entity_kind::FACE, false).GID(f));
      if (bf >= 0) {
        double zf = (mesh_->getFaceCentroid(f))[dim_ - 1];
        hh_bnd[0][bf] = u_bnd[0][bf] + rho_g * zf;
      }
    }
  }

  if (!is_scalar_) {
    rho_cv_->ScatterMasterToGhosted("cell");
  }

  PDE_DiffusionNLFVwithBndFaces::UpdateMatrices(flux, hh.ptr());

  // add gravity fluxes to the right-hand side.
  global_op_->rhs()->PutScalarGhosted(0.0);

  const std::vector<int>& bc_model = bcs_trial_[0]->bc_model();
  Epetra_MultiVector& rhs_cell = *global_op_->rhs()->ViewComponent("cell", true);
  Epetra_MultiVector& rhs_bnd = *global_op_->rhs()->ViewComponent("boundary_face");

  for (int f = 0; f < nfaces_owned; ++f) {
    WhetStone::DenseMatrix& Aface = local_op_->matrices[f];

    auto cells = mesh_->getFaceCells(f);
    int ncells = cells.size();

    if (ncells == 2) {
      WhetStone::DenseVector v(ncells), av(ncells);
      for (int n = 0; n < ncells; n++) {
        int c = cells[n];
        double rho_g = GetDensity(c) * fabs(g_[dim_ - 1]);
        double zc = (mesh_->getCellCentroid(c))[dim_ - 1];
        v(n) = zc * rho_g;
      }

      Aface.Multiply(v, av, false);

      for (int n = 0; n < ncells; n++) {
        rhs_cell[0][cells[n]] -= av(n);
      }
    } else if ((bc_model[f] == OPERATOR_BC_DIRICHLET) || (bc_model[f] == OPERATOR_BC_NEUMANN)) {
      int c = cells[0];
      double rho_g = GetDensity(c) * fabs(g_[dim_ - 1]);
      double zf = (mesh_->getFaceCentroid(f))[dim_ - 1];
      double zc = (mesh_->getCellCentroid(c))[dim_ - 1];

      //rhs_cell[0][c] -= Aface(0, 0) * (zc - zf) * rho_g;
      rhs_cell[0][c] -= (Aface(0, 0) * zc + Aface(0, 1) * zf) * rho_g;

      int bf = mesh_->getMap(AmanziMesh::Entity_kind::BOUNDARY_FACE, false)
                 .LID(mesh_->getMap(AmanziMesh::Entity_kind::FACE, false).GID(f));
      rhs_bnd[0][bf] -= (Aface(1, 0) * zc + Aface(1, 1) * zf) * rho_g;
    }
  }

  global_op_->rhs()->GatherGhostedToMaster();
}


/* ******************************************************************
* Calculate flux using cell-centered data.
* **************************************************************** */
void
PDE_DiffusionNLFVwithBndFacesGravity::UpdateFlux(const Teuchos::Ptr<const CompositeVector>& u,
                                                 const Teuchos::Ptr<CompositeVector>& flux)
{
  // Map field u for the local system. For Richards's equation, this
  // is equivalent to calculating the hydraulic head.
  Teuchos::RCP<CompositeVector> hh = Teuchos::rcp(new CompositeVector(*u));
  Epetra_MultiVector& hh_c = *hh->ViewComponent("cell");
  Epetra_MultiVector& hh_bnd = *hh->ViewComponent("boundary_face");
  const Epetra_MultiVector& u_c = *u->ViewComponent("cell");
  const Epetra_MultiVector& u_bnd = *u->ViewComponent("boundary_face");

  for (int c = 0; c < ncells_owned; ++c) {
    double rho_g = GetDensity(c) * fabs(g_[dim_ - 1]);
    double zc = (mesh_->getCellCentroid(c))[dim_ - 1];
    hh_c[0][c] = u_c[0][c] + rho_g * zc;
    auto faces = mesh_->getCellFaces(c);
    for (auto f : faces) {
      int bf = mesh_->getMap(AmanziMesh::Entity_kind::BOUNDARY_FACE, false)
                 .LID(mesh_->getMap(AmanziMesh::Entity_kind::FACE, false).GID(f));
      if (bf >= 0) {
        double zf = (mesh_->getFaceCentroid(f))[dim_ - 1];
        hh_bnd[0][bf] = u_bnd[0][bf] + rho_g * zf;
      }
    }
  }


  PDE_DiffusionNLFVwithBndFaces::UpdateFlux(hh.ptr(), flux);
}


} // namespace Operators
} // namespace Amanzi
