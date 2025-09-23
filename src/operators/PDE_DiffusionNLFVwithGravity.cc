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
#include "PDE_DiffusionNLFVwithGravity.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Populate face-based matrices.
****************************************************************** */
void
PDE_DiffusionNLFVwithGravity::UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& flux,
                                             const Teuchos::Ptr<const CompositeVector>& u)
{
  // affine map of u. It is equivalent to calculating hydraulic head.
  Teuchos::RCP<CompositeVector> hh = Teuchos::rcp(new CompositeVector(*u));
  Epetra_MultiVector& hh_c = *hh->ViewComponent("cell");
  const Epetra_MultiVector& u_c = *u->ViewComponent("cell");

  Teuchos::RCP<const Epetra_MultiVector> rho_c;
  if (!is_scalar_) rho_c = rho_cv_->ViewComponent("cell", false);
  double gnorm = norm(g_);

  for (int c = 0; c < ncells_owned; ++c) {
    double zc = (mesh_->getCellCentroid(c))[dim_ - 1];
    double rho = (is_scalar_) ? rho_ : (*rho_c)[0][c];
    hh_c[0][c] = u_c[0][c] + rho * gnorm * zc;
  }

  PDE_DiffusionNLFV::UpdateMatrices(flux, hh.ptr());

  // add gravity fluxes to the right-hand side.
  global_op_->rhs()->PutScalarGhosted(0.0);

  const std::vector<int>& bc_model = bcs_trial_[0]->bc_model();
  Epetra_MultiVector& rhs_cell = *global_op_->rhs()->ViewComponent("cell", true);

  for (int f = 0; f < nfaces_owned; ++f) {
    WhetStone::DenseMatrix& Aface = local_op_->matrices[f];

    auto cells = mesh_->getFaceCells(f);
    int ncells = cells.size();

    if (ncells == 2) {
      WhetStone::DenseVector v(ncells), av(ncells);
      for (int n = 0; n < ncells; n++) {
        int c = cells[n];
        double zc = (mesh_->getCellCentroid(c))[dim_ - 1];
        double rho = (is_scalar_) ? rho_ : (*rho_c)[0][c];
        v(n) = zc * rho * gnorm;
      }

      Aface.Multiply(v, av, false);

      for (int n = 0; n < ncells; n++) {
        rhs_cell[0][cells[n]] -= av(n);
      }
    } else if (bc_model[f] == OPERATOR_BC_DIRICHLET) {
      int c = cells[0];
      double zf = (mesh_->getFaceCentroid(f))[dim_ - 1];
      double zc = (mesh_->getCellCentroid(c))[dim_ - 1];
      double rho = (is_scalar_) ? rho_ : (*rho_c)[0][c];
      rhs_cell[0][c] -= Aface(0, 0) * (zc - zf) * rho * gnorm;
    }
  }

  global_op_->rhs()->GatherGhostedToMaster();
}


/* ******************************************************************
* Calculate flux using cell-centered data.
* **************************************************************** */
void
PDE_DiffusionNLFVwithGravity::UpdateFlux(const Teuchos::Ptr<const CompositeVector>& u,
                                         const Teuchos::Ptr<CompositeVector>& flux)
{
  // Map field u for the local system. For Richards's equation, this
  // is equivalent to calculating the hydraulic head.
  Teuchos::RCP<CompositeVector> hh = Teuchos::rcp(new CompositeVector(*u));
  Epetra_MultiVector& hh_c = *hh->ViewComponent("cell");
  const Epetra_MultiVector& u_c = *u->ViewComponent("cell");

  Teuchos::RCP<const Epetra_MultiVector> rho_c;
  if (!is_scalar_) rho_c = rho_cv_->ViewComponent("cell", false);
  double gnorm = norm(g_);

  for (int c = 0; c < ncells_owned; ++c) {
    double zc = (mesh_->getCellCentroid(c))[dim_ - 1];
    double rho = (is_scalar_) ? rho_ : (*rho_c)[0][c];
    hh_c[0][c] = u_c[0][c] + rho * gnorm * zc;
  }

  PDE_DiffusionNLFV::UpdateFlux(hh.ptr(), flux);
}


/* ******************************************************************
* BCs are typically given in base system and must be re-mapped.
* **************************************************************** */
double
PDE_DiffusionNLFVwithGravity::MapBoundaryValue_(int f, double u)
{
  int c = getFaceOnBoundaryInternalCell(*mesh_, f);
  double rho = (is_scalar_) ? rho_ : (*rho_cv_->ViewComponent("cell"))[0][c];
  double g = fabs(g_[dim_ - 1]);
  double zf = (mesh_->getFaceCentroid(f))[dim_ - 1];

  return u + rho * g * zf;
}

} // namespace Operators
} // namespace Amanzi
