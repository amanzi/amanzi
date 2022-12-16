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

  double rho_g = rho_ * norm(g_);
  for (int c = 0; c < ncells_owned; ++c) {
    double zc = (mesh_->cell_centroid(c))[dim_ - 1];
    hh_c[0][c] = u_c[0][c] + rho_g * zc;
  }

  PDE_DiffusionNLFV::UpdateMatrices(flux, hh.ptr());

  // add gravity fluxes to the right-hand side.
  global_op_->rhs()->PutScalarGhosted(0.0);

  const std::vector<int>& bc_model = bcs_trial_[0]->bc_model();
  Epetra_MultiVector& rhs_cell = *global_op_->rhs()->ViewComponent("cell", true);

  AmanziMesh::Entity_ID_List cells;

  for (int f = 0; f < nfaces_owned; ++f) {
    WhetStone::DenseMatrix& Aface = local_op_->matrices[f];

    mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
    int ncells = cells.size();

    if (ncells == 2) {
      WhetStone::DenseVector v(ncells), av(ncells);
      for (int n = 0; n < ncells; n++) {
        int c = cells[n];
        double zc = (mesh_->cell_centroid(c))[dim_ - 1];
        v(n) = zc * rho_g;
      }

      Aface.Multiply(v, av, false);

      for (int n = 0; n < ncells; n++) { rhs_cell[0][cells[n]] -= av(n); }
    } else if (bc_model[f] == OPERATOR_BC_DIRICHLET) {
      int c = cells[0];
      double zf = (mesh_->face_centroid(f))[dim_ - 1];
      double zc = (mesh_->cell_centroid(c))[dim_ - 1];
      rhs_cell[0][c] -= Aface(0, 0) * (zc - zf) * rho_g;
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

  double rho_g = rho_ * norm(g_);
  for (int c = 0; c < ncells_owned; ++c) {
    double zc = (mesh_->cell_centroid(c))[dim_ - 1];
    hh_c[0][c] = u_c[0][c] + rho_g * zc;
  }

  PDE_DiffusionNLFV::UpdateFlux(hh.ptr(), flux);
}


/* ******************************************************************
* BCs are typically given in base system and must be re-mapped.
* **************************************************************** */
double
PDE_DiffusionNLFVwithGravity::MapBoundaryValue_(int f, double u)
{
  double rho_g = rho_ * fabs(g_[dim_ - 1]);
  double zf = (mesh_->face_centroid(f))[dim_ - 1];
  return u + rho_g * zf;
}

} // namespace Operators
} // namespace Amanzi
