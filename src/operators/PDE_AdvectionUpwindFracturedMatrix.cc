/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
           Ethan Coon (ecoon@lanl.gov)
*/

/*
  Operators

*/

#include <vector>

#include "WhetStoneDefs.hh"

#include "OperatorDefs.hh"
#include "Operator_Cell.hh"
#include "Op_Face_Cell.hh"
#include "Op_SurfaceFace_SurfaceCell.hh"
#include "PDE_AdvectionUpwindFracturedMatrix.hh"
#include "UniqueLocalIndex.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Advection requires a velocity field.
****************************************************************** */
void
PDE_AdvectionUpwindFracturedMatrix::Setup(const CompositeVector& u)
{
  IdentifyUpwindCells_(u);
}


/* ******************************************************************
* A simple first-order transport method.
* Advection operator is of the form: div (u C), where u is the given
* velocity field and C is the advected field.
****************************************************************** */
void
PDE_AdvectionUpwindFracturedMatrix::UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& u)
{
  std::vector<WhetStone::DenseMatrix>& matrix = local_op_->matrices;

  const Epetra_MultiVector& uf = *u->ViewComponent("face");

  for (int f = 0; f < nfaces_owned; ++f) {
    int g = gmap_->FirstPointInElement(f);
    int c1 = (*upwind_cell_)[g];
    int c2 = (*downwind_cell_)[g];

    auto cells = mesh_->getFaceCells(f);
    int ncells = cells.size();
    WhetStone::DenseMatrix Aface(ncells, ncells);
    Aface.PutScalar(0.0);

    double umod = fabs(uf[0][g]);
    if (c1 < 0) {
      Aface(0, 0) = umod;
    } else if (c2 < 0) {
      Aface(0, 0) = umod;
    } else {
      int i = (cells[0] == c1) ? 0 : 1;
      Aface(i, i) = umod;
      Aface(1 - i, i) = -umod;
    }

    matrix[f] = Aface;
  }

  // removed matrices on faces where fracture is located
  for (int i = 0; i < fractures_.size(); ++i) {
    auto [block, vofs] = mesh_->getSetEntitiesAndVolumeFractions(
      fractures_[i], AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED);

    for (int n = 0; n < block.size(); ++n) { matrix[block[n]] *= 0.0; }
  }
}


/* ******************************************************************
* Add a simple first-order upwind method where the advected quantity
* is not the primary variable (used in Jacobians).
* Advection operator is of the form: div (q H(u))
*     q:    flux
*     H(u): advected quantity (i.e. enthalpy)
****************************************************************** */
void
PDE_AdvectionUpwindFracturedMatrix::UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& u,
                                                   const Teuchos::Ptr<const CompositeVector>& dhdT)
{
  std::vector<WhetStone::DenseMatrix>& matrix = local_op_->matrices;

  const Epetra_MultiVector& uf = *u->ViewComponent("face");

  dhdT->ScatterMasterToGhosted("cell");
  const Epetra_MultiVector& dh = *dhdT->ViewComponent("cell", true);

  for (int f = 0; f < nfaces_owned; ++f) {
    int g = gmap_->FirstPointInElement(f);
    int c1 = (*upwind_cell_)[g];
    int c2 = (*downwind_cell_)[g];

    auto cells = mesh_->getFaceCells(f);
    int ncells = cells.size();
    WhetStone::DenseMatrix Aface(ncells, ncells);
    Aface.PutScalar(0.0);

    double umod = fabs(uf[0][g]);
    if (c1 < 0) {
      Aface(0, 0) = umod * dh[0][c2];
    } else if (c2 < 0) {
      Aface(0, 0) = umod * dh[0][c1];
    } else {
      int i = (cells[0] == c1) ? 0 : 1;
      Aface(i, i) = umod * dh[0][c1];
      Aface(1 - i, i) = -umod * dh[0][c2];
    }

    matrix[f] = Aface;
  }

  // removed matrices fof faces where fracture is located
  for (int i = 0; i < fractures_.size(); ++i) {
    auto [block, vofs] = mesh_->getSetEntitiesAndVolumeFractions(
      fractures_[i], AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED);

    for (int n = 0; n < block.size(); ++n) { matrix[block[n]] *= 0.0; }
  }
}


/* *******************************************************************
* Apply boundary condition to the local matrices
* Recommended options: primary=true, eliminate=false, essential_eqn=true
*  - must deal with Dirichlet BC on inflow boundary
*  - Dirichlet on outflow boundary is ill-posed
*  - Neumann on inflow boundary is typically not used, since it is
*    equivalent to Dirichlet BC. We perform implicit conversion to
*    Dirichlet BC.
*
* Advection-diffusion problem.
* Recommended options: primary=false, eliminate=true, essential_eqn=true
*  - Dirichlet BC is treated as usual
*  - Neuman on inflow boundary: If diffusion takes care of the total
*    flux, then TOTAL_FLUX model must be used. If diffusion deals
*    with the diffusive flux only (NEUMANN model), value of the
*    advective flux is in general not available and negative value
*    is added to matrix diagonal. The discrete system may lose SPD
*    property.
*  - Neuman on outflow boundary: If diffusion takes care of the total
*    flux, then TOTAL_FLUX model must be used. Otherwise, do nothing.
*
* FIXME: So far we support the case bc_test = bc_trial
******************************************************************* */
void
PDE_AdvectionUpwindFracturedMatrix::ApplyBCs(bool primary, bool eliminate, bool essential_eqn)
{
  std::vector<WhetStone::DenseMatrix>& matrix = local_op_->matrices;
  Epetra_MultiVector& rhs_cell = *global_op_->rhs()->ViewComponent("cell");

  const std::vector<int>& bc_model = bcs_trial_[0]->bc_model();
  const std::vector<double>& bc_value = bcs_trial_[0]->bc_value();

  for (int f = 0; f < nfaces_owned; f++) {
    int g = gmap_->FirstPointInElement(f);
    int ndofs = gmap_->ElementSize(f);

    int c1 = (*upwind_cell_)[g];
    int c2 = (*downwind_cell_)[g];
    if (bc_model[f] == OPERATOR_BC_DIRICHLET) {
      if (c2 < 0) {
        // pass, the upwind cell is internal to the domain, so all is good
      } else if (c1 < 0) {
        // downwind cell is internal to the domain
        rhs_cell[0][c2] += matrix[f](0, 0) * bc_value[f];
        matrix[f] = 0.0;
      }
    }
    // coupling fluxes are separate object
    else if (ndofs == 2) {
      matrix[f] = 0.0;
    }
    // treat as essential inflow BC for pure advection
    else if (bc_model[f] == OPERATOR_BC_NEUMANN && primary) {
      if (c1 < 0) {
        rhs_cell[0][c2] += mesh_->getFaceArea(f) * bc_value[f];
        matrix[f] = 0.0;
      }
    }
    // leave in matrix for composite operator
    else if (bc_model[f] == OPERATOR_BC_NEUMANN && !primary) {
      if (c1 < 0) matrix[f] *= -1.0;
    }
    // total flux was processed by another operator -> remove here
    else if (bc_model[f] == OPERATOR_BC_TOTAL_FLUX && !primary) {
      matrix[f] = 0.0;
    }
    // do not know what to do
    else if (bc_model[f] != OPERATOR_BC_NONE) {
      AMANZI_ASSERT(false);
    }
  }
}


/* *******************************************************************
* Identify flux direction based on orientation of the face normal
* and sign of the  Darcy velocity.
******************************************************************* */
void
PDE_AdvectionUpwindFracturedMatrix::IdentifyUpwindCells_(const CompositeVector& u)
{
  u.ScatterMasterToGhosted("face");
  const Epetra_MultiVector& uf = *u.ViewComponent("face", true);
  gmap_ = u.ComponentMap("face", true);

  upwind_cell_ = Teuchos::rcp(new Epetra_IntVector(*gmap_));
  downwind_cell_ = Teuchos::rcp(new Epetra_IntVector(*gmap_));

  upwind_cell_->PutValue(-1); // negative value indicates boundary
  downwind_cell_->PutValue(-1);

  for (int c = 0; c < ncells_wghost; c++) {
    const auto& [faces, fdirs] = mesh_->getCellFacesAndDirections(c);

    for (int i = 0; i < faces.size(); i++) {
      int f = faces[i];
      int g = gmap_->FirstPointInElement(f);

      int ndofs = gmap_->ElementSize(f);
      if (ndofs == 2) g += UniqueIndexFaceToCells(*mesh_, f, c);

      if (uf[0][g] * fdirs[i] >= 0) {
        (*upwind_cell_)[g] = c;
      } else {
        (*downwind_cell_)[g] = c;
      }
    }
  }
}


/* ******************************************************************
* Initialize additional parameters
****************************************************************** */
void
PDE_AdvectionUpwindFracturedMatrix::InitAdvection_(Teuchos::ParameterList& plist)
{
  fractures_ = plist.get<Teuchos::Array<std::string>>("fracture").toVector();
}

} // namespace Operators
} // namespace Amanzi
