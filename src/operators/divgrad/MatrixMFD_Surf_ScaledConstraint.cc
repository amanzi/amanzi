/*
  License: BSD
  Ethan Coon (ecoon@lanl.gov)

  MatrixMFD_Surf provides for solving the p-lambda system where a subset of
  the lambdas are more densely coupled for overland flow.

*/
#include "Epetra_FECrsGraph.h"
#include "EpetraExt_RowMatrixOut.h"

#include "errors.hh"
#include "MatrixMFD_Surf_ScaledConstraint.hh"

#define DEBUG_FLAG 0

namespace Amanzi {
namespace Operators {


MatrixMFD_Surf_ScaledConstraint::MatrixMFD_Surf_ScaledConstraint(
        Teuchos::ParameterList& plist,
        const Teuchos::RCP<const AmanziMesh::Mesh> mesh,
        const Teuchos::RCP<const AmanziMesh::Mesh> surface_mesh) :
    MatrixMFD_Surf(plist, mesh, surface_mesh),
    MatrixMFD_ScaledConstraint(plist,mesh),
    MatrixMFD(plist,mesh) {}


MatrixMFD_Surf_ScaledConstraint::MatrixMFD_Surf_ScaledConstraint(
        const MatrixMFD_ScaledConstraint& other,
        const Teuchos::RCP<const AmanziMesh::Mesh> surface_mesh) :
    MatrixMFD_Surf(other, surface_mesh),
    MatrixMFD_ScaledConstraint(other),
    MatrixMFD(other) {}


// Assumes the Surface A was already assembled.
void MatrixMFD_Surf_ScaledConstraint::AssembleGlobalMatrices() {
  int ierr(0);

  // Get the standard MFD pieces.
  MatrixMFD_ScaledConstraint::AssembleGlobalMatrices();

  // Add the TPFA on the surface parts from surface_A.
  const Epetra_Map& surf_cmap_wghost = surface_mesh_->cell_map(true);
  const Epetra_Map& fmap_wghost = mesh_->face_map(true);
  const Epetra_FECrsMatrix& Spp = *surface_A_->TPFA();

  int entries = 0;
  int gentries = 0;
  double *values;
  values = new double[9];
  double *gvalues;
  gvalues = new double[9];
  int *surfindices;
  surfindices = new int[9];
  int *gsurfindices;
  gsurfindices = new int[9];
  int *indices;
  indices = new int[9];
  int *indices_global;
  indices_global = new int[9];

  // Loop over surface cells (subsurface faces)
  int nfaces_sub = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  int ncells_surf = surface_mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (AmanziMesh::Entity_ID sc=0; sc!=ncells_surf; ++sc) {
    // Access the row in Spp
    AmanziMesh::Entity_ID sc_global = surf_cmap_wghost.GID(sc);
    ierr = Spp.ExtractGlobalRowCopy(sc_global, 9, entries, values, gsurfindices);
    ASSERT(!ierr);

    // Convert Spp global cell numbers to Aff local face numbers
    AmanziMesh::Entity_ID frow = surface_mesh_->entity_get_parent(AmanziMesh::CELL,sc);
    ASSERT(frow < nfaces_sub);
    int frow_global = fmap_wghost.GID(frow);

    for (int m=0; m!=entries; ++m) {
      surfindices[m] = surf_cmap_wghost.LID(gsurfindices[m]);
      indices[m] = surface_mesh_->entity_get_parent(AmanziMesh::CELL,surfindices[m]);
      indices_global[m] = fmap_wghost.GID(indices[m]);

      // SCALING: -- divide by rel perm to keep constraint equation consistent
      if (std::abs(values[m]) > 0.) {
        ASSERT( (*Krel_)[frow] > 0.);
        values[m] /= (*Krel_)[frow];
      }
    }

    ierr = Aff_->SumIntoGlobalValues(frow_global, entries, values, indices_global);
    ASSERT(!ierr);
  }

  ierr = Aff_->GlobalAssemble();
  ASSERT(!ierr);

  delete[] gsurfindices;
  delete[] surfindices;
  delete[] indices;
  delete[] indices_global;
  delete[] gvalues;
  delete[] values;


  // Deal with RHS
  const Epetra_MultiVector& rhs_surf_cells =
      *surface_A_->rhs()->ViewComponent("cell",false);
  const Epetra_MultiVector& rhs_faces =
      *rhs_->ViewComponent("face",false);
  for (AmanziMesh::Entity_ID sc=0; sc!=ncells_surf; ++sc) {
    AmanziMesh::Entity_ID f = surface_mesh_->entity_get_parent(AmanziMesh::CELL,sc);
    if (std::abs(rhs_surf_cells[0][sc]) > 0.) {
      ASSERT( (*Krel_)[f] > 0. );
      rhs_faces[0][f] += rhs_surf_cells[0][sc] / (*Krel_)[f];
    }
  }
}


void MatrixMFD_Surf_ScaledConstraint::ApplyBoundaryConditions(
    const std::vector<MatrixBC>& bc_markers,
    const std::vector<double>& bc_values) {
  // Ensure that none of the surface faces have a BC in them.
  std::vector<MatrixBC> new_markers(bc_markers);

  int ncells_surf = surface_mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (AmanziMesh::Entity_ID sc=0; sc!=ncells_surf; ++sc) {
    AmanziMesh::Entity_ID f = surface_mesh_->entity_get_parent(AmanziMesh::CELL,sc);
    new_markers[f] = MATRIX_BC_NULL;
  }

  MatrixMFD_ScaledConstraint::ApplyBoundaryConditions(new_markers, bc_values);
}



void MatrixMFD_Surf_ScaledConstraint::ComputeSchurComplement(const std::vector<MatrixBC>& bc_markers,
        const std::vector<double>& bc_values) {
  std::vector<MatrixBC> new_markers(bc_markers);

  int ierr(0);
  int ncells_surf = surface_mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (AmanziMesh::Entity_ID sc=0; sc!=ncells_surf; ++sc) {
    AmanziMesh::Entity_ID f = surface_mesh_->entity_get_parent(AmanziMesh::CELL,sc);
    new_markers[f] = MATRIX_BC_NULL;
  }

  // Call base Schur
  std::cout << " Acc(99) = " << Acc_cells_[99] << std::endl;
  std::cout << " Aff(501,501) = " << Aff_cells_[99](5,5) << std::endl;
  MatrixMFD_ScaledConstraint::ComputeSchurComplement(new_markers, bc_values);

  //  dump the schur complement
  // std::stringstream filename_s;
  // filename_s << "schur_pre_surf_" << 0 << ".txt";
  // EpetraExt::RowMatrixToMatlabFile(filename_s.str().c_str(), *Sff_);

  // Add the TPFA on the surface parts from surface_A.
  const Epetra_Map& surf_cmap_wghost = surface_mesh_->cell_map(true);
  const Epetra_Map& fmap_wghost = mesh_->face_map(true);

  // ASSUMES surface_A_'s AssembleGlobalMatrices() has been called
  const Epetra_FECrsMatrix& Spp = *surface_A_->TPFA();

  int entries = 0;
  double *values;
  values = new double[9];
  int *indices;
  indices = new int[9];
  int *indices_global;
  indices_global = new int[9];

  // Loop over surface cells (subsurface faces)
  for (AmanziMesh::Entity_ID sc=0; sc!=ncells_surf; ++sc) {
    // Access the row in Spp
    AmanziMesh::Entity_ID sc_global = surf_cmap_wghost.GID(sc);
    ierr = Spp.ExtractGlobalRowCopy(sc_global, 9, entries, values, indices);
    ASSERT(!ierr);

    // Convert Spp local cell numbers to Sff local face numbers
    AmanziMesh::Entity_ID frow = surface_mesh_->entity_get_parent(AmanziMesh::CELL,sc);
    AmanziMesh::Entity_ID frow_global = fmap_wghost.GID(frow);
    for (int m=0; m!=entries; ++m) {
      indices[m] = surf_cmap_wghost.LID(indices[m]);
      indices[m] = surface_mesh_->entity_get_parent(AmanziMesh::CELL,indices[m]);
      indices_global[m] = fmap_wghost.GID(indices[m]);

      // additionally divide by rel perm to keep constraint equation
      // consistent
      if (std::abs(values[m]) > 0.) {
        ASSERT( (*Krel_)[frow] > 0.);
        values[m] /= (*Krel_)[frow];
      }
    }

    ierr = Sff_->SumIntoGlobalValues(frow_global, entries, values, indices_global);
    ASSERT(!ierr);
  }

  ierr = Sff_->GlobalAssemble();
  ASSERT(!ierr);

  delete[] indices;
  delete[] indices_global;
  delete[] values;

  // dump the schur complement
  // std::stringstream filename_s2;
  // filename_s2 << "schur_post_surf_" << 0 << ".txt";
  // EpetraExt::RowMatrixToMatlabFile(filename_s2.str().c_str(), *Sff_);
}


} // namespace
} // namespace
