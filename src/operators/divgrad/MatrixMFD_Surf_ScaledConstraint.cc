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
        const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) :
    MatrixMFD_Surf(plist,mesh),
    MatrixMFD_ScaledConstraint(plist,mesh),
    MatrixMFD(plist,mesh) {}


// Assumes the Surface A was already assembled.
void MatrixMFD_Surf_ScaledConstraint::AssembleAff_() const {
  int ierr(0);

  // Get the standard MFD pieces.
  MatrixMFD_ScaledConstraint::AssembleAff_();

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

  int subsurf_entries = 0;
  int *gsubsurfindices;
  gsubsurfindices = new int[9];
  double *gsubsurfvalues;
  gsubsurfvalues = new double[9];  

  // Loop over surface cells (subsurface faces)
  int nfaces_sub = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  int ncells_surf = surface_mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  for (AmanziMesh::Entity_ID sc=0; sc!=ncells_surf; ++sc) {
    // Access the row in Spp
    AmanziMesh::Entity_ID sc_global = surf_cmap_wghost.GID(sc);
    ierr = Spp.ExtractGlobalRowCopy(sc_global, 9, entries, values, gsurfindices);
    AMANZI_ASSERT(!ierr);

    // Convert Spp global cell numbers to Aff local face numbers
    AmanziMesh::Entity_ID frow = surface_mesh_->entity_get_parent(AmanziMesh::CELL,sc);
    AMANZI_ASSERT(frow < nfaces_sub);
    int frow_global = fmap_wghost.GID(frow);

    int diag = -1;
    for (int m=0; m!=entries; ++m) {
      surfindices[m] = surf_cmap_wghost.LID(gsurfindices[m]);
      indices[m] = surface_mesh_->entity_get_parent(AmanziMesh::CELL,surfindices[m]);
      indices_global[m] = fmap_wghost.GID(indices[m]);
      if (surfindices[m] == sc) diag = m;
    }
    AMANZI_ASSERT(diag > -1);

    // If frozen surface, and surface water, then we want the operator to just
    // have the surface terms.  Clobber Aff's row.
    bool surf_water = std::abs(values[diag]) > 1.e-11;
    if ((*Krel_)[frow] == 0. && surf_water) {
      ierr = Aff_->ExtractGlobalRowCopy(frow_global, 9, subsurf_entries,
              gsubsurfvalues, gsubsurfindices);
      AMANZI_ASSERT(!ierr);
      for (int m=0; m!=subsurf_entries; ++m)
        gsubsurfvalues[m] = 0.;
      ierr = Aff_->InsertGlobalValues(frow_global, subsurf_entries,
              gsubsurfvalues, gsubsurfindices);
      AMANZI_ASSERT(!ierr);
    }

    // If not frozen surface, Krel > 0, and we can use the scaled constraint.
    // Adjust for Krel factor.
    if ((*Krel_)[frow] > 0.) {
      for (int m=0; m!=entries; ++m) {
        // SCALING: -- divide by rel perm to keep constraint equation consistent
        if (std::abs(values[m]) > 0.) {
          values[m] /= (*Krel_)[frow];
        }
      }
    }

    // Add contributions into Aff
    ierr = Aff_->SumIntoGlobalValues(frow_global, entries, values, indices_global);
    AMANZI_ASSERT(!ierr);
  }

  ierr = Aff_->GlobalAssemble();
  AMANZI_ASSERT(!ierr);

  delete[] gsurfindices;
  delete[] surfindices;
  delete[] indices;
  delete[] indices_global;
  delete[] gvalues;
  delete[] values;
  delete[] gsubsurfindices;
  delete[] gsubsurfvalues;
}

void MatrixMFD_Surf_ScaledConstraint::AssembleRHS_() const {
  // MFD portion
  MatrixMFD::AssembleRHS_();

  // Add in surf portion
  const Epetra_MultiVector& rhs_surf_cells =
      *surface_A_->rhs()->ViewComponent("cell",false);
  const Epetra_MultiVector& rhs_faces =
      *rhs_->ViewComponent("face",false);

  int ncells_surf = surface_mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  for (AmanziMesh::Entity_ID sc=0; sc!=ncells_surf; ++sc) {
    AmanziMesh::Entity_ID f = surface_mesh_->entity_get_parent(AmanziMesh::CELL,sc);
    if (std::abs(rhs_surf_cells[0][sc]) > 0.) {
      AMANZI_ASSERT( (*Krel_)[f] > 0. );
      rhs_faces[0][f] += rhs_surf_cells[0][sc] / (*Krel_)[f];
    }
  }
}

void MatrixMFD_Surf_ScaledConstraint::ApplyBoundaryConditions(
    const std::vector<MatrixBC>& bc_markers,
    const std::vector<double>& bc_values,
    bool APPLY_BC_FLUX) {
  // Ensure that none of the surface faces have a BC in them.
  std::vector<MatrixBC> new_markers(bc_markers);

  int ncells_surf = surface_mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  for (AmanziMesh::Entity_ID sc=0; sc!=ncells_surf; ++sc) {
    AmanziMesh::Entity_ID f = surface_mesh_->entity_get_parent(AmanziMesh::CELL,sc);
    new_markers[f] = MATRIX_BC_NULL;
  }

  MatrixMFD_ScaledConstraint::ApplyBoundaryConditions(new_markers, bc_values, APPLY_BC_FLUX);
}



void MatrixMFD_Surf_ScaledConstraint::AssembleSchur_() const {
  int ierr(0);
  int ncells_surf = surface_mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int nfaces_sub = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);

  // Call base Schur
  MatrixMFD_ScaledConstraint::AssembleSchur_();

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

  int subsurf_entries = 0;
  int *gsubsurfindices;
  gsubsurfindices = new int[9];
  double *gsubsurfvalues;
  gsubsurfvalues = new double[9];  

  // Loop over surface cells (subsurface faces)
  for (AmanziMesh::Entity_ID sc=0; sc!=ncells_surf; ++sc) {
    // Access the row in Spp
    AmanziMesh::Entity_ID sc_global = surf_cmap_wghost.GID(sc);
    ierr = Spp.ExtractGlobalRowCopy(sc_global, 9, entries, values, indices);
    AMANZI_ASSERT(!ierr);

    // Convert Spp local cell numbers to Sff local face numbers
    AmanziMesh::Entity_ID frow = surface_mesh_->entity_get_parent(AmanziMesh::CELL,sc);
    AMANZI_ASSERT(frow < nfaces_sub);
    AmanziMesh::Entity_ID frow_global = fmap_wghost.GID(frow);

    // first loop sets up indices and determines the diagonal
    int diag = -1;
    for (int m=0; m!=entries; ++m) {
      indices[m] = surf_cmap_wghost.LID(indices[m]);
      if (indices[m] == sc) diag = m;
      indices[m] = surface_mesh_->entity_get_parent(AmanziMesh::CELL,indices[m]);
      indices_global[m] = fmap_wghost.GID(indices[m]);
    }
    AMANZI_ASSERT(diag > -1);

    // If frozen surface, and surface water, then we want the operator to just
    // have the surface terms.  Clobber Aff's row.
    bool surf_water = std::abs(values[diag]) > 1.e-11;
    if ((*Krel_)[frow] == 0. && surf_water) {
      ierr = Sff_->ExtractGlobalRowCopy(frow_global, 9, subsurf_entries,
              gsubsurfvalues, gsubsurfindices);
      AMANZI_ASSERT(!ierr);
      for (int m=0; m!=subsurf_entries; ++m)
        gsubsurfvalues[m] = 0.;
      ierr = Sff_->InsertGlobalValues(frow_global, subsurf_entries,
              gsubsurfvalues, gsubsurfindices);
      AMANZI_ASSERT(!ierr);
    }
    

    // If not frozen surface, Krel > 0, and we can use the scaled constraint.
    // Adjust for Krel factor.
    if ((*Krel_)[frow] > 0.) {
      for (int m=0; m!=entries; ++m) {
        // SCALING: -- divide by rel perm to keep constraint equation consistent
        if (std::abs(values[m]) > 0.) {
          values[m] /= (*Krel_)[frow];
        }
      }
    }

    // Add contributions into Sff
    ierr = Sff_->SumIntoGlobalValues(frow_global, entries, values, indices_global);
    AMANZI_ASSERT(!ierr);
  }

  ierr = Sff_->GlobalAssemble();
  AMANZI_ASSERT(!ierr);

  delete[] indices;
  delete[] indices_global;
  delete[] values;
  delete[] gsubsurfindices;
  delete[] gsubsurfvalues;

  // dump the schur complement
  // std::stringstream filename_s2;
  // filename_s2 << "schur_post_surf_" << 0 << ".txt";
  // EpetraExt::RowMatrixToMatlabFile(filename_s2.str().c_str(), *Sff_);
}


} // namespace
} // namespace
