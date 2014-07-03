/*
  License: BSD
  Ethan Coon (ecoon@lanl.gov)

  MatrixMFD_Surf provides for solving the p-lambda system where a subset of
  the lambdas are more densely coupled for overland flow.

*/
#include "Epetra_FECrsGraph.h"
#include "EpetraExt_RowMatrixOut.h"

#include "errors.hh"
#include "MatrixMFD_Surf.hh"

#define DEBUG_FLAG 0

namespace Amanzi {
namespace Operators {

MatrixMFD_Surf::MatrixMFD_Surf(Teuchos::ParameterList& plist,
			       const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) :
    MatrixMFD(plist,mesh) {
  // dump
  dump_schur_ = plist_.get<bool>("dump Schur complement", false);
}

void
MatrixMFD_Surf::SetSurfaceOperator(const Teuchos::RCP<MatrixMFD_TPFA>& surface_A) {
  surface_A_ = surface_A;
  surface_mesh_ = surface_A_->Mesh();
}



int MatrixMFD_Surf::Apply(const CompositeVector& X,
			  CompositeVector& Y) const {
  // Apply the base, subsurface matrix
  int ierr = MatrixMFD::Apply(X,Y);
  ASSERT(!ierr);

  // Manually copy data -- TRILINOS FAIL
  const Epetra_MultiVector& Xf = *X.ViewComponent("face", false);
  Epetra_MultiVector surf_X(surface_mesh_->cell_map(false),1);
  for (int sc=0; sc!=surf_X.MyLength(); ++sc) {
    surf_X[0][sc] = Xf[0][surface_mesh_->entity_get_parent(AmanziMesh::CELL, sc)];
  }
  
  // Apply the surface-only operators, blockwise
  Epetra_MultiVector surf_Y(surface_mesh_->cell_map(false),1);
  ierr |= surface_A_->Apply(surf_X, surf_Y);
  ASSERT(!ierr);

  // Add back into Y
  Epetra_MultiVector& Yf = *Y.ViewComponent("face",false);
  for (int sc=0; sc!=surf_X.MyLength(); ++sc) {
    Yf[0][surface_mesh_->entity_get_parent(AmanziMesh::CELL, sc)] += surf_Y[0][sc];
  }
  return ierr;
}


void MatrixMFD_Surf::FillMatrixGraphs_(const Teuchos::Ptr<Epetra_CrsGraph> cf_graph,
          const Teuchos::Ptr<Epetra_FECrsGraph> ff_graph) {
  // get the standard fill
  MatrixMFD::FillMatrixGraphs_(cf_graph, ff_graph);

  // additional face-to-face connections for the TPF on surface.
  const Epetra_Map& cmap = mesh_->cell_map(false);
  const Epetra_Map& fmap = mesh_->face_map(false);
  const Epetra_Map& fmap_wghost = mesh_->face_map(true);
  const Epetra_Map& surf_cmap_wghost = surface_mesh_->cell_map(true);
  int ierr(0);

  AmanziMesh::Entity_ID_List surf_cells;
  int equiv_face_LID[2];
  int equiv_face_GID[2];

  unsigned int nfaces_surf = surface_mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);

  for (unsigned int fs=0; fs!=nfaces_surf; ++fs) {
    surface_mesh_->face_get_cells(fs, AmanziMesh::USED, &surf_cells);
    int ncells_surf = surf_cells.size(); ASSERT(ncells_surf <= 2);

    // get the equivalent faces on the subsurface mesh
    for (int n=0; n!=ncells_surf; ++n) {
      equiv_face_LID[n] = surface_mesh_->entity_get_parent(AmanziMesh::CELL, surf_cells[n]);
      equiv_face_GID[n] = fmap_wghost.GID(equiv_face_LID[n]);
    }

    // insert the connection
    ierr = ff_graph->InsertGlobalIndices(ncells_surf, equiv_face_GID,
            ncells_surf, equiv_face_GID);
    ASSERT(!ierr);
  }
}


// Assumes the Surface A was already assembled.
void MatrixMFD_Surf::AssembleAff_() const {
  int ierr(0);

  // Get the standard MFD pieces.
  MatrixMFD::AssembleAff_();

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
}


void MatrixMFD_Surf::AssembleRHS_() const {
  // MFD portion
  MatrixMFD::AssembleRHS_();

  // Add in surf portion
  const Epetra_MultiVector& rhs_surf_cells =
      *surface_A_->rhs()->ViewComponent("cell",false);
  const Epetra_MultiVector& rhs_faces =
      *rhs_->ViewComponent("face",false);

  int ncells_surf = surface_mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (AmanziMesh::Entity_ID sc=0; sc!=ncells_surf; ++sc) {
    AmanziMesh::Entity_ID f = surface_mesh_->entity_get_parent(AmanziMesh::CELL,sc);
    rhs_faces[0][f] += rhs_surf_cells[0][sc];
  }
}


void MatrixMFD_Surf::ApplyBoundaryConditions(
    const std::vector<MatrixBC>& bc_markers,
    const std::vector<double>& bc_values,
    bool ADD_BC_FLUX) {
  // Ensure that none of the surface faces have a BC in them.
  std::vector<MatrixBC> new_markers(bc_markers);

  int ncells_surf = surface_mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (AmanziMesh::Entity_ID sc=0; sc!=ncells_surf; ++sc) {
    AmanziMesh::Entity_ID f = surface_mesh_->entity_get_parent(AmanziMesh::CELL,sc);
    new_markers[f] = MATRIX_BC_NULL;
  }

  MatrixMFD::ApplyBoundaryConditions(new_markers, bc_values, ADD_BC_FLUX);
}



void MatrixMFD_Surf::AssembleSchur_() const {
  int ierr(0);

  // Call base Schur
  MatrixMFD::AssembleSchur_();

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
  int ncells_surf = surface_mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
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
  if (dump_schur_) {
    std::stringstream filename_s2;
    filename_s2 << "schur_MatrixMFD_Surf_" << 0 << ".txt";
    EpetraExt::RowMatrixToMatlabFile(filename_s2.str().c_str(), *Sff_);
  }
}


} // namespace
} // namespace
