/*
  License: BSD
  Ethan Coon (ecoon@lanl.gov)

  MatrixMFD_Surf provides for solving the p-lambda system where a subset of
  the lambdas are more densely coupled for overland flow.

*/
#include "Epetra_FECrsGraph.h"

#include "errors.hh"
#include "matrix_mfd_surf.hh"

#define DEBUG_FLAG 0

namespace Amanzi {
namespace Operators {


MatrixMFD_Surf::MatrixMFD_Surf(Teuchos::ParameterList& plist,
        const Teuchos::RCP<const AmanziMesh::Mesh> mesh,
        const Teuchos::RCP<const AmanziMesh::Mesh> surface_mesh) :
    MatrixMFD(plist,mesh),
    surface_mesh_(surface_mesh) {}


MatrixMFD_Surf::MatrixMFD_Surf(const MatrixMFD& other,
        const Teuchos::RCP<const AmanziMesh::Mesh> surface_mesh) :
    MatrixMFD(other),
    surface_mesh_(surface_mesh) {}

void MatrixMFD_Surf::SymbolicAssembleGlobalMatrices() {
  const Epetra_Map& cmap = mesh_->cell_map(false);
  const Epetra_Map& fmap = mesh_->face_map(false);
  const Epetra_Map& fmap_wghost = mesh_->face_map(true);
  int ierr(0);

  int avg_entries_row = (mesh_->space_dimension() == 2) ? MFD_QUAD_FACES : MFD_HEX_FACES;
  Epetra_CrsGraph cf_graph(Copy, cmap, fmap_wghost, avg_entries_row, false);  // FIX (lipnikov@lanl.gov)
  Epetra_FECrsGraph ff_graph(Copy, fmap, 2*avg_entries_row);

  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;
  int faces_LID[MFD_MAX_FACES];  // Contigious memory is required.
  int faces_GID[MFD_MAX_FACES];

  // deals with the standard MFD connections
  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c=0; c!=ncells; ++c) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    for (int n=0; n!=nfaces; ++n) {
      faces_LID[n] = faces[n];
      faces_GID[n] = fmap_wghost.GID(faces_LID[n]);
    }
    ierr = cf_graph.InsertMyIndices(c, nfaces, faces_LID);
    ASSERT(!ierr);
    ierr = ff_graph.InsertGlobalIndices(nfaces, faces_GID, nfaces, faces_GID);
    ASSERT(!ierr);
  }

  // additional face-to-face connections for the TPF on surface.
  AmanziMesh::Entity_ID_List surf_cells;
  int equiv_face_GID[2];

  int nfaces_surf = surface_mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);

  if (mesh_->get_comm()->MyPID() != 0) mesh_->get_comm()->Barrier();

  for (int fs=0; fs!=nfaces_surf; ++fs) {
    surface_mesh_->face_get_cells(fs, AmanziMesh::USED, &surf_cells);
    int ncells_surf = surf_cells.size(); ASSERT(ncells_surf <= 2);

    // get the equivalent faces on the subsurface mesh
    for (int n=0; n!=ncells_surf; ++n) {
      surf_cells[n] = surface_mesh_->entity_get_parent(AmanziMesh::CELL, surf_cells[n]);
      equiv_face_GID[n] = fmap_wghost.GID(surf_cells[n]);
    }

    // insert the connection
    ierr = ff_graph.InsertGlobalIndices(ncells_surf, equiv_face_GID,
            ncells_surf, equiv_face_GID);
    ASSERT(!ierr);
  }

  if (mesh_->get_comm()->MyPID() == 0) mesh_->get_comm()->Barrier();

  ierr = cf_graph.FillComplete(fmap, cmap);
  ASSERT(!ierr);

  ierr = ff_graph.GlobalAssemble();  // Symbolic graph is complete.
  ASSERT(!ierr);

  // create global matrices
  Acc_ = Teuchos::rcp(new Epetra_Vector(cmap));
  Acf_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, cf_graph));
  Aff_ = Teuchos::rcp(new Epetra_FECrsMatrix(Copy, ff_graph));
  Sff_ = Teuchos::rcp(new Epetra_FECrsMatrix(Copy, ff_graph));
  Aff_->GlobalAssemble();
  Sff_->GlobalAssemble();

  // Symmetry
  if (flag_symmetry_) {
    Afc_ = Acf_;
  } else {
    Afc_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, cf_graph));
  }

  std::vector<std::string> names(2);
  names[0] = "cell";
  names[1] = "face";

  std::vector<AmanziMesh::Entity_kind> locations(2);
  locations[0] = AmanziMesh::CELL;
  locations[1] = AmanziMesh::FACE;

  std::vector<int> num_dofs(2,1);
  rhs_ = Teuchos::rcp(new CompositeVector(mesh_, names, locations, num_dofs, true));
  rhs_->CreateData();
}

// Assumes the Surface A was already assembled.
void MatrixMFD_Surf::AssembleGlobalMatrices() {
  int ierr(0);

  // Get the standard MFD pieces.
  MatrixMFD::AssembleGlobalMatrices();

  // Add the TPFA on the surface parts from surface_A.
  const Epetra_Map& fmap_wghost = mesh_->face_map(true);
  const Epetra_FECrsMatrix& Spp = *surface_A_->TPFA();

  int entries = 0;
  double *values;
  values = new double[9];
  int *indices;
  indices = new int[9];

  // Loop over surface cells (subsurface faces)
  int ncells_surf = surface_mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (AmanziMesh::Entity_ID sc=0; sc!=ncells_surf; ++sc) {
    // Access the row in Spp
    ierr = Spp.ExtractMyRowCopy(sc, 9, entries, values, indices);
    ASSERT(!ierr);

    // Convert Spp local cell numbers to Aff local face numbers
    AmanziMesh::Entity_ID frow = surface_mesh_->entity_get_parent(AmanziMesh::CELL,sc);
    int frow_global = fmap_wghost.GID(frow);

    for (int m=0; m!=entries; ++m) {
      indices[m] = surface_mesh_->entity_get_parent(AmanziMesh::CELL,indices[m]);
      indices[m] = fmap_wghost.GID(indices[m]);
    }

    ierr = Aff_->SumIntoGlobalValues(frow_global, entries, values, indices);
    ASSERT(!ierr);
  }

  ierr = Aff_->GlobalAssemble();
  ASSERT(!ierr);

  delete[] indices;
  delete[] values;


  // Deal with RHS
  const Epetra_MultiVector& rhs_surf_cells =
      *surface_A_->rhs()->ViewComponent("cell",false);
  const Epetra_MultiVector& rhs_faces =
      *rhs_->ViewComponent("face",false);
  for (AmanziMesh::Entity_ID sc=0; sc!=ncells_surf; ++sc) {
    AmanziMesh::Entity_ID f = surface_mesh_->entity_get_parent(AmanziMesh::CELL,sc);
    rhs_faces[0][f] += rhs_surf_cells[0][sc];
  }
}


void MatrixMFD_Surf::ApplyBoundaryConditions(
    const std::vector<Matrix_bc>& bc_markers,
    const std::vector<double>& bc_values) {
  // Ensure that none of the surface faces have a BC in them.
  std::vector<Matrix_bc> new_markers(bc_markers);

  int ncells_surf = surface_mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (AmanziMesh::Entity_ID sc=0; sc!=ncells_surf; ++sc) {
    AmanziMesh::Entity_ID f = surface_mesh_->entity_get_parent(AmanziMesh::CELL,sc);
    new_markers[f] = MATRIX_BC_NULL;
  }

  MatrixMFD::ApplyBoundaryConditions(new_markers, bc_values);
}



void MatrixMFD_Surf::ComputeSchurComplement(const std::vector<Matrix_bc>& bc_markers,
        const std::vector<double>& bc_values) {
  std::vector<Matrix_bc> new_markers(bc_markers);

  int ierr(0);
  int ncells_surf = surface_mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (AmanziMesh::Entity_ID sc=0; sc!=ncells_surf; ++sc) {
    AmanziMesh::Entity_ID f = surface_mesh_->entity_get_parent(AmanziMesh::CELL,sc);
    new_markers[f] = MATRIX_BC_NULL;
  }

  // Call base Schur
  MatrixMFD::ComputeSchurComplement(new_markers, bc_values);

  // Add the TPFA on the surface parts from surface_A.
  const Epetra_Map& fmap_wghost = mesh_->face_map(true);

  // ASSUMES surface_A_'s AssembleGlobalMatrices() has been called
  const Epetra_FECrsMatrix& Spp = *surface_A_->TPFA();

  int entries = 0;
  double *values;
  values = new double[9];
  int *indices;
  indices = new int[9];

  // Loop over surface cells (subsurface faces)
  for (AmanziMesh::Entity_ID sc=0; sc!=ncells_surf; ++sc) {
    // Access the row in Spp
    ierr = Spp.ExtractMyRowCopy(sc, 9, entries, values, indices);
    ASSERT(!ierr);

    // Convert Spp local cell numbers to Sff local face numbers
    AmanziMesh::Entity_ID frow = surface_mesh_->entity_get_parent(AmanziMesh::CELL,sc);
    AmanziMesh::Entity_ID frow_global = fmap_wghost.GID(frow);
    for (int m=0; m!=entries; ++m) {
      indices[m] = surface_mesh_->entity_get_parent(AmanziMesh::CELL,indices[m]);
      indices[m] = fmap_wghost.GID(indices[m]);
    }

    ierr = Sff_->SumIntoGlobalValues(frow_global, entries, values, indices);
    ASSERT(!ierr);

  }

  ierr = Sff_->GlobalAssemble();
  ASSERT(!ierr);

  delete[] indices;
  delete[] values;


}


} // namespace
} // namespace
