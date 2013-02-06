/*
  License: BSD
  Ethan Coon (ecoon@lanl.gov)

  MatrixMFD_Surf provides for solving the p-lambda system where a subset of
  the lambdas are more densely coupled for overland flow.

*/
#include "Epetra_FECrsGraph.h"

#include "errors.hh"
#include "matrix_mfd_surf.hh"


namespace Amanzi {
namespace Operators {


MatrixMFD_Surf::MatrixMFD_Surf(Teuchos::ParameterList& plist,
        const Teuchos::RCP<const AmanziMesh::Mesh> mesh,
        const Teuchos::RCP<const AmanziMesh::Mesh> surface_mesh) :
    MatrixMFD(plist,mesh),
    surface_mesh_(surface_mesh) {}

void MatrixMFD_Surf::SymbolicAssembleGlobalMatrices() {
  const Epetra_Map& cmap = mesh_->cell_map(false);
  const Epetra_Map& fmap = mesh_->face_map(false);
  const Epetra_Map& fmap_wghost = mesh_->face_map(true);

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
    cf_graph.InsertMyIndices(c, nfaces, faces_LID);
    ff_graph.InsertGlobalIndices(nfaces, faces_GID, nfaces, faces_GID);
  }

  // additional face-to-face connections for the TPF on surface.
  AmanziMesh::Entity_ID_List surf_cells;
  int equiv_face_GID[2];

  int nfaces_surf = surface_mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  for (int fs=0; fs!=nfaces_surf; ++fs) {
    surface_mesh_->face_get_cells(fs, AmanziMesh::USED, &surf_cells);

    // get the equivalent faces on the subsurface mesh
    int ncells_surf = surf_cells.size(); ASSERT(ncells_surf <= 2);
    for (int n=0; n!=ncells_surf; ++n) {
      int f = surface_mesh_->entity_get_parent(AmanziMesh::CELL, surf_cells[n]);
      equiv_face_GID[n] = fmap_wghost.GID(f);
    }

    // insert the connection
    ff_graph.InsertGlobalIndices(ncells_surf, equiv_face_GID,
            ncells_surf, equiv_face_GID);
  }

  cf_graph.FillComplete(fmap, cmap);
  ff_graph.GlobalAssemble();  // Symbolic graph is complete.

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


void MatrixMFD_Surf::AssembleGlobalMatrices() {
  // Get the standard MFD pieces.
  MatrixMFD::AssembleGlobalMatrices();
  surface_A_->AssembleGlobalMatrices();

  // Add the TPFA on the surface parts from surface_A.
  const Epetra_Map& fmap_wghost = mesh_->face_map(true);
  const Epetra_FECrsMatrix& Spp = *surface_A_->TPFA();

  int entries = 0;
  double *values;
  int *indices;

  // Loop over surface cells (subsurface faces)
  int ncells_surf = surface_mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (AmanziMesh::Entity_ID sc=0; sc!=ncells_surf; ++sc) {
    // Access the row in Spp
    int ierr = Spp.ExtractMyRowView(sc, entries, values, indices);

    // Convert Spp local cell numbers to Aff local face numbers
    AmanziMesh::Entity_ID frow = surface_mesh_->entity_get_parent(AmanziMesh::CELL,sc);
    for (int m=0; m!=entries; ++m) {
      indices[m] = surface_mesh_->entity_get_parent(AmanziMesh::CELL,indices[m]);

    }

    // Add into Aff
    Aff_->SumIntoMyValues(frow, entries, values, indices);
  }
  Aff_->GlobalAssemble();

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
    const std::vector<Matrix_bc>& subsurface_markers,
    const std::vector<double>& subsurface_values,
    const std::vector<Matrix_bc>& surface_markers,
    const std::vector<double>& surface_values) {
  // Ensure that none of the surface faces have a BC in them.
  int ncells_surf = surface_mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (AmanziMesh::Entity_ID sc=0; sc!=ncells_surf; ++sc) {
    AmanziMesh::Entity_ID f = surface_mesh_->entity_get_parent(AmanziMesh::CELL,sc);
    if (subsurface_markers[f] != MFD_BC_NULL) {
      Errors::Message msg("MatrixMFD_Surf::Subsurface's cell on the surface has a non-Null BC.");
      Exceptions::amanzi_throw(msg);
    }
  }

  MatrixMFD::ApplyBoundaryConditions(subsurface_markers, subsurface_values);
  surface_A_->ApplyBoundaryConditions(surface_markers, surface_values);
}



void MatrixMFD_Surf::ComputeSchurComplement(const std::vector<Matrix_bc>& bc_markers,
        const std::vector<double>& bc_values) {
  // Ensure that none of the surface faces have a BC in them.
  int ncells_surf = surface_mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (AmanziMesh::Entity_ID sc=0; sc!=ncells_surf; ++sc) {
    AmanziMesh::Entity_ID f = surface_mesh_->entity_get_parent(AmanziMesh::CELL,sc);
    if (subsurface_markers[f] != MFD_BC_NULL) {
      Errors::Message msg("MatrixMFD_Surf::Subsurface's cell on the surface has a non-Null BC.");
      Exceptions::amanzi_throw(msg);
    }
  }

  // Call base Schur
  MatrixMFD::ComputeSchurComplement(bc_markers, bc_values);
}


} // namespace
} // namespace
