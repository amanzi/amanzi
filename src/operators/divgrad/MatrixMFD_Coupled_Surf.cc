#include "EpetraExt_RowMatrixOut.h"

#include "MatrixMFD_Coupled_Surf.hh"

namespace Amanzi {
namespace Operators {

MatrixMFD_Coupled_Surf::MatrixMFD_Coupled_Surf(Teuchos::ParameterList& plist,
        const Teuchos::RCP<const AmanziMesh::Mesh> mesh) :
    MatrixMFD_Coupled(plist,mesh) {}

MatrixMFD_Coupled_Surf::MatrixMFD_Coupled_Surf(const MatrixMFD_Coupled_Surf& other) :
    MatrixMFD_Coupled(other),
    surface_mesh_(other.surface_mesh_),
    surface_A_(other.surface_A_),
    surface_B_(other.surface_B_) {}


void MatrixMFD_Coupled_Surf::ComputeSchurComplement() {

  // Base ComputeSchurComplement() gets the standard face parts
  MatrixMFD_Coupled::ComputeSchurComplement();

  // now need to add in the surf parts, this code mimics that in
  // MatrixMFD_Surf::ComputeSchurComplement().

  // Add the TPFA on the surface parts from surface_A.
  const Epetra_Map& surf_cmap_wghost = surface_mesh_->cell_map(true);
  const Epetra_Map& fmap_wghost = mesh_->face_map(true);

  // Grab the surface operators' TPFA cell-cell matrices
  const Epetra_FECrsMatrix& App = *surface_A_->TPFA();
  const Epetra_FECrsMatrix& Bpp = *surface_B_->TPFA();

  int entriesA = 0;
  int entriesB = 0;

  double *valuesA, *valuesB;
  valuesA = new double[9];
  valuesB = new double[9];

  int *indicesA, *indicesB;
  indicesA = new int[9];
  indicesB = new int[9];

  Epetra_SerialDenseMatrix block(2,2);

  // DEBUG dump
  // std::stringstream filename_sp;
  // filename_sp << "coupled_schur_pre_" << 0 << ".txt";
  // EpetraExt::RowMatrixToMatlabFile(filename_sp.str().c_str(), *P2f2f_);

  // Loop over surface cells (subsurface faces)
  int ierr(0);
  int ncells_surf = surface_mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (AmanziMesh::Entity_ID sc=0; sc!=ncells_surf; ++sc) {
    // Access the row from the surfaces
    AmanziMesh::Entity_ID sc_global = surf_cmap_wghost.GID(sc);

    ierr = App.ExtractGlobalRowCopy(sc_global, 9, entriesA, valuesA, indicesA);
    ASSERT(!ierr);
    ierr = Bpp.ExtractGlobalRowCopy(sc_global, 9, entriesB, valuesB, indicesB);
    ASSERT(!ierr);

    // ensure consistency of the sparsity structure, this can likely be
    // removed eventually.
    ASSERT(entriesA == entriesB);
    for (int m=0; m!=entriesA; ++m) {
      ASSERT(indicesA[m] == indicesB[m]);
    }

    // Convert local cell numbers to domain's local face numbers
    AmanziMesh::Entity_ID frow = surface_mesh_->entity_get_parent(AmanziMesh::CELL,sc);
    AmanziMesh::Entity_ID frow_global = fmap_wghost.GID(frow);
    for (int m=0; m!=entriesA; ++m) {
      indicesA[m] = surf_cmap_wghost.LID(indicesA[m]);
      indicesA[m] = surface_mesh_->entity_get_parent(AmanziMesh::CELL,indicesA[m]);
      indicesA[m] = fmap_wghost.GID(indicesA[m]);
    }

    // Add the entries
    ierr = P2f2f_->BeginSumIntoGlobalValues(frow_global, entriesA, indicesA);
    ASSERT(!ierr);

    for (int m=0; m!=entriesA; ++m) {
      if (indicesA[m] == 500) {
        std::cout << "Adding Value from TPFA on surface to subsurface." << std::endl;
        std::cout << "  val from A = " << valuesA[m] << std::endl;
        std::cout << "  val from B = " << valuesB[m] << std::endl;
      }

      block(0,0) = valuesA[m];
      block(1,1) = valuesB[m];

      // if (frow_global == indicesA[m]) {
      //   block(0,1) = (*Ccc_surf_)[0][sc];
      //   block(1,0) = (*Dcc_surf_)[0][sc];
      // }

      ierr = P2f2f_->SubmitBlockEntry(block.A(), block.LDA(), block.M(), block.N());
      ASSERT(!ierr);
    }

    ierr = P2f2f_->EndSubmitEntries();
    ASSERT(!ierr);
  }

  ierr = P2f2f_->GlobalAssemble();
  ASSERT(!ierr);

  // DEBUG dump
  // std::stringstream filename_s;
  // filename_s << "coupled_schur_" << 0 << ".txt";
  // EpetraExt::RowMatrixToMatlabFile(filename_s.str().c_str(), *P2f2f_);

  delete[] indicesA;
  delete[] indicesB;
  delete[] valuesA;
  delete[] valuesB;
}


} // namespace
} // namespace
