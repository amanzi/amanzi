#include "EpetraExt_RowMatrixOut.h"

#include "MatrixMFD_Permafrost.hh"

namespace Amanzi {
namespace Operators {

MatrixMFD_Permafrost::MatrixMFD_Permafrost(Teuchos::ParameterList& plist,
        const Teuchos::RCP<const AmanziMesh::Mesh> mesh) :
    MatrixMFD_Coupled(plist,mesh),
    scaling_(1.)
{}

MatrixMFD_Permafrost::MatrixMFD_Permafrost(const MatrixMFD_Permafrost& other) :
    MatrixMFD_Coupled(other),
    surface_mesh_(other.surface_mesh_),
    surface_A_(other.surface_A_),
    surface_B_(other.surface_B_),
    scaling_(other.scaling_)
{}


void MatrixMFD_Permafrost::ComputeSchurComplement() {
  // Base ComputeSchurComplement() gets the standard face parts
  MatrixMFD_Coupled::ComputeSchurComplement();

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

  int subsurf_entries = 0;
  int *gsubsurfindices;
  gsubsurfindices = new int[9];
  double *gsubsurfvalues;
  gsubsurfvalues = new double[9];  

  Epetra_SerialDenseMatrix block(2,2);

  int ierr(0);
  
  // First, clobber the surface face entries that are frozen, but above zero.
  Epetra_Vector scaling(P2f2f_->RowMap());
  Epetra_Vector App_diag(App.RowMap());
  ierr = App.ExtractDiagonalCopy(App_diag);
  AMANZI_ASSERT(!ierr);
  scaling.PutScalar(1.);

  int ncells_surf = surface_mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  for (AmanziMesh::Entity_ID sc=0; sc!=ncells_surf; ++sc) {
    AmanziMesh::Entity_ID frow = surface_mesh_->entity_get_parent(AmanziMesh::CELL,sc);
    if ((*blockA_sc_->Krel_)[frow] == 0 && std::abs(App_diag[sc]) > 1.e-11) {
      // Krel == 0 --> frozen, App_diag[sc] > 0 --> p_surf > p_atm
      scaling[2*frow] = 0.;
    }
  }
  ierr = P2f2f_->LeftScale(scaling);
  AMANZI_ASSERT(!ierr);

  // Now, add in the contributions from App
  // Loop over surface cells (subsurface faces)
  for (AmanziMesh::Entity_ID sc=0; sc!=ncells_surf; ++sc) {
    // Access the row from the surfaces
    AmanziMesh::Entity_ID sc_global = surf_cmap_wghost.GID(sc);

    ierr = App.ExtractGlobalRowCopy(sc_global, 9, entriesA, valuesA, indicesA);
    AMANZI_ASSERT(!ierr);
    ierr = Bpp.ExtractGlobalRowCopy(sc_global, 9, entriesB, valuesB, indicesB);
    AMANZI_ASSERT(!ierr);

    // ensure consistency of the sparsity structure, this can likely be
    // removed eventually.
    AMANZI_ASSERT(entriesA == entriesB);
    for (int m=0; m!=entriesA; ++m) {
      AMANZI_ASSERT(indicesA[m] == indicesB[m]);
    }

    // Convert local cell numbers to domain's local face numbers
    AmanziMesh::Entity_ID frow = surface_mesh_->entity_get_parent(AmanziMesh::CELL,sc);
    AmanziMesh::Entity_ID frow_global = fmap_wghost.GID(frow);
    int diag = -1;
    for (int m=0; m!=entriesA; ++m) {
      indicesA[m] = surf_cmap_wghost.LID(indicesA[m]);
      if (sc == indicesA[m]) diag = m; // save the diagonal index
      indicesA[m] = surface_mesh_->entity_get_parent(AmanziMesh::CELL,indicesA[m]);
      indicesA[m] = fmap_wghost.GID(indicesA[m]);
    }
    
    // Adjust for Krel factor.
    if ((*blockA_sc_->Krel_)[frow] > 0.) {
      // If not frozen surface, Krel > 0, and we can use the scaled constraint.
      for (int m=0; m!=entriesA; ++m) {
        // SCALING: -- divide by rel perm to keep constraint equation consistent
        if (std::abs(valuesA[m]) > 0.) {
          valuesA[m] /= (*blockA_sc_->Krel_)[frow];
        }
      }
    } else {
      // One of two cases is possible for Krel == 0
      //  - p_surf > p_atm, and we want the full values (hence the Row scaling done previously)
      //  - p_surf < p_atm, and so App's row is identically zero
    }

    // Add the entries
    ierr = P2f2f_->BeginSumIntoGlobalValues(frow_global, entriesA, indicesA);
    AMANZI_ASSERT(!ierr);

    for (int m=0; m!=entriesA; ++m) {
      // if (indicesA[m] == 500) {
      //   std::cout << "Adding Value from TPFA on surface to subsurface." << std::endl;
      //   std::cout << "  val from A = " << valuesA[m] << std::endl;
      //   std::cout << "  val from B = " << valuesB[m] << std::endl;
      // }

      block(0,0) = valuesA[m];
      block(1,1) = valuesB[m];

      // if (frow_global == indicesA[m]) {
      //   block(0,1) = (*Ccc_surf_)[0][sc] * scaling_;
      //   block(1,0) = (*Dcc_surf_)[0][sc] * scaling_;
      // }

      ierr = P2f2f_->SubmitBlockEntry(block.A(), block.LDA(), block.M(), block.N());
      AMANZI_ASSERT(!ierr);
    }

    ierr = P2f2f_->EndSubmitEntries();
    AMANZI_ASSERT(!ierr);
  }

  ierr = P2f2f_->GlobalAssemble();
  AMANZI_ASSERT(!ierr);

  // DEBUG dump
  // std::stringstream filename_s;
  // filename_s << "coupled_schur_" << 0 << ".txt";
  // EpetraExt::RowMatrixToMatlabFile(filename_s.str().c_str(), *P2f2f_);

  delete[] indicesA;
  delete[] indicesB;
  delete[] valuesA;
  delete[] valuesB;
  delete[] gsubsurfindices;
  delete[] gsubsurfvalues;
}


} // namespace
} // namespace
