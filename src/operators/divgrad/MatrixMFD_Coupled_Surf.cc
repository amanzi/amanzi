#include "EpetraExt_RowMatrixOut.h"
#include "Epetra_FECrsGraph.h"

#include "MatrixMFD_Coupled_Surf.hh"
#include "MatrixMFD_Surf.hh"

namespace Amanzi {
namespace Operators {

MatrixMFD_Coupled_Surf::MatrixMFD_Coupled_Surf(Teuchos::ParameterList& plist,
        const Teuchos::RCP<const AmanziMesh::Mesh> mesh) :
    MatrixMFD_Coupled(plist,mesh)
{}

MatrixMFD_Coupled_Surf::MatrixMFD_Coupled_Surf(const MatrixMFD_Coupled_Surf& other) :
    MatrixMFD_Coupled(other),
    surface_mesh_(other.surface_mesh_),
    surface_A_(other.surface_A_),
    surface_B_(other.surface_B_)
{}


void
MatrixMFD_Coupled_Surf::SetSurfaceOperators(const Teuchos::RCP<MatrixMFD_TPFA>& surface_A,
                         const Teuchos::RCP<MatrixMFD_TPFA>& surface_B) {
  surface_A_ = surface_A;
  surface_B_ = surface_B;
  ASSERT(surface_A_->Mesh() == surface_B_->Mesh());
  surface_mesh_ = surface_A_->Mesh();

  // Create the surface->subsurface importer and map
  int nsurf_cells = surface_mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  std::vector<int> surf_gids(nsurf_cells, -1);

  const Epetra_Map& face_map = blockA_->Mesh()->face_map(false);
  for (unsigned int sc=0; sc!=nsurf_cells; ++sc) {
    surf_gids[sc] = face_map.GID(surface_mesh_->entity_get_parent(AmanziMesh::CELL, sc));
  }

  surf_map_in_subsurf_ = Teuchos::rcp(new Epetra_Map(-1, nsurf_cells, &surf_gids[0], 0, *blockA_->Mesh()->get_comm()));
  surf_importer_ = Teuchos::rcp(new Epetra_Import(*surf_map_in_subsurf_,face_map));
}

void MatrixMFD_Coupled_Surf::SymbolicAssembleGlobalMatrices() {
  // Identical to MatrixMFD_Coupled, but does the FillMatrixGraph call
  // with a _Surf matrix, which must first be created.
  int ierr(0);
  const Epetra_BlockMap& cmap = mesh_->cell_map(false);
  const Epetra_BlockMap& fmap = mesh_->face_map(false);
  const Epetra_BlockMap& fmap_wghost = mesh_->face_map(true);

  // Make the double maps
  double_fmap_ = Teuchos::rcp(new Epetra_BlockMap(fmap.NumGlobalPoints(),
          fmap.NumMyPoints(), fmap.MyGlobalElements(), 2,
          fmap.IndexBase(), fmap.Comm() ));

  double_cmap_ = Teuchos::rcp(new Epetra_BlockMap(cmap.NumGlobalPoints(),
          cmap.NumMyPoints(), cmap.MyGlobalElements(), 2,
          cmap.IndexBase(), cmap.Comm() ));

  double_fmap_wghost_ = Teuchos::rcp(
      new Epetra_BlockMap(fmap_wghost.NumGlobalPoints(),
                          fmap_wghost.NumMyPoints(),
                          fmap_wghost.MyGlobalElements(), 2,
                          fmap_wghost.IndexBase(), fmap_wghost.Comm() ));

  // Make the matrix graphs
  int avg_entries_row = 6;
  Teuchos::RCP<Epetra_CrsGraph> cf_graph =
      Teuchos::rcp(new Epetra_CrsGraph(Copy, *double_cmap_, *double_fmap_wghost_,
              avg_entries_row, false));
  Teuchos::RCP<Epetra_FECrsGraph> ff_graph =
      Teuchos::rcp(new Epetra_FECrsGraph(Copy, *double_fmap_, 2*avg_entries_row - 1,
              false));

  // Fill the graphs from the symbolic patterns in matrix_mfd from one of the
  // sub-block matrices.

  // IT IS ASSUMED that the two sub-block matrices fill in identical ways!
  // Create a MatrixMFD_Surf operator for filling.
  MatrixMFD_Surf blockA_surf(blockA_->plist_, blockA_->mesh_);
  blockA_surf.SetSurfaceOperator(surface_A_);
  blockA_surf.FillMatrixGraphs_(cf_graph.ptr(), ff_graph.ptr());

  // Assemble the graphs
  ierr = cf_graph->FillComplete(*double_fmap_, *double_cmap_);
  ASSERT(!ierr);
  ierr = ff_graph->GlobalAssemble();  // Symbolic graph is complete.
  ASSERT(!ierr);

  // Create the matrices
  A2f2c_ = Teuchos::rcp(new Epetra_VbrMatrix(Copy, *cf_graph)); // stored in transpose
  A2c2f_ = Teuchos::rcp(new Epetra_VbrMatrix(Copy, *cf_graph));
  P2f2f_ = Teuchos::rcp(new Epetra_FEVbrMatrix(Copy, *ff_graph, false));

  P2f2f_ = Teuchos::rcp(new Epetra_FEVbrMatrix(Copy, *ff_graph, false));
  ierr = P2f2f_->GlobalAssemble();

  if (assemble_matrix_) {
    A2f2f_ = Teuchos::rcp(new Epetra_FEVbrMatrix(Copy, *ff_graph, false));
    ierr = A2f2f_->GlobalAssemble();
  }
  ASSERT(!ierr);
}


void MatrixMFD_Coupled_Surf::AssembleGlobalMatrices() {
  // Assemble Aff, without Surf
  ASSERT(assemble_matrix_);
  MatrixMFD_Coupled::AssembleGlobalMatrices();

  // Add in surf terms
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

  // Now, add in the contributions from App
  // Loop over surface cells (subsurface faces)
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
    int diag = -1;
    for (int m=0; m!=entriesA; ++m) {
      indicesA[m] = surf_cmap_wghost.LID(indicesA[m]);
      if (sc == indicesA[m]) diag = m; // save the diagonal index
      indicesA[m] = surface_mesh_->entity_get_parent(AmanziMesh::CELL,indicesA[m]);
      indicesA[m] = fmap_wghost.GID(indicesA[m]);
    }

    /*
    if (sc == 65) {
      std::cout << "MatrixMFD_Coupled_Surf sc(65), surf terms:\n"
		<< "     scaling=" << scaling_ << "\n"
		<< "  App_surf[65,65] = " << valuesA[diag] << "\n"
		<< "  Bpp_surf[65,65] = " << valuesB[diag] << "\n"
		<< "  Ccc_surf[65,65] = " << (*Ccc_surf_)[0][sc]*scaling_ << "\n"
		<< "  Dcc_surf[65,65] = " << (*Dcc_surf_)[0][sc]*scaling_ << std::endl;
    }
    */

    // Add the entries
    ierr = A2f2f_->BeginSumIntoGlobalValues(frow_global, entriesA, indicesA);
    ASSERT(!ierr);


    for (int m=0; m!=entriesA; ++m) {
      block(0,0) = valuesA[m];
      //      block(0,0) = std::max(valuesA[m], 1.e-10);
      block(1,1) = valuesB[m];

      if (frow_global == indicesA[m]) {
        block(0,1) = (*Ccc_surf_)[0][sc] * scaling_;
        block(1,0) = (*Dcc_surf_)[0][sc] * scaling_;
      } else {
        block(0,1) = 0.;
        block(1,0) = 0.;
      }

      ierr = A2f2f_->SubmitBlockEntry(block.A(), block.LDA(), block.M(), block.N());
      ASSERT(!ierr);
    }

    ierr = A2f2f_->EndSubmitEntries();
    ASSERT(!ierr);
  }

  ierr = A2f2f_->GlobalAssemble();
  ASSERT(!ierr);

  delete[] indicesA;
  delete[] indicesB;
  delete[] valuesA;
  delete[] valuesB;
  delete[] gsubsurfindices;
  delete[] gsubsurfvalues;
  
}


void MatrixMFD_Coupled_Surf::ComputeSchurComplement() {
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

  // Now, add in the contributions from App
  // Loop over surface cells (subsurface faces)
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
    int diag = -1;
    for (int m=0; m!=entriesA; ++m) {
      indicesA[m] = surf_cmap_wghost.LID(indicesA[m]);
      if (sc == indicesA[m]) diag = m; // save the diagonal index
      indicesA[m] = surface_mesh_->entity_get_parent(AmanziMesh::CELL,indicesA[m]);
      indicesA[m] = fmap_wghost.GID(indicesA[m]);
    }

    /*
    if (sc == 65) {
      std::cout << "MatrixMFD_Coupled_Surf sc(65), surf terms:\n"
		<< "     scaling=" << scaling_ << "\n"
		<< "  App_surf[65,65] = " << valuesA[diag] << "\n"
		<< "  Bpp_surf[65,65] = " << valuesB[diag] << "\n"
		<< "  Ccc_surf[65,65] = " << (*Ccc_surf_)[0][sc]*scaling_ << "\n"
		<< "  Dcc_surf[65,65] = " << (*Dcc_surf_)[0][sc]*scaling_ << std::endl;
    }
    */

    // Add the entries
    ierr = P2f2f_->BeginSumIntoGlobalValues(frow_global, entriesA, indicesA);
    ASSERT(!ierr);


    for (int m=0; m!=entriesA; ++m) {
      block(0,0) = valuesA[m];
      //      block(0,0) = std::max(valuesA[m], 1.e-10);
      block(1,1) = valuesB[m];

      if (frow_global == indicesA[m]) {
        block(0,1) = (*Ccc_surf_)[0][sc] * scaling_;
        block(1,0) = (*Dcc_surf_)[0][sc] * scaling_;
      } else {
        block(0,1) = 0.;
        block(1,0) = 0.;
      }

      ierr = P2f2f_->SubmitBlockEntry(block.A(), block.LDA(), block.M(), block.N());
      ASSERT(!ierr);
    }

    ierr = P2f2f_->EndSubmitEntries();
    ASSERT(!ierr);
  }

  ierr = P2f2f_->GlobalAssemble();
  ASSERT(!ierr);

  // DEBUG dump
  if (dump_schur_) {
    std::stringstream filename_s;
    filename_s << "schur_MatrixMFD_Coupled_Surf_" << 0 << ".txt";
    EpetraExt::RowMatrixToMatlabFile(filename_s.str().c_str(), *P2f2f_);
  }

  delete[] indicesA;
  delete[] indicesB;
  delete[] valuesA;
  delete[] valuesB;
  delete[] gsubsurfindices;
  delete[] gsubsurfvalues;
}


int
MatrixMFD_Coupled_Surf::Apply(const TreeVector& X,
        TreeVector& Y) const {
  // Apply the coupled Matrix, which gets the subsurface.
  int ierr = MatrixMFD_Coupled::Apply(X,Y);

  // Add in the Surface components
  Teuchos::RCP<const CompositeVector> XA = X.SubVector(0)->Data();
  Teuchos::RCP<const CompositeVector> XB = X.SubVector(1)->Data();
  Teuchos::RCP<CompositeVector> YA = Y.SubVector(0)->Data();
  Teuchos::RCP<CompositeVector> YB = Y.SubVector(1)->Data();

  // -- import X from subsurface faces into surface cells

  // NOTE --etc
  // This confuses the hell out of me.  Is it not possible to pull out an
  // arbitrary set of global indices into a vector from a previously created
  // map?  The obvious form of this errors because the target vector's map
  // (the surface cell map) is not the same as the importer's target map (the
  // subset of face GIDs).  To me that should be ok anyway... the nth face GID
  // corresponds to the nth surface cell.
  //
  // Note that the target MUST have the surface cell map, as it is used with
  // Matrices that use that map.
  //
  // Two potential options forward:
  // -- two imports are required (the send taking the face GID referenced data
  //    and copying it to the cell Vector), and the second one is a simple
  //    copy with no re-mapping.
  // -- fool Trilinos by making a new Epetra_Vector for the surface cell
  //    vector that shares data with the old one but has a new map.
  // We try the latter as it removes the copy...
  //
  Epetra_MultiVector surf_XA(surface_mesh_->cell_map(false),1);
  Epetra_MultiVector surf_XB(surface_mesh_->cell_map(false),1);

  double **surf_XA_data, **surf_XB_data;
  ierr |= surf_XA.ExtractView(&surf_XA_data);
  ierr |= surf_XB.ExtractView(&surf_XB_data);
  ASSERT(!ierr);

  Epetra_MultiVector surf_XA_face_GIDs(View, *surf_map_in_subsurf_, surf_XA_data, 1);
  Epetra_MultiVector surf_XB_face_GIDs(View, *surf_map_in_subsurf_, surf_XB_data, 1);

  // Now, surf_XA_face_GIDs and surf_XA have the same values, but different maps.
  // Importing into the face_GIDs version results in the values in surf_XA
  ierr |= surf_XA_face_GIDs.Import(*XA->ViewComponent("face",false), *surf_importer_, Insert);
  ierr |= surf_XB_face_GIDs.Import(*XB->ViewComponent("face",false), *surf_importer_, Insert);
  ASSERT(!ierr);

  // -- apply the surface-only operators, blockwise, repeating the "multiple views" trick
  Epetra_MultiVector surf_YA(surface_mesh_->cell_map(false),1);
  Epetra_MultiVector surf_YB(surface_mesh_->cell_map(false),1);

  double **surf_YA_data, **surf_YB_data;
  ierr |= surf_YA.ExtractView(&surf_YA_data);
  ierr |= surf_YB.ExtractView(&surf_YB_data);
  ASSERT(!ierr);

  Epetra_MultiVector surf_YA_face_GIDs(View, *surf_map_in_subsurf_, surf_YA_data, 1);
  Epetra_MultiVector surf_YB_face_GIDs(View, *surf_map_in_subsurf_, surf_YB_data, 1);

  // -- A_surf * lambda_p
  ierr |= surface_A_->Apply(surf_XA, surf_YA);

  // -- Ccc_surf * lambda_T
  ierr |= surf_YA.Multiply(scaling_, *Ccc_surf_, surf_XB, 1.0);

  // -- B_surf * lambda_T
  ierr |= surface_B_->Apply(surf_XB, surf_YB);

  // -- Dcc_surf * lambda_p
  ierr |= surf_YB.Multiply(scaling_, *Dcc_surf_, surf_XA, 1.0);
  ASSERT(!ierr);
  
  // -- export back into subsurface face vectors
  ierr |= YA->ViewComponent("face",false)->Export(surf_YA_face_GIDs, *surf_importer_, Add);
  ierr |= YB->ViewComponent("face",false)->Export(surf_YB_face_GIDs, *surf_importer_, Add);
  ASSERT(!ierr);

  return ierr;
}


void MatrixMFD_Coupled_Surf::UpdateConsistentFaceCorrection(const TreeVector& u,
        const Teuchos::Ptr<TreeVector>& Pu) {
  ASSERT(assemble_matrix_);

  // Aff solutions
  if (Aff_solver_ == Teuchos::null) {
    if (plist_.isSublist("consistent face solver")) {
      Teuchos::ParameterList Aff_plist = plist_.sublist("consistent face solver");
      Aff_op_ = Teuchos::rcp(new EpetraMatrixDefault<Epetra_FEVbrMatrix>(Aff_plist));
      Aff_op_->Update(A2f2f_);

      if (Aff_plist.isParameter("iterative method")) {
        AmanziSolvers::LinearOperatorFactory<EpetraMatrix,Epetra_Vector,Epetra_BlockMap> op_fac;
        Aff_solver_ = op_fac.Create(Aff_plist, Aff_op_);
      } else {
        Aff_solver_ = Aff_op_;
      }
    } else {
      Errors::Message msg("MatrixMFD::UpdateConsistentFaceConstraints was called, but no consistent face solver sublist was provided.");
      Exceptions::amanzi_throw(msg);
    }
  }


  // pull cell data
  Teuchos::RCP<CompositeVector> PuA = Pu->SubVector(0)->Data();
  Teuchos::RCP<CompositeVector> PuB = Pu->SubVector(1)->Data();
  const Epetra_MultiVector& PuA_c = *PuA->ViewComponent("cell", false);
  const Epetra_MultiVector& PuB_c = *PuB->ViewComponent("cell", false);
  Epetra_MultiVector& PuA_f = *PuA->ViewComponent("face", false);
  Epetra_MultiVector& PuB_f = *PuB->ViewComponent("face", false);

  Teuchos::RCP<const CompositeVector> uA = Pu->SubVector(0)->Data();
  Teuchos::RCP<const CompositeVector> uB = Pu->SubVector(1)->Data();
  const Epetra_MultiVector& uA_f = *uA->ViewComponent("face", false);
  const Epetra_MultiVector& uB_f = *uB->ViewComponent("face", false);

  // Temporary cell and face vectors.
  Epetra_MultiVector Pu_c(*double_cmap_, 1);
  Epetra_MultiVector update_f(*double_fmap_, 1);
  Epetra_MultiVector Pu_f(*double_fmap_, 1);

  int ierr(0);

  int ncells = PuA_c.MyLength();
  int nfaces = PuA_f.MyLength();
  double val = 0.;

  // update_f <-- -A2f2c * Pu_c
  for (int c=0; c!=ncells; ++c){
    Pu_c[0][2*c] = -PuA_c[0][c];
    Pu_c[0][2*c+1] = -PuB_c[0][c];
  }
  ierr = A2f2c_->Multiply(true, Pu_c, update_f);

  // update_f <-- u_f - A2f2c * Pu_c
  for (int f=0; f!=nfaces; ++f){
    update_f[0][2*f] += uA_f[0][f];
    update_f[0][2*f+1] += uB_f[0][f];
  }

  // u_f <-- A2f2f ^ -1 * update_f
  Aff_op_->Destroy();
  Aff_op_->Update(A2f2f_);
  ierr = Aff_solver_->ApplyInverse(*update_f(0), *Pu_f(0));
  ASSERT(!ierr);

  // put back into tree vectors
  for (int f=0; f!=nfaces; ++f){
    PuA_f[0][f] = Pu_f[0][2*f];
    PuB_f[0][f] = Pu_f[0][2*f+1];
  }
}


} // namespace
} // namespace
