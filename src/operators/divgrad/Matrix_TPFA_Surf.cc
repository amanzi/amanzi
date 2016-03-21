/*
  License: BSD
  Ethan Coon (ecoon@lanl.gov)
  Daniil Svyatskiy (dasvyat@lanl.gov)

  Matrix_TPFA_Surf provides for solving the cell-centered system with lambdas
  on the boundary system where a subset of
  the lambdas are more densely coupled for overland flow.

*/
#include "Epetra_FECrsGraph.h"
#include "EpetraExt_RowMatrixOut.h"

#include "errors.hh"
#include "MatrixMFD.hh"
#include "Matrix_TPFA.hh"
#include "Matrix_TPFA_Surf.hh"

#define DEBUG_FLAG 0

namespace Amanzi {
namespace Operators {

   Matrix_TPFA_Surf::Matrix_TPFA_Surf(Teuchos::ParameterList& plist,
				      const Teuchos::RCP<const AmanziMesh::Mesh>& mesh):
     Matrix_TPFA(plist,mesh) {fill_graph = false;};


void Matrix_TPFA_Surf::FillMatrixGraphs_(const Teuchos::Ptr<Epetra_CrsGraph> cf_graph,
          const Teuchos::Ptr<Epetra_FECrsGraph> fbfb_graph) {
  // get the standard fill
  Matrix_TPFA::FillMatrixGraphs_(cf_graph, fbfb_graph);

  // additional face-to-face connections for the TPF on surface.
  const Epetra_Map& cmap = mesh_->cell_map(false);
  const Epetra_Map& fmap = mesh_->face_map(false);
  const Epetra_Map& fmap_wghost = mesh_->face_map(true);
  const Epetra_Map& fbmap = mesh_->exterior_face_map(false);
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
    ierr = fbfb_graph->InsertGlobalIndices(ncells_surf, equiv_face_GID,
            ncells_surf, equiv_face_GID);
    ASSERT(!ierr);
  }
  fill_graph = true;
}


// Assumes the Surface A was already assembled.

void Matrix_TPFA_Surf::AssembleSchur_() const{

  int ierr(0);
  // std::vector<MatrixBC> new_markers(bc_markers);

  int ncells_surf = surface_mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  // for (AmanziMesh::Entity_ID sc=0; sc!=ncells_surf; ++sc) {
  //   AmanziMesh::Entity_ID f = surface_mesh_->entity_get_parent(AmanziMesh::CELL,sc);
  //   new_markers[f] = MATRIX_BC_NULL;
  // }

  // // Call Assemble TPFA matrices + boundary lambda

  // Matrix_TPFA::ApplyBoundaryConditions(new_markers, bc_values);
  // Get the standard TPFA pieces.
  Matrix_TPFA::AssembleSchur_();

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
  //int ncells_surf = surface_mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
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


  // Deal with RHS
  // const Epetra_MultiVector& rhs_surf_cells =
  //     *surface_A_->rhs()->ViewComponent("cell",false);
  // const Epetra_MultiVector& rhs_faces =
  //     *rhs_->ViewComponent("face",false);

  // for (AmanziMesh::Entity_ID sc=0; sc!=ncells_surf; ++sc) {
  //   AmanziMesh::Entity_ID f = surface_mesh_->entity_get_parent(AmanziMesh::CELL,sc);
  //   rhs_faces[0][f] += rhs_surf_cells[0][sc];
  //}
}

void Matrix_TPFA_Surf::AssembleRHS_() const{

  Matrix_TPFA::AssembleRHS_();

  int ncells_surf = surface_mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  // Deal with RHS
  const Epetra_MultiVector& rhs_surf_cells =
      *surface_A_->rhs()->ViewComponent("cell",false);
  const Epetra_MultiVector& rhs_faces =
      *rhs_->ViewComponent("boundary_face",false);
  const Epetra_Map& fb_map = mesh_->exterior_face_map(false);
  const Epetra_Map& f_map = mesh_->face_map(false);
  
  for (AmanziMesh::Entity_ID sc=0; sc!=ncells_surf; ++sc) {
    AmanziMesh::Entity_ID f = surface_mesh_->entity_get_parent(AmanziMesh::CELL, sc);
    int f_gid  = f_map.GID(f);
    int fb_lid = fb_map.LID(f_gid);
    rhs_faces[0][fb_lid] += rhs_surf_cells[0][sc];
  }


}


void Matrix_TPFA_Surf::ApplyBoundaryConditions(
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

  Matrix_TPFA::ApplyBoundaryConditions(new_markers, bc_values, ADD_BC_FLUX);
}

void Matrix_TPFA_Surf::ComputeNegativeResidual(const CompositeVector& solution,
			     const Teuchos::Ptr<CompositeVector>& residual) const{

  if (!assembled_rhs_) 
    AssembleRHS_();
  if (!assembled_schur_) {
    AssembleSchur_();    
    UpdatePreconditioner_();
  }

  Apply(solution, (*residual));


  //std::cout<<"solut face\n"<<*solution.ViewComponent("boundary_face", false)<<"\n";
  // //std::cout<<"resid\n"<<*(*residual).ViewComponent("cell", false)<<"\n";
  // std::cout<<"resid\n"<<*(*residual).ViewComponent("cell", false)<<"\n";

  residual->Update(-1.0, *rhs_, 1.0);

  //std::cout<<"rhs\n"<<*(*rhs_).ViewComponent("boundary_face", false)<<"\n";
  //std::cout<<"resid\n"<<*(*residual).ViewComponent("boundary_face", false)<<"\n";
  //exit(0);

}



void Matrix_TPFA_Surf::ComputeSchurComplement(const std::vector<MatrixBC>& bc_markers,
        const std::vector<double>& bc_values) {
  std::vector<MatrixBC> new_markers(bc_markers);

  int ierr(0);
  int ncells_surf = surface_mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (AmanziMesh::Entity_ID sc=0; sc!=ncells_surf; ++sc) {
    AmanziMesh::Entity_ID f = surface_mesh_->entity_get_parent(AmanziMesh::CELL,sc);
    new_markers[f] = MATRIX_BC_NULL;
  }

  // Call Assemble TPFA matrices + boundary lambda

  Matrix_TPFA::ApplyBoundaryConditions(new_markers, bc_values);
  Matrix_TPFA::AssembleSchur_();

  // dump the schur complement
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
    }

    ierr = Aff_->SumIntoGlobalValues(frow_global, entries, values, indices_global);
    ASSERT(!ierr);

  }

  ierr = Aff_->GlobalAssemble();
  ASSERT(!ierr);

  delete[] indices;
  delete[] indices_global;
  delete[] values;

  // dump the schur complement
  // std::stringstream filename_s2;
  // filename_s2 << "schur_post_surf_" << 0 << ".txt";
  // EpetraExt::RowMatrixToMatlabFile(filename_s2.str().c_str(), *Sff_);
}

// int Matrix_TPFA_Surf::Apply(const CompositeVector& X,
// 			     CompositeVector& Y){

//   return Matrix_TPFA::Apply(X,Y);

// }

// int Matrix_TPFA_Surf::ApplyInverse(const CompositeVector& X,
// 			     CompositeVector& Y){

//   return Matrix_TPFA::ApplyInverse(X,Y);

// }


} // namespace
} // namespace
