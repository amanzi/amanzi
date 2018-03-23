/*
  This is the flow component of the Amanzi code.

  Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
  Daniil Svyatskiy (dasvyat@lanl.gov)
  Ethan Coon (ATS version) (ecoon@lanl.gov)
*/

#include <vector>

#include "Epetra_FECrsGraph.h"

#include "errors.hh"
#include "Matrix_TPFA.hh"
#include "BlockMatrix.hh"
#include "LinearOperatorFactory.hh"



namespace Amanzi {
namespace Operators {

Matrix_TPFA::Matrix_TPFA(Teuchos::ParameterList& plist,
                         const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) :
    MatrixMFD(plist,mesh) {

  // create the space
  cells_only_ = plist.get<bool>("TPFA use cells only", false);
  if (cells_only_) {
    space_ = Teuchos::rcp(new CompositeVectorSpace());
    space_->SetMesh(mesh_)->SetComponent("cell", AmanziMesh::CELL, 1);

  } else {
    std::vector<std::string> names(2);
    names[0] = "cell";
    names[1] = "boundary_face";
    std::vector<AmanziMesh::Entity_kind> locations(2);
    locations[0] = AmanziMesh::CELL;
    locations[1] = AmanziMesh::BOUNDARY_FACE;

    std::vector<int> num_dofs(2,1);
    space_ = Teuchos::rcp(new CompositeVectorSpace());
    space_->SetMesh(mesh_)->SetGhosted()->SetComponents(names,locations,num_dofs);
  }
}


/* ******************************************************************
 * Update face transmisibilities.
 ****************************************************************** */
void Matrix_TPFA::CreateMFDstiffnessMatrices(
    const Teuchos::Ptr<const CompositeVector>& Krel) {
  // tag global matrices as invalid
  assembled_schur_ = false;
  assembled_operator_ = false;
  assembled_rhs_ = false;

  // communicate as necessary
  if (Krel.get() && Krel->HasComponent("face")) Krel->ScatterMasterToGhosted("face");

  int dim = mesh_->space_dimension();
  AmanziMesh::Entity_ID_List faces;
  AmanziMesh::Entity_ID_List cells;

  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int nfaces = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);

  Epetra_MultiVector& Dff_f = *Dff_->ViewComponent("boundary_face",false);
  const Epetra_Map& fb_map = mesh_->exterior_face_map(false);
  const Epetra_Map& f_map = mesh_->face_map(false);

  if (Afc_cells_.size() != nfaces) {
    Afc_cells_.resize(static_cast<size_t>(nfaces));
  }
  if (Aff_cells_.size() != nfaces) {
    Aff_cells_.resize(static_cast<size_t>(nfaces));
  }

  for (int f=0; f!=nfaces; ++f) {
    if (Krel == Teuchos::null ||
        (!Krel->HasComponent("cell") && !Krel->HasComponent("face"))) {
      (*rel_perm_transmissibility_)[f] = (*transmissibility_)[f];
    }
    else if (Krel->HasComponent("cell") && !Krel->HasComponent("face")) {
      Errors::Message msg("Matrix_TPFA: finite volume discretization methods doesn't work with this rel_perm");
      Exceptions::amanzi_throw(msg);
    }
    else if (!Krel->HasComponent("cell") && Krel->HasComponent("face")) {
      const Epetra_MultiVector& Krel_f = *Krel->ViewComponent("face",true);
      (*rel_perm_transmissibility_)[f] = (*transmissibility_)[f] * Krel_f[0][f];
    }
    else if (Krel->HasComponent("cell") && Krel->HasComponent("face")) {
      //      const Epetra_MultiVector& Krel_f = *Krel->ViewComponent("face",true);
      //      (*rel_perm_transmissibility_)[f] = (*transmissibility_)[f] * Krel_f[0][f];
      Errors::Message msg("Matrix_TPFA: finite volume discretization methods doesn't work with this rel_perm");
      Exceptions::amanzi_throw(msg);
    }

    mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
    int mcells = cells.size();
    if (mcells == 1) {
      Teuchos::SerialDenseMatrix<int, double> Bff(1,1);
      Epetra_SerialDenseVector Bfc(1);
      int fb_lid = fb_map.LID(f_map.GID(f));

      Bff(0,0) =  (*rel_perm_transmissibility_)[f];
      Bfc(0)   =  -(*rel_perm_transmissibility_)[f];

      Aff_cells_[f] = Bff;
      Afc_cells_[f] = Bfc;
      Dff_f[0][fb_lid] = Bff(0,0);
    }
  }
}

void Matrix_TPFA::FillMatrixGraphs_(const Teuchos::Ptr<Epetra_CrsGraph> cf_graph,
        const Teuchos::Ptr<Epetra_FECrsGraph> ff_graph) {

  const Epetra_Map& cmap_wghost = mesh_->cell_map(true);
  const Epetra_Map& fbmap = mesh_->exterior_face_map(false);
  const Epetra_Map& fmap = mesh_->face_map(false);

  AmanziMesh::Entity_ID_List cells;
  int cell_GID, face_GID;
  int nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);

  for (int f=0; f < nfaces_owned; f++) {
    mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
    if (cells.size() == 1) {
      cell_GID = cmap_wghost.GID(cells[0]);
      face_GID = fmap.GID(f);
      cf_graph->InsertGlobalIndices(cell_GID, 1, &face_GID);
    }
  }

  // assemble the graphs
  int ierr = cf_graph->FillComplete();
  ASSERT(!ierr);
};


void Matrix_TPFA::UpdatePreconditioner_() const {
  /// additional block preconditioner
  if (Aff_pc_ == Teuchos::null) {
    if (plist_.isSublist("preconditioner")) {
      Teuchos::ParameterList pc_list = plist_.sublist("preconditioner");
      if (pc_list.get<std::string>("preconditioner type") == (std::string)("boomer amg")) {
        Errors::Message msg("Matrix_TPFA:: boomer amg solver doesn't work with boundary faces set.");
        Exceptions::amanzi_throw(msg);
      }
      AmanziPreconditioners::PreconditionerFactory pc_fac;
      Aff_pc_ = pc_fac.Create(pc_list);
    }
  }

  if ((S_pc_ == Teuchos::null)||(Aff_pc_ == Teuchos::null)) {
    Errors::Message msg("Matrix_TPFA::ApplyInverse called but no preconditioner sublist was provided");
    Exceptions::amanzi_throw(msg);
  }

  S_pc_->Destroy();
  S_pc_->Update(Spp_);

  Aff_pc_->Destroy();
  Aff_pc_->Update(Aff_);
}


/* ******************************************************************
 * Initialize Trilinos matrices. It must be called only once.
 * If matrix is non-symmetric, we generate transpose of the matrix
 ****************************************************************** */
void Matrix_TPFA::SymbolicAssembleGlobalMatrices() {
  // get the standard matrices
  const Epetra_Map& cmap = mesh_->cell_map(false);
  const Epetra_Map& cmap_wghost = mesh_->cell_map(true);
  const Epetra_Map& fbmap = mesh_->exterior_face_map(false);
  const Epetra_Map& fmap = mesh_->face_map(false);

  int avg_entries_row = (mesh_->space_dimension() == 2) ? MFD_QUAD_FACES : MFD_HEX_FACES;

  Epetra_FECrsGraph pp_graph(Copy, cmap, avg_entries_row + 1);
  Teuchos::RCP<Epetra_FECrsGraph> fbfb_graph =
      Teuchos::rcp(new Epetra_FECrsGraph(Copy, fbmap, avg_entries_row + 1));


  // allocate the graphs
  Teuchos::RCP<Epetra_CrsGraph> cf_graph =
      Teuchos::rcp(new Epetra_CrsGraph(Copy, cmap, fbmap, avg_entries_row, false));


  FillMatrixGraphs_(cf_graph.ptr(), fbfb_graph.ptr());
  // Acc_ = Teuchos::rcp(new Epetra_Vector(cmap));
  // Acf_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *cf_graph));
  // if (symmetric()) {
  //   Afc_ = Acf_;
  // } else {
  //   Afc_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *cf_graph));
  // }

  AmanziMesh::Entity_ID_List cells;
  int nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);

  // create the RHS
  std::vector<std::string> names(2);
  names[0] = "cell";
  names[1] = "boundary_face";

  std::vector<AmanziMesh::Entity_kind> locations(2);
  locations[0] = AmanziMesh::CELL;
  locations[1] = AmanziMesh::BOUNDARY_FACE;

  std::vector<int> num_dofs(2,1);
  CompositeVectorSpace space2;
  space2.SetMesh(mesh_)->SetGhosted()->SetComponents(names,locations,num_dofs);
  rhs_ = Teuchos::rcp(new CompositeVector(space2));


  // also create a cell-cell matrix for the TPFA

  int cells_GID[2], face_GID;

  int nfaces = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  for (int f = 0; f < nfaces; f++) {
    mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
    int ncells = cells.size();

    for (int n = 0; n < ncells; n++) cells_GID[n] = cmap_wghost.GID(cells[n]);
    pp_graph.InsertGlobalIndices(ncells, cells_GID, ncells, cells_GID);

    if (ncells == 1) {
      face_GID = fmap.GID(f);
      fbfb_graph->InsertGlobalIndices(1, &face_GID, 1, &face_GID);
    }
    // face_GID = fmap.GID(f);
    // ff_graph.InsertGlobalIndices(1, &face_GID, 1, &face_GID);

  }
  pp_graph.GlobalAssemble();  // Symbolic graph is complete.
  fbfb_graph->GlobalAssemble();
  //ff_graph.GlobalAssemble();

  // create global matrices
  std::vector<std::string> names1(1,"boundary_face");
  std::vector<AmanziMesh::Entity_kind> locations1(1,AmanziMesh::BOUNDARY_FACE);
  std::vector<int> ndofs1(1,1);

  CompositeVectorSpace space;
  space.SetMesh(mesh_)->SetGhosted()->SetComponents(names1,locations1,ndofs1);
  Dff_ = Teuchos::rcp(new CompositeVector(space));
  Dff_->PutScalar(0.);

  Spp_ = Teuchos::rcp(new Epetra_FECrsMatrix(Copy, pp_graph));
  Spp_->GlobalAssemble();

  Aff_ = Teuchos::rcp(new Epetra_FECrsMatrix(Copy, *fbfb_graph));
  Aff_->GlobalAssemble();


}


void Matrix_TPFA::AssembleSchur_() const {
  std::vector<std::string> names_c(1,"cell");
  std::vector<AmanziMesh::Entity_kind> locations_c(1,AmanziMesh::CELL);
  std::vector<std::string> names_f(1,"face");
  std::vector<AmanziMesh::Entity_kind> locations_f(1,AmanziMesh::FACE);
  std::vector<int> ndofs(1,1);

  AmanziMesh::Entity_ID_List faces;
  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  const Epetra_Map& fmap_wghost = mesh_->face_map(true);
  const Epetra_Map& fb_map = mesh_->exterior_face_map(false);
  const Epetra_Map& cmap = mesh_->cell_map(false);
  const Epetra_Map& cmap_wghost = mesh_->cell_map(true);
  AmanziMesh::Entity_ID_List cells;
  int cells_GID[2];
  int face_GID;
  // assemble rhs

  Epetra_MultiVector& rhs_cells = *rhs_->ViewComponent("cell",false);
  Epetra_MultiVector& rhs_bf = *rhs_->ViewComponent("boundary_face",false);
  Epetra_MultiVector& Dff_f = *Dff_->ViewComponent("boundary_face",false);


  // std::cout<<"rhs_cells\n"<<rhs_cells<<"\n";
  // std::cout<<"rhs_bf\n"<<rhs_bf<<"\n";
  // exit(0);

  Spp_->PutScalar(0.0);
  Aff_->PutScalar(0.0);

  for (AmanziMesh::Entity_ID f=0; f!=nfaces_owned; ++f) {
    mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
    int mcells = cells.size();

    // populate face-based matrix.
    Teuchos::SerialDenseMatrix<int, double> Bpp(mcells, mcells);
    if (mcells == 2) {
      for (int i=0; i<mcells; i++) {
        int c = cells[i];
        cells_GID[i] = cmap_wghost.GID(c);
        // double tij;
        double tij = (*rel_perm_transmissibility_)[f];
        Bpp(i,i) = tij;
        for (int j=i+1; j<mcells; j++) {
          Bpp(i,j) = -tij; //Bpp(i,j) =0.;
          Bpp(j,i) = -tij; //Bpp(j,i) = 0.;
        }
      }
    }
    else if (mcells == 1) {
      int c = cells[0];
      cells_GID[0] = cmap_wghost.GID(c);
      double val = Afc_cells_[f](0);
      face_GID = fmap_wghost.GID(f);
      //(*Acf_).ReplaceGlobalValues(cells_GID[0], 1, &val, &face_GID);
      int fb_lid = fb_map.LID(face_GID);

      /// Afc_cells_[f](0) = 0               if DIRICHLET boundary
      /// Afc_cells_[f](0) = rel_perm_trans  else
      /// Aff_cells_[f](0,0) = rel_perm_tran

      if (fabs(Afc_cells_[f](0)) < 1e-23) {  ///
        Bpp(0,0) = Aff_cells_[f](0,0);
        Dff_f[0][fb_lid] = 1.;
      }
      else {
        Bpp(0,0) = Aff_cells_[f](0,0);
        Dff_f[0][fb_lid] = Aff_cells_[f](0,0);
      }


      (*Aff_).SumIntoGlobalValues(face_GID, 1,  &(Dff_f[0][fb_lid]), &face_GID);

      /// Schur comlement contribution for rhs_cells (old verion)
      //  rhs_cells[0][c] += rhs_bf[0][fb_lid]*Afc_cells_[f](0) / Dff_f[0][fb_lid];

    }
    (*Spp_).SumIntoGlobalValues(mcells, cells_GID, Bpp.values());

    // double val=100;
    // face_GID = fmap_wghost.GID(f);
    // Att_->SumIntoGlobalValues(face_GID, 1,  &val, &face_GID);

  }

  (*Spp_).GlobalAssemble();
  (*Aff_).GlobalAssemble();

  //(*Att_).GlobalAssemble();
  //(*Acf_).FillComplete();


  //std::cout<< (*Spp_);
  // std::cout<<"rhs_cells\n"<<rhs_cells<<"\n";
  //std::cout<< (*Aff_);
  //std::cout<< Dff_f;
  //exit(0);
  // tag matrices as assembled
  assembled_operator_ = true;
  assembled_schur_ = true;
  Sff_ = Spp_;
}


void
Matrix_TPFA::ComputeSchurComplement(const std::vector<MatrixBC>& bc_markers,
        const std::vector<double>& bc_values) {
  //AssertAssembledOperator_or_die_();
  //assembled_schur_ = true;
  //Sff_ = Spp_;
}



void Matrix_TPFA::AssembleRHS_() const {

  const Epetra_Map& fmap_wghost = mesh_->face_map(true);
  const Epetra_Map& fb_map = mesh_->exterior_face_map(false);

  Epetra_MultiVector& rhs_cells = *rhs_->ViewComponent("cell",false);
  Epetra_MultiVector& rhs_bf = *rhs_->ViewComponent("boundary_face",false);

  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  AmanziMesh::Entity_ID_List faces;
  /// Add Gravity to RHS
  ///
  for (int c=0; c != ncells_owned; ++c) {
    rhs_cells[0][c] -= Fc_cells_[c];
    mesh_->cell_get_faces(c, &faces);
    for (int n=0; n<faces.size(); n++) {
      int f = faces[n];
      int f_gid = fmap_wghost.GID(f);
      int fb_lid = fb_map.LID(f_gid);
      if (fb_lid >= 0) {
        rhs_bf[0][fb_lid] -= Ff_cells_[c][n];
      }
    }
  }
  assembled_rhs_ = true;

}




/* ******************************************************************
 * Parallel matvec product Spp * Xc.
 ****************************************************************** */
int Matrix_TPFA::Apply(const CompositeVector& X,
                       CompositeVector& Y) const {

  if (!assembled_schur_) {
    AssembleSchur_();
    UpdatePreconditioner_();
  }

  int ierr = Spp_->Multiply(false, *X.ViewComponent("cell",false),
                            *Y.ViewComponent("cell",false));

  ASSERT(!ierr);


  const Epetra_Map& fb_map = mesh_->exterior_face_map(false);
  const Epetra_Map& f_map = mesh_->face_map(false);

  const  Epetra_MultiVector& Xc  = *X.ViewComponent("cell", false);
  const  Epetra_MultiVector& Xfb = *X.ViewComponent("boundary_face", false);
  Epetra_MultiVector& Yc  = *Y.ViewComponent("cell", false);
  Epetra_MultiVector& Yfb = *Y.ViewComponent("boundary_face", false);
  Epetra_MultiVector& Dff_f = *Dff_->ViewComponent("boundary_face",false);

  AmanziMesh::Entity_ID_List cells;
  int nb = fb_map.NumMyElements();

  ierr = Aff_->Multiply(false, Xfb, Yfb);
  ASSERT(!ierr);

  for (int fb=0; fb!=nb; ++fb) {
    int face_lid = f_map.LID(fb_map.GID(fb));
    Epetra_SerialDenseVector Bfc = Afc_cells_[face_lid];
    Teuchos::SerialDenseMatrix<int, double> Bff = Aff_cells_[face_lid];
    mesh_->face_get_cells(face_lid, AmanziMesh::Parallel_type::ALL, &cells);

    Yc[0][cells[0]] += Bfc(0) * Xfb[0][fb];
    Yfb[0][fb] += Bfc(0) * Xc[0][cells[0]];// + Dff_f[0][fb] * Xfb[0][fb];
    //Yfb[0][fb] = Dff_f[0][fb] * Xfb[0][fb];
  }


  if (ierr) {
    Errors::Message msg("Matrix_TPFA::Apply has failed to calculate y = A*x.");
    Exceptions::amanzi_throw(msg);
  }
  return ierr;
}

// Apply, for the cell-only system
// int Matrix_TPFA::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {
//   AssertAssembledOperator_or_die_();

//   int ierr = Spp_->Multiply(false, X, Y);
//   if (ierr) {
//     Errors::Message msg("Matrix_TPFA::Apply has failed to calculate y = A*x.");
//     Exceptions::amanzi_throw(msg);
//   }
//   return ierr;
// }

void Matrix_TPFA::ComputeNegativeResidual(const CompositeVector& solution,
        const Teuchos::Ptr<CompositeVector>& residual) const
{

  //Apply(solution, *residual);

  if (!assembled_rhs_) AssembleRHS_();

  const Epetra_MultiVector& uc = *solution.ViewComponent("cell", true);
  const Epetra_MultiVector& ub = *solution.ViewComponent("boundary_face");
  Epetra_MultiVector& rc  = *residual->ViewComponent("cell");
  Epetra_MultiVector& rbf = *residual->ViewComponent("boundary_face");

  AmanziMesh::Entity_ID_List cells;


  //std::cout << rc;

  //exit(0);

  // matvec product A*u
  solution.ScatterMasterToGhosted("cell");

  residual->PutScalar(0.0);
  int nfaces_wghost = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);
  int ncells_owned  = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

  // right-hand side
  Epetra_MultiVector& rhs_cell = *rhs_->ViewComponent("cell");
  Epetra_MultiVector& rhs_bf   = *rhs_->ViewComponent("boundary_face");
  Epetra_MultiVector& Dff_f = *Dff_->ViewComponent("boundary_face",false);

  const Epetra_Map& fb_map = mesh_->exterior_face_map(false);
  const Epetra_Map& f_map = mesh_->face_map(false);



  for (int f = 0; f < nfaces_wghost; f++) {
    mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
    int ncells = cells.size();

    if (ncells == 2) {
      int c1 = cells[0];
      int c2 = cells[1];

      double tmp = (*rel_perm_transmissibility_)[f] * (uc[0][c1] - uc[0][c2]);

      if (c1 < ncells_owned) rc[0][c1] += tmp;
      if (c2 < ncells_owned) rc[0][c2] -= tmp;

      // if ((c1==0)||(c2==0)){
      //    std::cout<<"c "<<c1<<" "<<c2<<" trans "<<(*rel_perm_transmissibility_)[f]<<" tmp "<<tmp<<" "<<uc[0][c1]<<" "<<uc[0][c2]<<"\n";
      // }

    } else {
      int c = cells[0];
      double bc_value = BoundaryValue(solution, f);
      Epetra_SerialDenseVector Bfc = Afc_cells_[f];
      Teuchos::SerialDenseMatrix<int, double> Bff = Aff_cells_[f];
      int face_gid = f_map.GID(f);
      int face_lbid = fb_map.LID(face_gid);

      double tmp = Bff(0,0) * uc[0][c] + Bfc(0) * bc_value;
      if (c < ncells_owned) rc[0][c] += tmp;

      // Residual on the boundary
      double rhs_val = BoundaryValue(*rhs_, f);
      double res_bf_val =  Dff_f[0][face_lbid] * bc_value + Bfc(0) * uc[0][c] - rhs_val;

      rbf[0][face_lbid] += res_bf_val;
      //std::cout<<Dff_f[0][face_lbid]<<" "<< bc_value<<" "<< Bfc(0) <<" "<< uc[0][c]<<"\n";
      //std::cout<<"rbf "<<" "<<face_lbid<<" "<<rbf[0][face_lbid]<<"\n";

    }
    //std::cout<<"f "<<f<<" "<<(*gravity_term_)[f]<<" "<<(*rel_perm_transmissibility_)[f]<<"\n";
  }

  AmanziMesh::Entity_ID_List faces;
  //std::cout<<"\n";
  for (int c = 0; c < ncells_owned; c++) {
    rc[0][c] -= rhs_cell[0][c];
  }

  //exit(0);

}

void Matrix_TPFA::DeriveFlux(const CompositeVector& solution,
                             const Teuchos::Ptr<CompositeVector>& mass_flux) const {


  solution.ScatterMasterToGhosted("cell");
  const Epetra_MultiVector& p = *solution.ViewComponent("cell", true);
  const Epetra_MultiVector& lb = *solution.ViewComponent("boundary_face", false);
  // const Epetra_MultiVector& Krel_face = *rel_perm->ViewComponent("face", true);
  Epetra_MultiVector& flux = *mass_flux->ViewComponent("face", true);

  const Epetra_Map& fb_map = mesh_->exterior_face_map(false);
  const Epetra_Map& f_map = mesh_->face_map(false);

  int ncells_owned  = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int nfaces_wghost = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);
  int nfaces_owned  = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);

  AmanziMesh::Entity_ID_List cells, faces;
  std::vector<int> dirs;

  std::vector<int> flag(nfaces_wghost, 0);

  for (int c = 0; c < ncells_owned; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];


      if (f < nfaces_owned && !flag[f]) {
        mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
        if (cells.size() == 1) {
          int face_gid = f_map.GID(f);
          int face_lbid = fb_map.LID(face_gid);
          if (face_lbid >= 0) {
            double value = lb[0][face_lbid];
            flux[0][f] = dirs[n] * (*rel_perm_transmissibility_)[f] * (p[0][c] - value);// + (*gravity_term_)[f];
            //flux[0][f] *= Krel_face[0][f];
            flag[f] = 1;
          }
        }
        else {
          int c1 = cells[0];
          int c2 = cells[1];
          if (c == c1) {
            flux[0][f] = dirs[n] * (*rel_perm_transmissibility_)[f] * (p[0][c1] - p[0][c2]);// + (*gravity_term_)[f];
          } else {
            flux[0][f] = dirs[n] * (*rel_perm_transmissibility_)[f] * (p[0][c2] - p[0][c1]);// + (*gravity_term_)[f];
          }
          //flux[0][f] *= Krel_face[0][f];
          flag[f] = 1;
        }
      }
    }
  }
}




int Matrix_TPFA::ApplyInverse(const CompositeVector& X,
        CompositeVector& Y) const {
  //AssertAssembledOperator_or_die_();
  //AssertAssembledSchur_or_die_();
  if (!assembled_schur_) {
    AssembleSchur_();
    UpdatePreconditioner_();
  }

  // Solve the Schur complement system Spp * Yc = Xc.
  int ierr = 0;
  const Epetra_MultiVector& Xc = *X.ViewComponent("cell",false);
  const Epetra_MultiVector& Xb = *X.ViewComponent("boundary_face",false);
  Epetra_MultiVector Tc(Xc);
  int nb_faces = Xb.MyLength();

  const Epetra_Map& fb_map = mesh_->exterior_face_map(false);
  const Epetra_Map& f_map = mesh_->face_map(false);
  AmanziMesh::Entity_ID_List cells;


  Teuchos::ParameterList  plist;
  Teuchos::ParameterList& pre_list = plist.sublist("gmres");
  Teuchos::ParameterList& slist = pre_list.sublist("gmres parameters");

  pre_list.set<std::string>("iterative method", "gmres");
  slist.set<double>("error tolerance", 1e-12);
  slist.set<int>("maximum number of iterations", 10000);
  Teuchos::ParameterList& vlist = slist.sublist("VerboseObject");
  //vlist.set("Verbosity Level", "extreme");
  vlist.set("Verbosity Level", "high");

  // delegating preconditioning to the base operator
  Teuchos::RCP<const Matrix_TPFA> op_matrix = Teuchos::rcp(this, false);
  Teuchos::RCP<BlockMatrix> op_prec   = Teuchos::rcp(new BlockMatrix(space_));

  op_prec -> SetNumberofBlocks(2);
  op_prec -> SetBlock(0, Spp_);
  op_prec -> SetBlock(1, Aff_);
  op_prec -> SetPrec(0, S_pc_);
  op_prec -> SetPrec(1, Aff_pc_);
  op_prec -> prec_list = plist_;

  AmanziSolvers::LinearOperatorFactory< CompositeMatrix, CompositeVector, CompositeVectorSpace> factory;
  Teuchos::RCP<AmanziSolvers::LinearOperator< CompositeMatrix, CompositeVector, CompositeVectorSpace> >
      solver = factory.Create("gmres", plist, op_matrix, op_prec);


  Y.PutScalar(0.0);
  ierr = solver->ApplyInverse(X, Y);
  // op_matrix -> Apply(X,Y);
  // op_prec -> ApplyInverse(X,Y);

  //exit(0);
  return ierr;
}


void Matrix_TPFA::UpdateConsistentFaceConstraints(
    const Teuchos::Ptr<CompositeVector>& u) {

  //    Teuchos::RCP<const Epetra_MultiVector> uc = u->ViewComponent("cell", false);
  //   Epetra_MultiVector& uf = *u->ViewComponent("face", false);

  //   Epetra_MultiVector& Dff_f = *Dff_->ViewComponent("face",false);
  //   Epetra_MultiVector& rhs_f = *rhs_->ViewComponent("face", false);
  //   Epetra_MultiVector update_f(rhs_f);

  //   Afc_->Multiply(true,*uc, update_f);  // Afc is kept in the transpose form.
  //   update_f.Update(1.0, rhs_f, -1.0);

  //   int nfaces = rhs_f.MyLength();
  //   for (int f=0; f!=nfaces; ++f) {
  //     uf[0][f] = update_f[0][f] / Dff_f[0][f];
  //   }
}

void Matrix_TPFA::UpdateConsistentFaceCorrection(const CompositeVector& u,
        const Teuchos::Ptr<CompositeVector>& Pu) {
  //   //AssertAssembledOperator_or_die_();

  //   Teuchos::RCP<const Epetra_MultiVector> Pu_c = Pu->ViewComponent("cell", false);
  //   Epetra_MultiVector& Pu_f = *Pu->ViewComponent("face", false);
  //   const Epetra_MultiVector& u_f = *u.ViewComponent("face", false);

  //   Epetra_MultiVector& Dff_f = *Dff_->ViewComponent("face",false);

  //   Afc_->Multiply(true, *Pu_c, Pu_f);  // Afc is kept in the transpose form.
  //   Pu_f.Update(1., u_f, -1.);

  //   int nfaces = Pu_f.MyLength();
  //   for (int f=0; f!=nfaces; ++f) {
  //     Pu_f[0][f] /= Dff_f[0][f];
  //   }
}


// MANY ASSUMPTIONS!  THIS IS NOT GENERAL!
// -- K_abs = 1
// -- upwinding is potential difference, with height overlap
void Matrix_TPFA::AnalyticJacobian(const Upwinding& upwinding,
        const Teuchos::Ptr<State>& S,
        std::string potential_key,
        const CompositeVector& dconductivity,
        const std::vector<MatrixBC>& bc_markers,
        const std::vector<double>& bc_values) {

  //AssertAssembledOperator_or_die_();

  // maps and counts
  AmanziMesh::Entity_ID_List faces;
  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  const Epetra_Map& cmap_wghost = mesh_->cell_map(true);

  // local work arrays
  AmanziMesh::Entity_ID_List cells;
  int cells_GID[2];
  int ierr = 0;

  // Get the derivatives
  std::vector<Teuchos::RCP<Teuchos::SerialDenseMatrix<int, double> > > Jpp_faces;
  upwinding.UpdateDerivatives(S, potential_key, dconductivity, bc_markers, bc_values, &Jpp_faces);
  ASSERT(Jpp_faces.size() == nfaces_owned);

  // Assemble into Spp
  for (unsigned int f=0; f!=nfaces_owned; ++f) {
    mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);

    int mcells = cells.size();
    for (int n=0; n!=mcells; ++n) {
      cells_GID[n] = cmap_wghost.GID(cells[n]);
    }
    ierr = (*Spp_).SumIntoGlobalValues(mcells, cells_GID, Jpp_faces[f]->values());
    ASSERT(!ierr);
  }

  // finish assembly
  ierr = Spp_->GlobalAssemble();
  ASSERT(!ierr);
}

/* ******************************************************************
 * Compute transmissibilities on faces
 ****************************************************************** */
void Matrix_TPFA::ComputeTransmissibilities_(const Teuchos::Ptr<std::vector<WhetStone::Tensor> >& K)
{
  transmissibility_->PutScalar(0.0);
  rel_perm_transmissibility_->PutScalar(0.0);

  //double rho_ = *S_->GetScalarData("fluid_density");
  //double mu_ = *S_->GetScalarData("fluid_viscosity");
  //const Epetra_Vector& gravity_ = *S_->GetConstantVectorData("gravity");

  int dim = mesh_->space_dimension();
  AmanziGeometry::Point gravity(dim);

  WhetStone::Tensor Kc;
  if (K == Teuchos::null) {
    Kc.Init(mesh_->space_dimension(), 1);
    Kc(0,0) = 1.0;
  }

  // Gravity tems are computed with respect to vector (0,0,-1)^T
  // and should be rescaled later by real gravity
  for (int k = 0; k < dim-1; k++) gravity[k] = 0;
  gravity[dim-1] = -1.;


  AmanziMesh::Entity_ID_List faces;
  AmanziMesh::Entity_ID_List cells;
  AmanziGeometry::Point a_dist, a[2];
  double h[2], perm[2], perm_test[2], h_test[2], beta[2];
  double trans_f;

  int nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  int nfaces_wghost = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);

  for (int f = 0; f < nfaces_owned; f++) {
    mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
    int ncells = cells.size();
    //const AmanziGeometry::Point& normal = mesh_->face_normal(f, false, cells[0]);
    const AmanziGeometry::Point& normal = mesh_->face_normal(f, false);
    const AmanziGeometry::Point& face_centr = mesh_->face_centroid(f);
    double area = mesh_->face_area(f);

    if (ncells == 2) {
      a_dist = mesh_->cell_centroid(cells[1]) - mesh_->cell_centroid(cells[0]);
    } else if (ncells == 1) {
      a_dist = face_centr - mesh_->cell_centroid(cells[0]);
    }

    a_dist *= 1./norm(a_dist);

    for (int i=0; i<ncells; i++) {
      a[i] = face_centr - mesh_->cell_centroid(cells[i]);
      h[i] = norm(a[i]);
      double s = area / h[i];
      if (K != Teuchos::null) {
        Kc = (*K)[cells[i]];
      }
      perm[i] =  ((Kc * a[i]) * normal) * s;
      //perm_test[i] = (rho_/mu_) * (((*K_)[cells[i]] * normal) * a_dist);
      //h_test[i] = pow(-1.0, i)*((face_centr - mesh_->cell_centroid(cells[i]))*normal) / area;
      double dxn = a[i]*normal;

      beta[i] = fabs(perm[i] / dxn);
    }


    double  grav;
    grav = (gravity * normal) / area;

    trans_f = 0.0;

    if (ncells == 2) {
      grav *= (h[0] + h[1]);
      trans_f = (beta[0]*beta[1]) / (beta[0] + beta[1]);
    } else if (ncells == 1) {
      grav *= h[0];
      trans_f = beta[0];
    }

    (*transmissibility_)[f] = trans_f;
    (*gravity_term_)[f] = (*transmissibility_)[f] * grav;

  }

  // parallelization using CV capability
  CompositeVectorSpace cvs;
  cvs.SetMesh(mesh_);
  cvs.SetGhosted(true);
  cvs.SetComponent("face", AmanziMesh::FACE, 1);

  CompositeVector tmp(cvs, true);
  Epetra_MultiVector& data = *tmp.ViewComponent("face", true);

  data = *transmissibility_;
  tmp.ScatterMasterToGhosted("face", true);
  for (int f = nfaces_owned; f < nfaces_wghost; f++) {
    (*transmissibility_)[f] = data[0][f];
  }

  data = *gravity_term_;
  tmp.ScatterMasterToGhosted("face", true);
  for (int f = nfaces_owned; f < nfaces_wghost; f++) {
    (*gravity_term_)[f] = data[0][f];
  }
}

void Matrix_TPFA::CreateMFDmassMatrices(const Teuchos::Ptr<std::vector<WhetStone::Tensor> >& K) {

  int ok;
  nokay_ = npassed_ = 0;

  assembled_schur_ = false;
  assembled_operator_ = false;

  const Epetra_BlockMap& fmap_wghost = mesh_->face_map(true);
  transmissibility_ = Teuchos::rcp(new Epetra_Vector(fmap_wghost));
  rel_perm_transmissibility_ = Teuchos::rcp(new Epetra_Vector(fmap_wghost));
  gravity_term_ = Teuchos::rcp(new Epetra_Vector(fmap_wghost));


  ComputeTransmissibilities_(K);

  //std::cout<<*transmissibility_<<endl;

  // std::cout<<"End of CreateMFDmassMatrices in TPFA\n";
  // exit(0);

}

// void Matrix_TPFA::CreateMFDrhsVectors(){

//   int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

//   if (Fc_cells_.size() != ncells) {
//     Fc_cells_.resize(static_cast<size_t>(ncells));
//   }
//   for (int c=0; c!=ncells; ++c) Fc_cells_[c] = 0.;
// }



void Matrix_TPFA::ApplyBoundaryConditions(const std::vector<MatrixBC>& bc_markers,
        const std::vector<double>& bc_values, bool ADD_BC_FLUX) {

  // tag global matrices as invalid
  assembled_schur_ = false;
  assembled_operator_ = false;

  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int nfaces = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  const Epetra_Map& fb_map = mesh_->exterior_face_map(false);
  const Epetra_Map& f_map = mesh_->face_map(false);
  AmanziMesh::Entity_ID_List faces;
  AmanziMesh::Entity_ID_List cells;

  Epetra_MultiVector& rhs_cells = *rhs_->ViewComponent("cell",false);
  Epetra_MultiVector& rhs_bc_faces = *rhs_->ViewComponent("boundary_face",false);
  Epetra_MultiVector& Dff_f = *Dff_->ViewComponent("boundary_face",false);

  rhs_cells.PutScalar(0.0);
  rhs_bc_faces.PutScalar(0.0);

  for (int f=0; f != nfaces; ++f) {

    if (bc_markers[f] == MATRIX_BC_DIRICHLET) {
      int face_gid = f_map.GID(f);
      int face_lbid = fb_map.LID(face_gid);
      mesh_->face_get_cells(f, AmanziMesh::Parallel_type::OWNED, &cells);

      //Epetra_SerialDenseVector& Bfc = Afc_cells_[f];

      rhs_cells[0][cells[0]] -= Afc_cells_[f](0)*bc_values[f];
      Afc_cells_[f](0) = 0.;
      Dff_f[0][face_lbid] = 1;

      //Aff_cells_[f](0,0) = 1.;
      //rhs_bc_faces[0][face_lbid] =  Aff_cells_[f](0,0)*bc_values[f];
      rhs_bc_faces[0][face_lbid] = bc_values[f];

    }
    else if  ((bc_markers[f] == MATRIX_BC_FLUX)&&(ADD_BC_FLUX)) {
      int face_gid = f_map.GID(f);
      int face_lbid = fb_map.LID(face_gid);
      rhs_bc_faces[0][face_lbid] = bc_values[f] * mesh_->face_area(f);
    }

  }




}


double Matrix_TPFA::BoundaryValue(const Amanzi::CompositeVector& solution, int face_id) const
{
  double value=0.;

  if (solution.HasComponent("face")) {
    const Epetra_MultiVector& pres = *solution.ViewComponent("face",false);
    value = pres[0][face_id];
  }
  else if  (solution.HasComponent("boundary_face")) {
    const Epetra_MultiVector& pres = *solution.ViewComponent("boundary_face",false);
    const Epetra_Map& fb_map = mesh_->exterior_face_map(false);
    const Epetra_Map& f_map = mesh_->face_map(false);

    int face_gid = f_map.GID(face_id);
    int face_lbid = fb_map.LID(face_gid);

    value =  pres[0][face_lbid];
  }
  else {
    Errors::Message msg("No component is defined for boundary faces\n");
    Exceptions::amanzi_throw(msg);
  }

  return value;
}

void Matrix_TPFA::SetBoundaryValue(Amanzi::CompositeVector& solution, int face_id, double value)
{

  if (solution.HasComponent("face")) {
    Epetra_MultiVector& pres = *solution.ViewComponent("face",false);
    pres[0][face_id] = value;
  }
  else if  (solution.HasComponent("boundary_face")) {
    Epetra_MultiVector& pres = *solution.ViewComponent("boundary_face",false);
    const Epetra_Map& fb_map = mesh_->exterior_face_map(false);
    const Epetra_Map& f_map = mesh_->face_map(false);

    int face_gid = f_map.GID(face_id);
    int face_lbid = fb_map.LID(face_gid);

    pres[0][face_lbid] = value;
  }
  else {
    Errors::Message msg("No component is defined for boundary faces\n");
    Exceptions::amanzi_throw(msg);
  }

}

}  // namespace AmanziFlow
}  // namespace Amanzi
