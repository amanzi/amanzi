/*
This is the flow component of the Amanzi code. 

Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided in the top-level COPYRIGHT file.

Authors: Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
         Daniil Svyatskiy (dasvyat@lanl.gov)
*/

#include <vector>

#include "Epetra_FECrsGraph.h"
#include "AztecOO.h"
#include "mfd3d_diffusion.hh"

#include "Flow_constants.hh"
#include "Matrix_MFD_TPFA.hh"


namespace Amanzi {
namespace AmanziFlow {

template<class T>
int FindPosition(const std::vector<T>& v, const T& value) {
  for (int i = 0; i < v.size(); i++)
     if (v[i] == value) return i;
  return -1;
}


/* ******************************************************************
* Constructor.                                           
****************************************************************** */
Matrix_MFD_TPFA::Matrix_MFD_TPFA(Teuchos::RCP<Flow_State> FS, const Epetra_Map& map) 
   :  Matrix_MFD(FS, map)
{
  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED); 
  Acc_cells_.resize(ncells);
  Fc_cells_.resize(ncells);
}


/* ******************************************************************
* Calculate elemental stiffness matrices.                                            
****************************************************************** */
void Matrix_MFD_TPFA::CreateMFDstiffnessMatrices(RelativePermeability& rel_perm)
{
  int dim = mesh_->space_dimension();
  WhetStone::MFD3D_Diffusion mfd(mesh_);
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  Aff_cells_.clear();
  Afc_cells_.clear();
  Acf_cells_.clear();
  Acc_cells_.clear();

  Epetra_Vector& Krel_cells = rel_perm.Krel_cells();
  Epetra_Vector& Krel_faces = rel_perm.Krel_faces();

  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c = 0; c < ncells; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    Teuchos::SerialDenseMatrix<int, double>& Mff = Mff_cells_[c];
    Teuchos::SerialDenseMatrix<int, double> Bff(nfaces, nfaces);
    Epetra_SerialDenseVector Bcf(nfaces), Bfc(nfaces);

    if (Krel_cells[c] == 1.0) {
      for (int n = 0; n < nfaces; n++) Bff(n, n) = Mff(n, n) * Krel_faces[faces[n]];
    } else {
      for (int n = 0; n < nfaces; n++) Bff(n, n) = Mff(n, n) * Krel_cells[c];
    }

    double matsum = 0.0;  // elimination of mass matrix
    for (int n = 0; n < nfaces; n++) {
      double rowsum = Bff(n, n), colsum = Bff(n, n);
      Bcf(n) = -colsum;
      Bfc(n) = -rowsum;
      matsum += colsum;
    }

    Aff_cells_.push_back(Bff);  // This the only place where memory can be allocated.
    Afc_cells_.push_back(Bfc);
    Acf_cells_.push_back(Bcf);
    Acc_cells_.push_back(matsum);
  }
}


/* ******************************************************************
* Initialize Trilinos matrices. It must be called only once. 
* If matrix is non-symmetric, we generate transpose of the matrix 
* block Afc to reuse cf_graph; otherwise, pointer Afc = Acf.   
****************************************************************** */
void Matrix_MFD_TPFA::SymbolicAssembleGlobalMatrices(const Epetra_Map& super_map)
{
  Matrix_MFD::SymbolicAssembleGlobalMatrices(super_map);

  const Epetra_Map& cmap = mesh_->cell_map(false);
  const Epetra_Map& cmap_wghost = mesh_->cell_map(true);
  const Epetra_Map& fmap_wghost = mesh_->face_map(true);

  int avg_entries_row = (mesh_->space_dimension() == 2) ? FLOW_QUAD_FACES : FLOW_HEX_FACES;
  Epetra_FECrsGraph pp_graph(Copy, cmap, avg_entries_row + 1);

  AmanziMesh::Entity_ID_List cells;
  int cells_GID[2];

  int nfaces = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  for (int f = 0; f < nfaces; f++) {
    mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
    int ncells = cells.size();

    for (int n = 0; n < ncells; n++)
        cells_GID[n] = cmap_wghost.GID(cells[n]);

    pp_graph.InsertGlobalIndices(ncells, cells_GID, ncells, cells_GID);
  }
  pp_graph.GlobalAssemble();  // Symbolic graph is complete.

  // create global matrices
  Dff_ = Teuchos::rcp(new Epetra_Vector(fmap_wghost));
  Spp_ = Teuchos::rcp(new Epetra_FECrsMatrix(Copy, pp_graph));
  Spp_->GlobalAssemble();
}


/* ******************************************************************
* Convert elemental mass matrices into stiffness matrices and 
* assemble them into four global matrices. 
* We need an auxiliary GHOST-based vector to assemble the RHS.
****************************************************************** */
// void Matrix_MFD_TPFA::AssembleGlobalMatrices()
// {
//   Matrix_MFD::AssembleGlobalMatrices();

//   AmanziMesh::Entity_ID_List faces;
//   std::vector<int> dirs;
//   int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);

//   Dff_->PutScalar(0.0);
//   for (int c = 0; c < ncells_owned; c++) {
//     mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
//     int mfaces = faces.size();

//     for (int n = 0; n < mfaces; n++) {
//       int f = faces[n];
//       (*Dff_)[f] += Aff_cells_[c](n, n);
//     }
//   }
//   FS_->CombineGhostFace2MasterFace(*Dff_, Add);

//   // convert right-hand side to a cell-based vector
//   const Epetra_Map& cmap = mesh_->cell_map(false);
//   Epetra_Vector Tc(cmap);
//   int nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);

//   for (int f = 0; f < nfaces_owned; f++) (*rhs_faces_)[f] /= (*Dff_)[f];
//   (*Acf_).Multiply(false, *rhs_faces_, Tc);
//   for (int c = 0; c < ncells_owned; c++) (*rhs_cells_)[c] -= Tc[c];

//   rhs_faces_->PutScalar(0.0);

//   // create a with-ghost copy of Acc
//   const Epetra_Map& cmap_wghost = mesh_->cell_map(true);
//   Epetra_Vector Dcc(cmap_wghost);

//   for (int c = 0; c < ncells_owned; c++) Dcc[c] = (*Acc_)[c];
//   FS_->CopyMasterCell2GhostCell(Dcc);

//   AmanziMesh::Entity_ID_List cells;
//   int cells_GID[2];
//   double Acf_copy[2];

//   // create auxiliaty with-ghost copy of Acf_cells
//   const Epetra_Map& fmap_wghost = mesh_->face_map(true);
//   Epetra_Vector Acf_parallel(fmap_wghost);

//   for (int c = 0; c < ncells_owned; c++) {
//     mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
//     int mfaces = faces.size();

//     for (int i = 0; i < mfaces; i++) {
//       int f = faces[i];
//       if (f >= nfaces_owned) Acf_parallel[f] = Acf_cells_[c][i];
//     }
//   }
//   FS_->CombineGhostFace2MasterFace(Acf_parallel, Add);

//   // populate the global matrix
//   Spp_->PutScalar(0.0);
//   for (AmanziMesh::Entity_ID f = 0; f < nfaces_owned; f++) {
//     mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
//     int mcells = cells.size();

//     // populate face-based matrix.
//     Teuchos::SerialDenseMatrix<int, double> Bpp(mcells, mcells);
//     for (int n = 0; n < mcells; n++) {
//       int c = cells[n];
//       cells_GID[n] = cmap_wghost.GID(c);

//       mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
//       int i = FindPosition<AmanziMesh::Entity_ID>(faces, f);
//       Bpp(n, n) = Dcc[c] / faces.size();
//       if (c < ncells_owned) {
//         int i = FindPosition<AmanziMesh::Entity_ID>(faces, f);
//         Acf_copy[n] = Acf_cells_[c][i];
//       } else {
//         Acf_copy[n] = Acf_parallel[f];
//       }
//     }

//     for (int n = 0; n < mcells; n++) {
//       for (int m = n; m < mcells; m++) {
//         Bpp(n, m) -= Acf_copy[n] * Acf_copy[m] / (*Dff_)[f];
//         Bpp(m, n) = Bpp(n, m);
//       }
//     }

//     (*Spp_).SumIntoGlobalValues(mcells, cells_GID, Bpp.values());
//   }
//   (*Spp_).GlobalAssemble();

//   // std::cout<<(*Spp_)<<endl;
//   // std::cout<<(*rhs_cells_)<<endl;
//   // exit(0);

//}

  void Matrix_MFD_TPFA::AssembleGlobalMatrices( const Epetra_Vector& Krel_faces, const Epetra_Vector& Trans_faces)
{
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;
  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  int nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);

  const Epetra_Map& cmap_wghost = mesh_->cell_map(true);
  AmanziMesh::Entity_ID_List cells;

  Spp_->PutScalar(0.0);

  int cells_LID[2], cells_GID[2];
  Teuchos::SerialDenseMatrix<int, double> Spp_local(2, 2);

  for (int f = 0; f < nfaces_owned; f++){
    mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
    int mcells = cells.size();

    for (int n = 0; n < mcells; n++) {
      cells_LID[n] = cells[n];
      cells_GID[n] = cmap_wghost.GID(cells_LID[n]);     
    }
    double tij;
    for (int i=0; i < mcells; i++){
      for (int j=0; j < mcells; j++){
	tij = Trans_faces[f]*Krel_faces[f];
	if (i==j) Spp_local(i,j) = tij;
	else Spp_local(i,j) = -tij;
      }
    }
    //cout<<"face "<<f<<": "<<tij<<" ukvr "<<Krel_faces[f]<<endl;
    (*Spp_).SumIntoGlobalValues(mcells, cells_GID, Spp_local.values());

  }

  for (int c = 0; c <= ncells_owned; c++){
    //cout<<"Acc_cells "<<Acc_cells_[c]<<" Fc_cells_ "<<Fc_cells_[c]<<endl;
    cells_GID[0] = cmap_wghost.GID(c);
    int mcells = 1;
    Spp_local(0,0) = Acc_cells_[c];
    (*Spp_).SumIntoGlobalValues(mcells, cells_GID, Spp_local.values());

    (*rhs_cells_)[c] += Fc_cells_[c];
  }

  Spp_->GlobalAssemble();
  // int tmp;
  // cout<<"From AssembleGlobalMatrices\n";
  // std::cout<<(*Spp_)<<endl;
  // std::cout<<(*rhs_cells_)<<endl;
  //cin >> tmp;
  //exit(0);
    
}


/* ******************************************************************
* Assembles preconditioner. It has same set of parameters as matrix.
****************************************************************** */
  void Matrix_MFD_TPFA::AssembleSchurComplement(const Epetra_Vector& Krel_faces, const Epetra_Vector& Trans_faces)
{
  AssembleGlobalMatrices(Krel_faces, Trans_faces);
}

 
/* ******************************************************************
* Computation of the part of the jacobian which depends on
* analytical derivatives of relative permeabilities
****************************************************************** */
void Matrix_MFD_TPFA::AnalyticJacobian(
   const Epetra_Vector& solution, int dim,
   std::vector<int>& bc_models, std::vector<bc_tuple>& bc_values,
   const Epetra_Vector& trans_faces,
   const Epetra_Vector& grav_term_faces,
   RelativePermeability& rel_perm)
{
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;
  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  int nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);

  const Epetra_Map& cmap_wghost = mesh_->cell_map(true);
  AmanziMesh::Entity_ID_List cells;

  int cells_LID[2], cells_GID[2];
  double perm_abs_vert[2];
  double perm_abs_horz[2];
  double k_rel[2];
  double dk_dp[2];
  double pres[2];
  AmanziGeometry::Point cntr_cell[2];
  double dist;

  Epetra_Vector pres_gh(cmap_wghost);

  // Epetra_Vector pressure_vec(solution);
  FS_->CopyMasterCell2GhostCell(solution, pres_gh);

  // Epetra_Vector perm_vert_gh(cmap_wghost);
  // Epetra_Vector& perm_vert_vec = FS_->ref_vertical_permeability();
  // FS_->CopyMasterCell2GhostCell(perm_vert_vec, perm_vert_gh);

  // Epetra_Vector perm_horz_gh(cmap_wghost);
  // Epetra_Vector& perm_horz_vec = FS_->ref_horizontal_permeability();
  // FS_->CopyMasterCell2GhostCell(perm_horz_vec, perm_horz_gh);

  Epetra_Vector& Krel_cells = rel_perm.Krel_cells();
  Epetra_Vector& Krel_faces = rel_perm.Krel_faces();
  const Epetra_Vector& dKdP_cells = rel_perm.dKdP_cells();
  const Epetra_Vector& dKdP_faces = rel_perm.dKdP_faces();
  int method = rel_perm.method();

  for (int f = 0; f < nfaces_owned; f++){
    mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
    int mcells = cells.size();
    Teuchos::SerialDenseMatrix<int, double> Jpp(mcells, mcells);
    //AmanziGeometry::Point face_cntr = mesh_->face_centroid(f);

    for (int n = 0; n < mcells; n++) {
      cells_LID[n] = cells[n];
      cells_GID[n] = cmap_wghost.GID(cells_LID[n]);
      //perm_abs_vert[n] = perm_vert_gh[cells_LID[n]];
      //perm_abs_horz[n] = perm_horz_gh[cells_LID[n]];
      pres[n] = pres_gh[cells_LID[n]];
      //k_rel[n] = Krel_cells[cells_LID[n]];
      dk_dp[n] = dKdP_cells[cells_LID[n]];
      //cntr_cell[n] = mesh_->cell_centroid(cells_LID[n]);
    }

    // if (mcells == 2) {
    //   //dist = norm(cntr_cell[0] - cntr_cell[1]);
    // } else 
    if (mcells == 1) {
      //dist = norm(cntr_cell[0] - face_cntr);
      //k_rel[0] = Krel_faces[f];
      dk_dp[0] = dKdP_faces[f];
    }

    const AmanziGeometry::Point& normal = mesh_->face_normal(f, false, cells[0]);

    // if (fabs(normal[2]) > 0.3) {
    //   cout<<f<<" Krel "<<Krel_faces[f]<<" "<<dKdP_faces[f]<<" ";
    //   if (mcells == 1) cout<<cells_GID[0]<<" ";
    //   else cout<<cells_GID[0]<<" "<<cells_GID[1];
    //   cout<<endl;
    // }

    // AmanziGeometry::Point normal = mesh_->face_normal(f, false, cells_LID[0]);
    // normal *= 1./ mesh_->face_area(f);

    // ComputeJacobianLocal(mcells, f, dim, method, bc_models, bc_values, dist, pres,
    // 			 perm_abs_vert, perm_abs_horz, k_rel, dk_dp, normal, trans_faces, grav_term_faces, Jpp);

    ComputeJacobianLocal(mcells, f, method, bc_models, bc_values, pres, dk_dp, trans_faces, grav_term_faces, Jpp);

    // if (((cells_GID[0]>5)||(cells_GID[1]>5))&& (mcells > 1)){
    // if (f==107){
    //   cout<<"GID "<<cells_GID[0]<<" "<<cells_GID[1]<<endl;
    //   cout<<"dist "<<dist<<endl;
    //   cout<<"pres "<<pres[0]<<" "<<pres[1]<<endl;
    //   cout<<"perm v "<< perm_abs_vert[0]<<" "<< perm_abs_vert[1]<<endl;
    //   cout<<"perm h "<< perm_abs_horz[0]<<" "<< perm_abs_horz[1]<<endl;
    //   cout<<"k_rel "<<k_rel[0]<<" "<<k_rel[1]<<endl;
    //   cout<<"dk_dp "<<dk_dp[0]<<" "<<dk_dp[1]<<endl;
    //   cout<<"normal "<<normal<<endl;
    //   cout<<Jpp<<endl;
    // }

    (*Spp_).SumIntoGlobalValues(mcells, cells_GID, Jpp.values());
  }

  Spp_->GlobalAssemble();

  //cout<<(*Spp_)<<endl;
}


/* ******************************************************************
* Computation of a local submatrix of 
* Analytical Jacobian (nonlinear part) on a particular face.
****************************************************************** */
void Matrix_MFD_TPFA::ComputeJacobianLocal(int mcells,
                                           int face_id,
                                           //int dim,
                                           int Krel_method,
                                           std::vector<int>& bc_models,
                                           std::vector<bc_tuple>& bc_values,
                                           //double dist,
                                           double *pres,
                                           // double *perm_abs_vert,
                                           // double *perm_abs_horz,
                                           // double *k_rel,
                                           double *dk_dp_cell,
                                           //AmanziGeometry::Point& normal,
					   const Epetra_Vector& trans_faces,
					   const Epetra_Vector& grav_term_faces,
                                           Teuchos::SerialDenseMatrix<int, double>& Jpp)
{
  // double K[2];
  double dKrel_dp[2];

  // double rho_w = FS_->ref_fluid_density();

  // for (int c = 0; c < mcells; c++) {
  //   K[c] = 0.0;
  //   for (int i = 0; i < dim-1; i++) K[c] += normal[i]*normal[i]*perm_abs_horz[c];
  //   K[c] += normal[dim-1]*normal[dim-1]*perm_abs_vert[c];
  //   K[c] /= FS_->ref_fluid_viscosity();
  // }

  // double Kabs_dist;

  // if (mcells == 1) {
  //   Kabs_dist = K[0]/dist;
  // } else {
  //   Kabs_dist = 2*K[0]*K[1]/(dist*(K[0] + K[1]));
  // }

  double rho_w = FS_->ref_fluid_density();
  // AmanziGeometry::Point gravity(dim);
  // for (int i = 0; i < dim; i++) gravity[i] = (*(FS_->gravity()))[i] * rho;

  // double grn = dist*gravity*normal;

  //cout<<Kabs_dist*rho_w*mesh_->face_area(face_id)<<" "<<trans_faces[face_id]<<endl;

  //cout<<"grn "<<grn*Kabs_dist<<" "<<grav_term_faces[face_id]/(rho_w*mesh_->face_area(face_id))<<endl;
  //if (abs(grn) > 1e-8) exit(0);
  // for (int i = 0; i < dim; i++) cout<<"dist "<<(*(FS->gravity()))[i]<<endl;

  double dpres;
  if (mcells == 2) {
    dpres = pres[0] - pres[1];// + grn;
		// cout<<"pres[0] "<<pres[0]<<" pres[1] "<<pres[1]<<" grv "<<grn<<endl;
    if (Krel_method == FLOW_RELATIVE_PERM_UPWIND_GRAVITY) {  // Define K and Krel_faces
      double cos_angle = grav_term_faces[face_id]/(rho_w*mesh_->face_area(face_id));
      if (cos_angle > FLOW_RELATIVE_PERM_TOLERANCE) {  // Upwind
        dKrel_dp[0] = dk_dp_cell[0];
        dKrel_dp[1] = 0.0;
      } else if (cos_angle < -FLOW_RELATIVE_PERM_TOLERANCE) {  // Upwind
        dKrel_dp[0] = 0.0;
        dKrel_dp[1] = dk_dp_cell[1];
      } else if (fabs(cos_angle) < FLOW_RELATIVE_PERM_TOLERANCE) {  // Upwind
        dKrel_dp[0] = 0.5*dk_dp_cell[0];
        dKrel_dp[1] = 0.5*dk_dp_cell[1];
      }
    } else if (Krel_method == FLOW_RELATIVE_PERM_UPWIND_DARCY_FLUX) {
      if (trans_faces[face_id]*dpres + grav_term_faces[face_id] > FLOW_RELATIVE_PERM_TOLERANCE) {  // Upwind
        dKrel_dp[0] = dk_dp_cell[0];
        dKrel_dp[1] = 0.0;
    } else if (trans_faces[face_id]*dpres + grav_term_faces[face_id] < -FLOW_RELATIVE_PERM_TOLERANCE) {  // Upwind
        dKrel_dp[0] = 0.0;
        dKrel_dp[1] = dk_dp_cell[1];
    } else if (fabs(trans_faces[face_id]*dpres + grav_term_faces[face_id]) < FLOW_RELATIVE_PERM_TOLERANCE) {  // Upwind
        dKrel_dp[0] = 0.5*dk_dp_cell[0];
        dKrel_dp[1] = 0.5*dk_dp_cell[1];
    }
  } else if (Krel_method == FLOW_RELATIVE_PERM_ARITHMETIC_MEAN) {
    dKrel_dp[0] = 0.5*dk_dp_cell[0];
    dKrel_dp[1] = 0.5*dk_dp_cell[1];
  }

    // if (face_id == 107){
    //   cout<<"trans "<<trans_faces[face_id]<<endl;
    //   cout<<"dpres "<<dpres<<endl;
    //   cout<<"grav  "<<grav_term_faces[face_id]<<endl;
    //   cout<<"dKrel_dp "<<dKrel_dp[0]<<" "<<dKrel_dp[1]<<endl;
    // }

    Jpp(0, 0) = (trans_faces[face_id]*dpres + grav_term_faces[face_id])*dKrel_dp[0];
    Jpp(0, 1) = (trans_faces[face_id]*dpres + grav_term_faces[face_id])*dKrel_dp[1];
    Jpp(1, 0) = -Jpp(0, 0);
    Jpp(1, 1) = -Jpp(0, 1);

  } else if (mcells == 1) {
    if (bc_models[face_id] == FLOW_BC_FACE_PRESSURE) {                   
      pres[1] = bc_values[face_id][0];

      dpres = pres[0] - pres[1];// + grn;
      //cout<<"dk_dp_cell[0] "<<dk_dp_cell[0]<<endl;
      Jpp(0,0) = (trans_faces[face_id]*dpres + grav_term_faces[face_id]) * dk_dp_cell[0];
    } else {
      Jpp(0,0) = 0.0;
    }
  }
}


/* ******************************************************************
* TBW
****************************************************************** */
void Matrix_MFD_TPFA::AddCol2NumJacob(int cell, Epetra_Vector& r)
{
  WhetStone::MFD3D_Diffusion mfd3d(mesh_);

  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  mesh_->cell_get_faces_and_dirs(cell, &faces, &dirs);
  int mfaces = faces.size();
  int indices;
  double values;

  indices = cell;
  values = r[cell];

  NumJac_->ReplaceGlobalValues(cell, 1, &values, &indices);

  for (int i = 0; i < mfaces; i++) {
    int c_adj = mfd3d.cell_get_face_adj_cell(cell, faces[i]);
    // cout<<cell<<" - c_adj "<<c_adj<<endl;
    if (c_adj >= 0) {
      values = r[c_adj];
      NumJac_->ReplaceGlobalValues(c_adj, 1, &values, &indices);
    }
  }
  // cout<<"NonZeros "<<NumJac_->NumGlobalEntries(cell)<<endl;;

  NumJac_->GlobalAssemble();
}


/* ******************************************************************
* TBW
****************************************************************** */
void Matrix_MFD_TPFA::CompareJacobians()
{
  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);

  double ErrMax = 0;
  double ErrL2 = 0;
  double NormL2 = 0;
  double NormMax = 0;

  for (int cell = 0; cell < ncells_owned; cell++) {
    int nonzeros = NumJac_->NumGlobalEntries(cell);

    for (int i = 0; i < nonzeros; i++) {
      ErrL2 += ((*NumJac_)[cell][i] - (*Spp_)[cell][i])*((*NumJac_)[cell][i] - (*Spp_)[cell][i]);
      ErrMax = max(fabs((*NumJac_)[cell][i] - (*Spp_)[cell][i]), ErrMax);
      NormL2 += (*NumJac_)[cell][i] * (*NumJac_)[cell][i];
      NormMax = max((*NumJac_)[cell][i], NormMax);
    }
  }

  ErrL2 = sqrt(ErrL2/NormL2);

  cout<<"Error L2  "<<ErrL2<<endl;
  cout<<"Error Max "<<ErrMax/NormMax<<endl;
}


/* ******************************************************************
* Initialization of the preconditioner                                                 
****************************************************************** */
void Matrix_MFD_TPFA::InitPreconditioner(int method, Teuchos::ParameterList& prec_list)
{
  method_ = method;

  if (method_ == FLOW_PRECONDITIONER_TRILINOS_ML) {
    ML_list = prec_list;
    MLprec = Teuchos::rcp(new ML_Epetra::MultiLevelPreconditioner(*Spp_, ML_list, false));
  } else if (method_ == FLOW_PRECONDITIONER_HYPRE_AMG) {
#ifdef HAVE_HYPRE
    // read some boomer amg parameters
    hypre_ncycles = prec_list.get<int>("cycle applications", 5);
    hypre_nsmooth = prec_list.get<int>("smoother sweeps", 3);
    hypre_tol = prec_list.get<double>("tolerance", 0.0);
    hypre_strong_threshold = prec_list.get<double>("strong threshold", 0.0);
#endif
  } else if (method_ == FLOW_PRECONDITIONER_TRILINOS_BLOCK_ILU) {
    ifp_plist_ = prec_list;
  }
}


/* ******************************************************************
* Rebuild the preconditioner.                                                 
****************************************************************** */
void Matrix_MFD_TPFA::UpdatePreconditioner()
{
  //cout<<"In update preconditioner\n"<<*Spp_;

  if (method_ == FLOW_PRECONDITIONER_TRILINOS_ML) {
    if (MLprec->IsPreconditionerComputed()) MLprec->DestroyPreconditioner();
    MLprec->SetParameterList(ML_list);
    MLprec->ComputePreconditioner();
  } else if (method_ == FLOW_PRECONDITIONER_HYPRE_AMG) {
#ifdef HAVE_HYPRE
    IfpHypre_Spp_ = Teuchos::rcp(new Ifpack_Hypre(&*Spp_));
    Teuchos::RCP<FunctionParameter> functs[8];
    functs[0] = Teuchos::rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetCoarsenType, 0));
    functs[1] = Teuchos::rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetPrintLevel, 0));
    functs[2] = Teuchos::rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetNumSweeps, hypre_nsmooth));
    functs[3] = Teuchos::rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetMaxIter, hypre_ncycles));
    functs[4] = Teuchos::rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetRelaxType, 6));
    functs[5] = Teuchos::rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetStrongThreshold, hypre_strong_threshold));
    functs[6] = Teuchos::rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetTol, hypre_tol));
    functs[7] = Teuchos::rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetCycleType, 1));

    Teuchos::ParameterList hypre_list;
    hypre_list.set("Solver", PCG);
    hypre_list.set("Preconditioner", BoomerAMG);
    hypre_list.set("SolveOrPrecondition", Preconditioner);
    hypre_list.set("SetPreconditioner", true);
    hypre_list.set("NumFunctions", 8);
    hypre_list.set<Teuchos::RCP<FunctionParameter>*>("Functions", functs);

    IfpHypre_Spp_->SetParameters(hypre_list);
    IfpHypre_Spp_->Initialize();
    IfpHypre_Spp_->Compute();
#endif
  } else if (method_ == FLOW_PRECONDITIONER_TRILINOS_BLOCK_ILU) {
    Ifpack factory;
    std::string prectype("ILU");
    int ovl = ifp_plist_.get<int>("overlap", 0);
    ifp_plist_.set<std::string>("schwarz: combine mode", "Add");
    ifp_prec_ = Teuchos::rcp(factory.Create(prectype, &*Spp_, ovl));
    ifp_prec_->SetParameters(ifp_plist_);
    ifp_prec_->Initialize();
    ifp_prec_->Compute();
  }
}


/* ******************************************************************
* Parallel matvec product Spp * Xc.                                              
****************************************************************** */
int Matrix_MFD_TPFA::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  int nvectors = X.NumVectors();

  const Epetra_Map& cmap = mesh_->cell_map(false);
  const Epetra_Map& fmap = mesh_->face_map(false);

  // Create views Xc into the cell segments of X.
  double **cvec_ptrs = X.Pointers();
  double **fvec_ptrs = new double*[nvectors];
  for (int i = 0; i < nvectors; i++) fvec_ptrs[i] = cvec_ptrs[i] + ncells;

  Epetra_MultiVector Xc(View, cmap, cvec_ptrs, nvectors);
  Epetra_MultiVector Xf(View, fmap, fvec_ptrs, nvectors);

  // Create views Yc and Yf into the cell and face segments of Y.
  cvec_ptrs = Y.Pointers();
  for (int i = 0; i < nvectors; i++) fvec_ptrs[i] = cvec_ptrs[i] + ncells;

  Epetra_MultiVector Yc(View, cmap, cvec_ptrs, nvectors);
  Epetra_MultiVector Yf(View, fmap, fvec_ptrs, nvectors);

  int ierr = (*Spp_).Multiply(false, Xc, Yc);

  // cout<<"Xc"<<Xc<<endl;
  
  //cout<<"Apply TPFA"<<(*Spp_)<<endl;

  // cout<<"Yc "<<Yc<<endl;

  if (ierr) {
    Errors::Message msg("Matrix_MFD_TPFA::Apply has failed to calculate y = A*x.");
    Exceptions::amanzi_throw(msg);
  }

  // Yf = Xf;
  Yf.PutScalar(0.0);

  delete [] fvec_ptrs;
  return 0;
}


/* ******************************************************************
* The OWNED cell-based and face-based d.o.f. are packed together into 
* the X and Y Epetra vectors, with the cell-based in the first part.                                           
****************************************************************** */
int Matrix_MFD_TPFA::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  int nvectors = X.NumVectors();

  const Epetra_Map& cmap = mesh_->cell_map(false);
  const Epetra_Map& fmap = mesh_->face_map(false);

  // Create views Xc into the cell segments of X.
  double **cvec_ptrs = X.Pointers();
  double **fvec_ptrs = new double*[nvectors];
  for (int i = 0; i < nvectors; i++) fvec_ptrs[i] = cvec_ptrs[i] + ncells;

  Epetra_MultiVector Xc(View, cmap, cvec_ptrs, nvectors);
  Epetra_MultiVector Xf(View, fmap, fvec_ptrs, nvectors);

  // Create views Yc and Yf into the cell and face segments of Y.
  cvec_ptrs = Y.Pointers();
  for (int i = 0; i < nvectors; i++) fvec_ptrs[i] = cvec_ptrs[i] + ncells;

  Epetra_MultiVector Yc(View, cmap, cvec_ptrs, nvectors);
  Epetra_MultiVector Yf(View, fmap, fvec_ptrs, nvectors);

  // cout<<"Xc\n";
  // cout<<Xc<<endl;

  // Solve the Schur complement system Spp * Yc = Xc. Since AztecOO may
  // use the same memory for X and Y, we introduce auxiliaty vector Tc.
//   int ierr = 0;
//   Epetra_Vector Tc(cmap);
//   if (method_ == FLOW_PRECONDITIONER_TRILINOS_ML) {
//     MLprec->ApplyInverse(Xc, Tc);
//   } else if (method_ == FLOW_PRECONDITIONER_HYPRE_AMG) {
// #ifdef HAVE_HYPRE
//     ierr = IfpHypre_Spp_->ApplyInverse(Xc, Tc);
// #endif
//   } else if (method_ == FLOW_PRECONDITIONER_TRILINOS_BLOCK_ILU) {
//     ifp_prec_->ApplyInverse(Xc, Tc);
//   }
//   Yc = Tc;


  Epetra_LinearProblem problem(&*Spp_, &Yc, &Xc);

  AztecOO solver(problem);

  solver.SetAztecOption(AZ_solver, AZ_gmres);
  solver.SetAztecOption(AZ_output, AZ_none);
  solver.SetAztecOption(AZ_conv, AZ_rhs);
 
  int max_itrs_linear = 100;
  double convergence_tol_linear = 1e-10;

  solver.Iterate(max_itrs_linear, convergence_tol_linear);
  int num_itrs = solver.NumIters();

  // if (ierr) {
  //   Errors::Message msg("Matrix_MFD_TPFA::ApplyInverse has failed in calculating y = inv(A)*x.");
  //   Exceptions::amanzi_throw(msg);
  // }

  Yf = Xf;

  delete [] fvec_ptrs;
  return 0;
}


void  Matrix_MFD_TPFA::ApplyBoundaryConditions(std::vector<int>& bc_model, 
					       std::vector<bc_tuple>& bc_values,
					       const Epetra_Vector& Krel_faces,
					       Epetra_Vector& Trans_faces,
					       Epetra_Vector& grav_term_faces){

  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  //cout<<"rhs_cell\n"<<*rhs_cells_<<endl;

  rhs_cells_ -> PutScalar(0.);

  for (int c = 0; c < ncells; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];
      double value = bc_values[f][0];


      if (bc_model[f] == FLOW_BC_FACE_PRESSURE) {
	//if (f==84) cout<<"value*Trans_faces[f]*Krel_faces[f] "<<value*Trans_faces[f]*Krel_faces[f]<<" Trans_faces[f]*Krel_faces[f] "<< Trans_faces[f]*Krel_faces[f]<<endl;
	(*rhs_cells_)[c] += value*Trans_faces[f]*Krel_faces[f];
      } else if (bc_model[f] == FLOW_BC_FACE_FLUX) {
        (*rhs_cells_)[c] -= value * mesh_->face_area(f);
	//if (value < 0)	cout<<"FLOW_BC_FACE_FLUX "<<value * mesh_->face_area(f)<<" "<<value<<" "<<mesh_->face_area(f)<<endl;
	Trans_faces[f] = 0.0;
	grav_term_faces[f] =0.0;
      } else if (bc_model[f] == FLOW_BC_FACE_MIXED) {
	Errors::Message msg;
	msg << "Mixed boundary conditions are not supported in TPFA mode\n";
	Exceptions::amanzi_throw(msg);
      }

    }

  }

  //cout<<"Trans_faces\n"<<*Trans_faces<<endl;
  //cout<<"rhs_cell after ApplyBoundaryConditions\n"<<*rhs_cells_<<endl;
  //exit(0);
}


// /* ******************************************************************
// * Linear algebra operations with matrices: r = A * x - f                                                 
// ****************************************************************** */
double Matrix_MFD_TPFA::ComputeNegativeResidual(const Epetra_Vector& solution, 
						const Epetra_Vector& Krel_faces, 
						const Epetra_Vector& trans_faces, 
						Epetra_Vector& residual)
{
 
  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  int nfaces_wghost = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::USED);

  AmanziMesh::Entity_ID_List faces;
  AmanziMesh::Entity_ID_List cells;
  std::vector<int> dirs;

  const Epetra_Map& cmap_wghost = mesh_->cell_map(true);

  residual.PutScalar(0.0);
  Epetra_Vector sol_gh(cmap_wghost);

  FS_->CopyMasterCell2GhostCell(solution, sol_gh);


  //cout<<(*rhs_cells_)<<endl;

  for (int f = 0; f < nfaces_wghost; f++) {
    mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
    int ncells = cells.size();
    if (ncells > 1){
      for (int i = 0; i < ncells; i++){
	int c = cells[i];
	if (c >= ncells_owned) continue;
	residual[c] += pow(-1.0, i)*Krel_faces[f]*trans_faces[f]*(sol_gh[cells[0]] - sol_gh[cells[1]]);  
	// if (c==11) {
	//   cout<<"Negative Res2 "<<c<<" "<<Krel_faces[f]*trans_faces[f]*(sol_gh[cells[0]] - sol_gh[cells[1]])<<" "<<Krel_faces[f]*trans_faces[f]<<" "<<(sol_gh[cells[0]] - sol_gh[cells[1]])<<endl;
	// } 
      }
    }
    else if (ncells == 1){
      int c = cells[0];
      residual[c] += Krel_faces[f]*trans_faces[f]*(sol_gh[c]);
      //if (c==11) cout<<"Negative Res1 "<<c<<" "<<Krel_faces[f]*trans_faces[f]*(sol_gh[c])<<endl;
    }							
  } 
  


  for (int c = 0; c < ncells_owned; c++) {    
    residual[c] -= (*rhs_cells_)[c];
    //cout<<"residual "<<c<<": "<< residual[c]<<" "<<sol_gh[c]<<" "<<(*rhs_cells_)[c]<<endl;
  }
    
  //cout<<"residual "<<c<<": "<< residual[c]<<endl;
  //if (c >1) break;
  

  //exit(0);
  double norm_residual;
  residual.Norm2(&norm_residual);
  return norm_residual;
}

  void Matrix_MFD_TPFA::DeriveDarcyMassFlux(const Epetra_Vector& solution_cells,
					  const Epetra_Vector& Krel_faces,
					  const Epetra_Vector& Trans_faces,
					  const Epetra_Vector& Grav_term,
					  std::vector<int>& bc_model, 
					  std::vector<bc_tuple>& bc_values,
					  Epetra_Vector& darcy_mass_flux)
{
//   Teuchos::RCP<Epetra_Vector> solution_cell = Teuchos::rcp(FS_->CreateCellView(solution));
// #ifdef HAVE_MPI
//   Epetra_Vector solution_cell_wghost(mesh_->cell_map(true));
//   solution_cell_wghost.Import(*solution_cell, cell_importer, Insert);
// #else
//   Epetra_Vector& solution_cell_wghost = *solution_cell;
// #endif
  
#ifdef HAVE_MPI
  Epetra_Vector solution_cell_wghost(mesh_->cell_map(true));
  FS_->CopyMasterCell2GhostCell(solution_cells, solution_cell_wghost);
#else
  Epetra_Vector& solution_cell_wghost = *solution_cells;
#endif

  AmanziMesh::Entity_ID_List faces;
  std::vector<double> dp;
  std::vector<int> dirs;

  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  int nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  int nfaces_wghost = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::USED);
  //  std::vector<int> flag(nfaces_wghost, 0);
  //cout.precision(10);

  AmanziMesh::Entity_ID_List cells;

  darcy_mass_flux.PutScalar(0.);;

  for (int c = 0; c < ncells_owned; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];
      if (bc_model[f] == FLOW_BC_FACE_PRESSURE) {
	double value = bc_values[f][0];
	darcy_mass_flux[f] = dirs[n]*Krel_faces[f]*(Trans_faces[f]*(solution_cell_wghost[c] - value) + Grav_term[f]);
	// cout<<"T "<<Krel_faces[f]*Trans_faces[f]<<endl;
	// cout<<"Boundary flux/grav "<<dirs[n]*Krel_faces[f]*(Trans_faces[f]*(solution[c] - value))<<endl;
	// cout<<"Boundary gravity "<< dirs[n]*Krel_faces[f]*Grav_term[f]<<endl;
	// cout<<"Dirichler BC flux "<<darcy_mass_flux[f] <<endl<<endl;
      }
      else if (bc_model[f] == FLOW_BC_FACE_FLUX) {
	double value = bc_values[f][0];
	darcy_mass_flux[f] = value;
      }
      else {
	if (f < nfaces_owned) {
	  mesh_->face_get_cells(f,  AmanziMesh::USED, &cells);
	  if (cells.size() < 2){
	    Errors::Message msg("Flow PK: Matrix_MFD_TPFA. These boundary conditions are not supported by TPFA discratization.");
	    Exceptions::amanzi_throw(msg);
	  }

	  // if ((fabs(Grav_term[f]) > 1e-12 )&&(dirs[n]>0)){
	  //   cout<<"FACE "<<f<<endl;
	  //   cout<<"T "<<Krel_faces[f]*Trans_faces[f]<<" "<<" Krel "<<Krel_faces[f]<<endl;
	  //   //cout<<"Boundary flux/grav "<<dirs[n]*Krel_faces[f]*(Trans_faces[f]*(solution[c] - value))<<endl;
	  //   cout<<"gravity "<< dirs[n]*Krel_faces[f]*Grav_term[f]<<endl;
	  //   cout<<"Grav_term "<<Grav_term[f]<<endl;
	  //   //if (f == 49) {
	  //   if (dirs[n]>0) {
	  //     cout<<"cells: "<<c<<" "<<c+1<<" --- "<<solution[c]<<" "<<solution[c+1]<<endl;
	  //     cout<<"Diff "<< Krel_faces[f]*Trans_faces[f]*(solution[c] - solution[c+1])<<endl;
	  //   }
	  //   //}
	  // }


	  // double s = Trans_faces[f]*solution[c];
	  // darcy_mass_flux[f] += s * dirs[n];
	  // if (cells[0] == c) darcy_mass_flux[f] += dirs[n]*Grav_term[f]*0.5;
	  // else darcy_mass_flux[f] -= dirs[n]*Grav_term[f]*0.5;  
	  // darcy_mass_flux[f] *= Krel_faces[f];
	  // //if (f==49) cout<<"f=49 "<<darcy_mass_flux[f]<<endl;
	  // //if (f==49) cout<<"f=49 "<<darcy_mass_flux[f]<<endl;
	  // //exit(0);
	  if (dirs[n] > 0){
	    if (c == cells[0]){
	      darcy_mass_flux[f] = Trans_faces[f]*(solution_cell_wghost[cells[0]] - solution_cell_wghost[cells[1]]) + Grav_term[f];
	    }
	    else {
	      darcy_mass_flux[f] = Trans_faces[f]*(solution_cell_wghost[cells[1]] - solution_cell_wghost[cells[0]]) + Grav_term[f];
	    }	    
	    darcy_mass_flux[f] *= Krel_faces[f];
	  }
	}
      }
    }
  }

  //cout<<"Darcy\n"<<darcy_mass_flux<<endl;


}



}  // namespace AmanziFlow
}  // namespace Amanzi


