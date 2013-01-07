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

#include "Flow_constants.hpp"
#include "Matrix_MFD_TPFA.hpp"

namespace Amanzi {
namespace AmanziFlow {

template<class T>
int FindPosition(const std::vector<T>& v, const T& value) {
  for (int i = 0; i < v.size(); i++)
     if (v[i] == value) return i;
  return -1;
}


/* ******************************************************************
* Calculate elemental stiffness matrices.                                            
****************************************************************** */
void Matrix_MFD_TPFA::CreateMFDstiffnessMatrices(Epetra_Vector& Krel_cells,
                                                 Epetra_Vector& Krel_faces)
{
  int dim = mesh_->space_dimension();
  WhetStone::MFD3D mfd(mesh_);
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  Aff_cells_.clear();
  Afc_cells_.clear();
  Acf_cells_.clear();
  Acc_cells_.clear();

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
void Matrix_MFD_TPFA::AssembleGlobalMatrices()
{
  Matrix_MFD::AssembleGlobalMatrices();

  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;
  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);

  Dff_->PutScalar(0.0);
  for (int c = 0; c < ncells_owned; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int mfaces = faces.size();

    for (int n = 0; n < mfaces; n++) {
      int f = faces[n];
      (*Dff_)[f] += Aff_cells_[c](n, n);
    }
  }
  FS->CombineGhostFace2MasterFace(*Dff_, Add);

  // convert right-hand side to a cell-based vector
  const Epetra_Map& cmap = mesh_->cell_map(false);
  Epetra_Vector Tc(cmap);
  int nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);

  for (int f = 0; f < nfaces_owned; f++) (*rhs_faces_)[f] /= (*Dff_)[f];
  (*Acf_).Multiply(false, *rhs_faces_, Tc);
  for (int c = 0; c < ncells_owned; c++) (*rhs_cells_)[c] -= Tc[c];

  rhs_faces_->PutScalar(0.0);

  // create a with-ghost copy of Acc
  const Epetra_Map& cmap_wghost = mesh_->cell_map(true);
  Epetra_Vector Dcc(cmap_wghost);

  for (int c = 0; c < ncells_owned; c++) Dcc[c] = (*Acc_)[c];
  FS->CopyMasterCell2GhostCell(Dcc);

  AmanziMesh::Entity_ID_List cells;
  int cells_GID[2];
  double Acf_copy[2];

  // create auxiliaty with-ghost copy of Acf_cells
  const Epetra_Map& fmap_wghost = mesh_->face_map(true);
  Epetra_Vector Acf_parallel(fmap_wghost);

  for (int c = 0; c < ncells_owned; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int mfaces = faces.size();

    for (int i = 0; i < mfaces; i++) {
      int f = faces[i];
      if (f >= nfaces_owned) Acf_parallel[f] = Acf_cells_[c][i];
    }
  }
  FS->CombineGhostFace2MasterFace(Acf_parallel, Add);

  // populate the global matrix
  Spp_->PutScalar(0.0);
  for (AmanziMesh::Entity_ID f = 0; f < nfaces_owned; f++) {
    mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
    int mcells = cells.size();

    // populate face-based matrix.
    Teuchos::SerialDenseMatrix<int, double> Bpp(mcells, mcells);
    for (int n = 0; n < mcells; n++) {
      int c = cells[n];
      cells_GID[n] = cmap_wghost.GID(c);

      mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
      int i = FindPosition<AmanziMesh::Entity_ID>(faces, f);
      Bpp(n, n) = Dcc[c] / faces.size();
      if (c < ncells_owned) {
        int i = FindPosition<AmanziMesh::Entity_ID>(faces, f);
        Acf_copy[n] = Acf_cells_[c][i];
      } else {
        Acf_copy[n] = Acf_parallel[f];
      }
    }

    for (int n = 0; n < mcells; n++) {
      for (int m = n; m < mcells; m++) {
        Bpp(n, m) -= Acf_copy[n] * Acf_copy[m] / (*Dff_)[f];
        Bpp(m, n) = Bpp(n, m);
      }
    }

    (*Spp_).SumIntoGlobalValues(mcells, cells_GID, Bpp.values());
  }
  (*Spp_).GlobalAssemble();
}


/* ******************************************************************
* Computation of the part of the jacobian which depends on
* analytical derivatives of relative permeabilities
****************************************************************** */
void Matrix_MFD_TPFA::AnalyticJacobian(const Epetra_Vector& solution,
                                       int dim,
                                       int Krel_method,
                                       std::vector<int>& bc_models,
                                       Epetra_Vector& Krel_cells,
                                       Epetra_Vector& dK_dP_cells,
                                       Epetra_Vector& Krel_faces,
                                       Epetra_Vector& dK_dP_faces)
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

     //     cout<<"Before Analytic\n";
     //     std::cout<<(*Spp_)<<endl;

     Epetra_Vector pres_gh(cmap_wghost);

//      Epetra_Vector pressure_vec(solution);
     FS->CopyMasterCell2GhostCell(solution, pres_gh);

     Epetra_Vector perm_vert_gh(cmap_wghost);
     Epetra_Vector& perm_vert_vec = FS->ref_vertical_permeability();
     FS->CopyMasterCell2GhostCell(perm_vert_vec, perm_vert_gh);

     Epetra_Vector perm_horz_gh(cmap_wghost);
     Epetra_Vector& perm_horz_vec = FS->ref_horizontal_permeability();
     FS->CopyMasterCell2GhostCell(perm_horz_vec, perm_horz_gh);

     Epetra_Vector Krel_gh(cmap_wghost);
     FS->CopyMasterCell2GhostCell(Krel_cells, Krel_gh);

     Epetra_Vector dK_dP_gh(cmap_wghost);
     FS->CopyMasterCell2GhostCell(dK_dP_cells, dK_dP_gh);



     for (int f = 0; f < nfaces_owned; f++){
                mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
                int mcells = cells.size();
                Teuchos::SerialDenseMatrix<int, double> Jpp(mcells, mcells);
                AmanziGeometry::Point face_cntr = mesh_->face_centroid(f);

                for (int n = 0; n < mcells; n++) {
                        cells_LID[n] = cells[n];
                        cells_GID[n] = cmap_wghost.GID(cells_LID[n]);
                        perm_abs_vert[n] = perm_vert_gh[cells_LID[n]];
                        perm_abs_horz[n] = perm_horz_gh[cells_LID[n]];
                        pres[n] = pres_gh[cells_LID[n]];
                        k_rel[n] = Krel_gh[cells_LID[n]];
                        dk_dp[n] = dK_dP_gh[cells_LID[n]];
                        cntr_cell[n] = mesh_->cell_centroid(cells_LID[n]);
                }
                if (mcells == 2){
                  dist = norm(cntr_cell[0] - cntr_cell[1]);
                }
                else if (mcells == 1)
                {
                  dist = norm(cntr_cell[0] - face_cntr);
                  k_rel[0] = Krel_faces[f];
                  dk_dp[0] = dK_dP_faces[f];
                }

                AmanziGeometry::Point normal = mesh_->face_normal(f, false, cells_LID[0]);
                normal *= 1./ mesh_->face_area(f);

                ComputeJacobianLocal(mcells, f, dim, Krel_method, bc_models, dist,  pres,
                                     perm_abs_vert, perm_abs_horz, k_rel, dk_dp, normal, Jpp);


                (*Spp_).SumIntoGlobalValues(mcells, cells_GID, Jpp.values());
     }
     Spp_->GlobalAssemble();
 //     cout << "AnalyticJacobian Spp\n";
 //             std::cout<<(*Spp_)<<endl;
 //               cout << "Matrix_MFD_TPFA:: ComputeJacobian\n";
 //     exit(0);
}


/* ******************************************************************
* Computation of a local submatrix of 
* Analytical Jacobian (nonlinear part) on a particular face.
****************************************************************** */
void Matrix_MFD_TPFA::ComputeJacobianLocal(int mcells,
                                           int face_id,
                                           int dim,
                                           int Krel_method,
                                           std::vector<int>& bc_models,
                                           double dist,
                                           double *pres,
                                           double *perm_abs_vert,
                                           double *perm_abs_horz,
                                           double *k_rel,
                                           double *dk_dp_cell,
                                           AmanziGeometry::Point& normal,
                                           Teuchos::SerialDenseMatrix<int, double>& Jpp
                                          ){

        double K[2];
        double dKrel_dp[2];

        double rho_w = FS->ref_fluid_density();


        for (int c = 0; c < mcells; c++) {
                K[c] = 0.;
                for (int i = 0;i < dim-1; i++)  K[c] += normal[i]*normal[i]*perm_abs_horz[c];
                K[c] += normal[dim-1]*normal[dim-1]*perm_abs_vert[c];
                K[c] /= FS->ref_fluid_viscosity();
        }

        double Kabs_dist;

        if (mcells == 1) {
                Kabs_dist = K[0]/dist;
        }
        else
        {
                Kabs_dist = 2*K[0]*K[1]/(dist*(K[0] + K[1]));
        }

        double rho = FS->ref_fluid_density();
        AmanziGeometry::Point gravity(dim);
        for (int i = 0; i < dim; i++) gravity[i] = (*(FS->gravity()))[i] * rho;

        double grn = dist*gravity*normal;

//         cout<<"grn "<<grn<<endl;
//         for (int i = 0; i < dim; i++) cout<<"dist "<<(*(FS->gravity()))[i]<<endl;

        double dphi = pres[0] - pres[1] + grn;

        if (mcells == 2){
          //                cout<<"pres[0] "<<pres[0]<<" pres[1] "<<pres[1]<<" grv "<<grn<<endl;
                if (Krel_method == FLOW_RELATIVE_PERM_UPWIND_GRAVITY) {  // Define K and Krel_faces
                        if (grn*Kabs_dist > FLOW_RELATIVE_PERM_TOLERANCE){    // Upwind
                                dKrel_dp[0] = dk_dp_cell[0];
                                dKrel_dp[1] = 0.;
                        }
                        else if (grn*Kabs_dist < -FLOW_RELATIVE_PERM_TOLERANCE){    // Upwind
                                dKrel_dp[0] = 0.;
                                dKrel_dp[1] = dk_dp_cell[1];
                        }
                        else if (fabs(grn*Kabs_dist) < FLOW_RELATIVE_PERM_TOLERANCE){   // Upwind
                                dKrel_dp[0] = 0.5*dk_dp_cell[0];
                                dKrel_dp[1] = 0.5*dk_dp_cell[1];
                        }
                } else if (Krel_method == FLOW_RELATIVE_PERM_UPWIND_DARCY_FLUX) {
                        if (dphi > FLOW_RELATIVE_PERM_TOLERANCE){   // Upwind
                                dKrel_dp[0] = dk_dp_cell[0];
                                dKrel_dp[1] = 0.;
                        }
                        else if (dphi < -FLOW_RELATIVE_PERM_TOLERANCE){   // Upwind
                                dKrel_dp[0] = 0.;
                                dKrel_dp[1] = dk_dp_cell[1];
                        }
                        else if (fabs(dphi) < FLOW_RELATIVE_PERM_TOLERANCE){   // Upwind
                                dKrel_dp[0] = 0.5*dk_dp_cell[0];
                                dKrel_dp[1] = 0.5*dk_dp_cell[1];
                        }
                } else if (Krel_method == FLOW_RELATIVE_PERM_ARITHMETIC_MEAN) {
                        dKrel_dp[0] = 0.5*dk_dp_cell[0];
                        dKrel_dp[1] = 0.5*dk_dp_cell[1];
                }
                // cout<<Kabs_dist*dKrel_dp[0]<<" "<<dphi<<" "<<rho_w*mesh_->face_area(face_id)<<endl;
                // cout<<"dKrel_dp[0] "<<dKrel_dp[0] <<" dKrel_dp[1] "<<dKrel_dp[1]<<endl;

                Jpp(0, 0) = Kabs_dist*dphi*dKrel_dp[0]*rho_w*mesh_->face_area(face_id);
                Jpp(0, 1) = Kabs_dist*dphi*dKrel_dp[1]*rho_w*mesh_->face_area(face_id);
                Jpp(1, 0) = -Jpp(0, 0);
                Jpp(1, 1) = -Jpp(0, 1);

                //                 cout<<"Jpp_local"<<endl;
                 // for (int i=0;i<2;i++){
                 //   for (int j=0;j<2;j++) cout<<Jpp(i,j)<<" ";
                 //   cout<<endl;
                 // }
                 // cout<<endl;
//                 exit(0);
        }
        else if (mcells == 1){
                 if ((bc_models[face_id] == FLOW_BC_FACE_PRESSURE) ||
                         (bc_models[face_id] == FLOW_BC_FACE_PRESSURE_SEEPAGE)) {
                                Jpp(0,0) = -Kabs_dist*dphi*dk_dp_cell[0]*rho_w*mesh_->face_area(face_id);
                 }
                 else
                        Jpp(0,0) = 0.;
        }

}


/* ******************************************************************
* TBW
****************************************************************** */
void Matrix_MFD_TPFA::AddCol2NumJacob(int cell, Epetra_Vector& r){

        WhetStone::MFD3D mfd3d(mesh_);

        AmanziMesh::Entity_ID_List faces;
        std::vector<int> dirs;

        mesh_->cell_get_faces_and_dirs(cell, &faces, &dirs);
        int mfaces = faces.size();
        int indices;
        double values;

        indices = cell;
        values = r[cell];

        NumJac_->ReplaceGlobalValues(cell, 1, &values, &indices);

        for (int i=0; i<mfaces; i++){
                int c_adj = mfd3d.cell_get_face_adj_cell(cell, faces[i]);
        //        cout<<cell<<" - c_adj "<<c_adj<<endl;
                if (c_adj >= 0) {
                    values = r[c_adj];
                    NumJac_->ReplaceGlobalValues(c_adj, 1, &values, &indices);
                }

        }
//         cout<<"NonZeros "<<NumJac_->NumGlobalEntries(cell)<<endl;;

        NumJac_->GlobalAssemble();

 //       std::cout<<"NumJac_ "<<r[cell]<<"\n";
 //       std::cout<<(*NumJac_)<<endl;

}


/* ******************************************************************
* TBW
****************************************************************** */
void Matrix_MFD_TPFA::CompareJacobians(){

        int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);

        double ErrMax = 0;
        double ErrL2 = 0;
        double NormL2 = 0;
        double NormMax = 0;


        for (int cell=0; cell<ncells_owned; cell++){
                int nonzeros = NumJac_->NumGlobalEntries(cell);


                for (int i=0;i<nonzeros;i++){
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
    // hypre_list.set("Solver", PCG);
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

  // Solve the Schur complement system Spp * Yc = Xc. Since AztecOO may
  // use the same memory for X and Y, we introduce auxiliaty vector Tc.
  int ierr = 0;
  Epetra_Vector Tc(cmap);
  if (method_ == FLOW_PRECONDITIONER_TRILINOS_ML) {
    MLprec->ApplyInverse(Xc, Tc);
  } else if (method_ == FLOW_PRECONDITIONER_HYPRE_AMG) {
#ifdef HAVE_HYPRE
    ierr = IfpHypre_Spp_->ApplyInverse(Xc, Tc);
#endif
  } else if (method_ == FLOW_PRECONDITIONER_TRILINOS_BLOCK_ILU) {
    ifp_prec_->ApplyInverse(Xc, Tc);
  }
  Yc = Tc;

  if (ierr) {
    Errors::Message msg("Matrix_MFD_TPFA::ApplyInverse has failed in calculating y = inv(A)*x.");
    Exceptions::amanzi_throw(msg);
  }

  Yf = Xf;

  delete [] fvec_ptrs;
  return 0;
}

}  // namespace AmanziFlow
}  // namespace Amanzi


