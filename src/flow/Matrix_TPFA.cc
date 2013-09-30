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
// #include "AztecOO.h"
#include "mfd3d_diffusion.hh"

#include "FlowDefs.hh"
#include "Matrix_TPFA.hh"


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
Matrix_MFD_TPFA::Matrix_MFD_TPFA(Teuchos::RCP<Flow_State> FS, Teuchos::RCP<const Epetra_Map> map) 
   :  Matrix_MFD(FS, map)
{
  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED); 
  Acc_cells_.resize(ncells);
  Fc_cells_.resize(ncells);
}

  Matrix_MFD_TPFA::Matrix_MFD_TPFA(Teuchos::RCP<Flow_State> FS, Teuchos::RCP<const Epetra_Map> map,
				   Teuchos::RCP<Epetra_Vector> Krel_faces,
				   Teuchos::RCP<Epetra_Vector> Trans_faces,
				   Teuchos::RCP<Epetra_Vector> Grav_faces) 
   :  Matrix_MFD(FS, map)
{
  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED); 
  Acc_cells_.resize(ncells);
  Fc_cells_.resize(ncells);

  Krel_faces_ = Krel_faces;
  trans_on_faces_ = Trans_faces;
  grav_on_faces_ = Grav_faces;
}



/* ******************************************************************
* Calculate elemental stiffness matrices.                                            
****************************************************************** */
void Matrix_MFD_TPFA::CreateMFDstiffnessMatrices(RelativePermeability& rel_perm)
{

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

void Matrix_MFD_TPFA::AssembleGlobalMatrices()
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
	tij = (*trans_on_faces_)[f] * (*Krel_faces_)[f];
	if (i==j) Spp_local(i,j) = tij;
	else Spp_local(i,j) = -tij;
      }
    }
    //cout<<"face "<<f<<": "<<tij<<" ukvr "<<Krel_faces[f]<<endl;
    (*Spp_).SumIntoGlobalValues(mcells, cells_GID, Spp_local.values());

  }

  for (int c = 0; c <= ncells_owned; c++){
    cells_GID[0] = cmap_wghost.GID(c);
    int mcells = 1;
    Spp_local(0,0) = Acc_cells_[c];
    (*Spp_).SumIntoGlobalValues(mcells, cells_GID, Spp_local.values());
    (*rhs_cells_)[c] += Fc_cells_[c];
  }

  Spp_->GlobalAssemble();

    
}


/* ******************************************************************
* Assembles preconditioner. It has same set of parameters as matrix.
****************************************************************** */
void Matrix_MFD_TPFA::AssembleSchurComplement(std::vector<int>& bc_model, std::vector<bc_tuple>& bc_values)
{
  AssembleGlobalMatrices();
}

 
/* ******************************************************************
* Computation of the part of the jacobian which depends on
* analytical derivatives of relative permeabilities
****************************************************************** */
void Matrix_MFD_TPFA::AnalyticJacobian(
   const Epetra_Vector& solution, 
   std::vector<int>& bc_models, std::vector<bc_tuple>& bc_values,
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

  FS_->CopyMasterCell2GhostCell(solution, pres_gh);

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
      pres[n] = pres_gh[cells_LID[n]];
      dk_dp[n] = dKdP_cells[cells_LID[n]];
    }

    if (mcells == 1) {
      dk_dp[0] = dKdP_faces[f];
    }

    const AmanziGeometry::Point& normal = mesh_->face_normal(f, false, cells[0]);


    ComputeJacobianLocal(mcells, f, method, bc_models, bc_values, pres, dk_dp, Jpp);

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
                                           int Krel_method,
                                           std::vector<int>& bc_models,
                                           std::vector<bc_tuple>& bc_values,
                                           double *pres,
                                           double *dk_dp_cell,
                                           Teuchos::SerialDenseMatrix<int, double>& Jpp)
{
  // double K[2];
  double dKrel_dp[2];



  double rho_w = FS_->ref_fluid_density();


  double dpres;
  if (mcells == 2) {
    dpres = pres[0] - pres[1];// + grn;
		// cout<<"pres[0] "<<pres[0]<<" pres[1] "<<pres[1]<<" grv "<<grn<<endl;
    if (Krel_method == FLOW_RELATIVE_PERM_UPWIND_GRAVITY) {  // Define K and Krel_faces
      double cos_angle = (*grav_on_faces_)[face_id]/(rho_w*mesh_->face_area(face_id));
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
      if ((*trans_on_faces_)[face_id]*dpres + (*grav_on_faces_)[face_id] > FLOW_RELATIVE_PERM_TOLERANCE) {  // Upwind
        dKrel_dp[0] = dk_dp_cell[0];
        dKrel_dp[1] = 0.0;
    } else if ((*trans_on_faces_)[face_id]*dpres + (*grav_on_faces_)[face_id] < -FLOW_RELATIVE_PERM_TOLERANCE) {  // Upwind
        dKrel_dp[0] = 0.0;
        dKrel_dp[1] = dk_dp_cell[1];
    } else if (fabs((*trans_on_faces_)[face_id]*dpres + (*grav_on_faces_)[face_id]) < FLOW_RELATIVE_PERM_TOLERANCE) {  // Upwind
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

    Jpp(0, 0) = ((*trans_on_faces_)[face_id]*dpres + (*grav_on_faces_)[face_id])*dKrel_dp[0];
    Jpp(0, 1) = ((*trans_on_faces_)[face_id]*dpres + (*grav_on_faces_)[face_id])*dKrel_dp[1];
    Jpp(1, 0) = -Jpp(0, 0);
    Jpp(1, 1) = -Jpp(0, 1);

  } else if (mcells == 1) {
    if (bc_models[face_id] == FLOW_BC_FACE_PRESSURE) {                   
      pres[1] = bc_values[face_id][0];

      dpres = pres[0] - pres[1];// + grn;
      //cout<<"dk_dp_cell[0] "<<dk_dp_cell[0]<<endl;
      Jpp(0,0) = ((*trans_on_faces_)[face_id]*dpres + (*grav_on_faces_)[face_id]) * dk_dp_cell[0];
    } else {
      Jpp(0,0) = 0.0;
    }
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
  preconditioner_->ApplyInverse(Xc, Tc);
  Yc = Tc;

  if (ierr) {
    Errors::Message msg("Matrix_MFD_TPFA::ApplyInverse has failed in calculating y = inv(A)*x.");
    Exceptions::amanzi_throw(msg);
  }

  Yf = Xf;

  delete [] fvec_ptrs;
  return 0;
}


void  Matrix_MFD_TPFA::ApplyBoundaryConditions(std::vector<int>& bc_model, 
					       std::vector<bc_tuple>& bc_values){

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
	(*rhs_cells_)[c] += value * (*trans_on_faces_)[f] * (*Krel_faces_)[f];
      } else if (bc_model[f] == FLOW_BC_FACE_FLUX) {
        (*rhs_cells_)[c] -= value * mesh_->face_area(f);
	(*trans_on_faces_)[f] = 0.0;
	(*grav_on_faces_)[f] = 0.0;
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
	residual[c] += pow(-1.0, i)*(*Krel_faces_)[f]*(*trans_on_faces_)[f]*(sol_gh[cells[0]] - sol_gh[cells[1]]);  
      }
    }
    else if (ncells == 1){
      int c = cells[0];
      residual[c] += (*Krel_faces_)[f] * (*trans_on_faces_)[f] * sol_gh[c];
    }							
  } 
  


  for (int c = 0; c < ncells_owned; c++) {    
    residual[c] -= (*rhs_cells_)[c];
  }
    
  double norm_residual;
  residual.Norm2(&norm_residual);
  return norm_residual;
}

  void Matrix_MFD_TPFA::DeriveDarcyMassFlux(const Epetra_Vector& solution_cells,		
					    const Epetra_Import& face_importer,
					    std::vector<int>& bc_model, 
					    std::vector<bc_tuple>& bc_values,
					    Epetra_Vector& darcy_mass_flux)
{
  
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

  AmanziMesh::Entity_ID_List cells;
  std::vector<int> flag(nfaces_wghost, 0);


  for (int c = 0; c < ncells_owned; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];


      //int GID = face_wghost.GID(f);


      if (bc_model[f] == FLOW_BC_FACE_PRESSURE) {
	double value = bc_values[f][0];
	darcy_mass_flux[f] = dirs[n]*(*Krel_faces_)[f]*((*trans_on_faces_)[f]*(solution_cell_wghost[c] - value) + (*grav_on_faces_)[f]);
      }
      else if (bc_model[f] == FLOW_BC_FACE_FLUX) {
	double value = bc_values[f][0];
	double area = mesh_->face_area(f);
	darcy_mass_flux[f] = value*area ;
      }
      else {
	if (f < nfaces_owned && !flag[f]) {
	//if (f < nfaces_owned) {
	  
	  mesh_->face_get_cells(f,  AmanziMesh::USED, &cells);
	  if (cells.size() <= 1 ){
	    Errors::Message msg("Flow PK: Matrix_MFD_TPFA. These boundary conditions are not supported by TPFA discratization.");
	    Exceptions::amanzi_throw(msg);
	  }

	    if (c == cells[0]){
	      darcy_mass_flux[f] = dirs[n]*(*trans_on_faces_)[f]*(solution_cell_wghost[cells[0]] - solution_cell_wghost[cells[1]]) + (*grav_on_faces_)[f];
	    }
	    else {
	      darcy_mass_flux[f] = dirs[n]*(*trans_on_faces_)[f]*(solution_cell_wghost[cells[1]] - solution_cell_wghost[cells[0]]) + (*grav_on_faces_)[f];
	    }	    
	    darcy_mass_flux[f] *= (*Krel_faces_)[f];
	    flag[f] = 1;
	 
	}
      }
    }
  }

  //cout<<"Darcy\n"<<darcy_mass_flux<<endl;


}



}  // namespace AmanziFlow
}  // namespace Amanzi


