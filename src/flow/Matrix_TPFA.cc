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

#include "FlowDefs.hh"
#include "Matrix_TPFA.hh"

#include "LinearOperatorFactory.hh"

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
Matrix_TPFA::Matrix_TPFA(Teuchos::RCP<State> S,
                         Teuchos::RCP<RelativePermeability> rel_perm)
    : Matrix(S, rel_perm)
{
  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED); 
  Acc_cells_.resize(ncells);
  Fc_cells_.resize(ncells);
}


/* ******************************************************************
* Constructor.                                           
****************************************************************** */
void Matrix_TPFA::Init() 
{
  const Epetra_BlockMap& fmap_wghost = mesh_->face_map(true);
  Transmis_faces = Teuchos::rcp(new Epetra_Vector(fmap_wghost));
  Grav_term_faces = Teuchos::rcp(new Epetra_Vector(fmap_wghost));
  ComputeTransmissibilities(*Transmis_faces, *Grav_term_faces);
}


/* ******************************************************************
* Initialize Trilinos matrices. It must be called only once. 
* If matrix is non-symmetric, we generate transpose of the matrix 
* block Afc to reuse cf_graph; otherwise, pointer Afc = Acf.   
****************************************************************** */
void Matrix_TPFA::SymbolicAssemble()
{
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
* Add gravity fluxes to RHS of TPFA approximation                                            
****************************************************************** */
void Matrix_TPFA::AddGravityFluxes(const Epetra_Vector& Krel_faces, 
                                   const Epetra_Vector& Grav_term)
{
  /* TODO
  AmanziMesh::Entity_ID_List cells;
  std::vector<int> dirs;
  Epetra_MultiVector& rhs_cells = *rhs_->ViewComponent("cell");

  int nfaces_wghost = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::USED);

  for (int f = 0; f < nfaces_wghost; f++) {
    if (bc_model[f] == FLOW_BC_FACE_FLUX) continue;
    mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
    int ncells = cells.size();
    for (int i = 0; i < ncells; i++){
      int c = cells[i];
      if (c >= ncells_owned) continue;
      rhs_cells[0][c] -= pow(-1.0, i)*Grav_term[f]*Krel_faces[f];  
    }
  }
  */
}


/* ******************************************************************
* 
****************************************************************** */
void Matrix_TPFA::Assemble()
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
    (*Spp_).SumIntoGlobalValues(mcells, cells_GID, Spp_local.values());
  }

  Epetra_MultiVector& rhs_cells = *rhs_->ViewComponent("cell", true);

  for (int c = 0; c <= ncells_owned; c++){
    cells_GID[0] = cmap_wghost.GID(c);
    int mcells = 1;
    Spp_local(0,0) = Acc_cells_[c];
    (*Spp_).SumIntoGlobalValues(mcells, cells_GID, Spp_local.values());
    rhs_cells[0][c] += Fc_cells_[c];
  }

  Spp_->GlobalAssemble();
}


/* ******************************************************************
* Computation of the part of the jacobian which depends on
* analytical derivatives of relative permeabilities
****************************************************************** */
void Matrix_TPFA::AnalyticJacobian(
   const Epetra_Vector& solution, 
   std::vector<int>& bc_models, std::vector<bc_tuple>& bc_values)
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

  /* TODO
  FS_->CopyMasterCell2GhostCell(solution, pres_gh);
  */

  const Epetra_MultiVector& Krel_cells = *rel_perm_->Krel().ViewComponent("cell");
  Epetra_MultiVector& Krel_faces = *rel_perm_->Krel().ViewComponent("face", true);

  const Epetra_MultiVector& dKdP_cells = *rel_perm_->dKdP().ViewComponent("cell");
  Epetra_MultiVector& dKdP_faces = *rel_perm_->dKdP().ViewComponent("face", true);
  int method = rel_perm_->method();

  for (int f = 0; f < nfaces_owned; f++){
    mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
    int mcells = cells.size();
    Teuchos::SerialDenseMatrix<int, double> Jpp(mcells, mcells);
    //AmanziGeometry::Point face_cntr = mesh_->face_centroid(f);

    for (int n = 0; n < mcells; n++) {
      cells_LID[n] = cells[n];
      cells_GID[n] = cmap_wghost.GID(cells_LID[n]);
      pres[n] = pres_gh[cells_LID[n]];
      dk_dp[n] = dKdP_cells[0][cells_LID[n]];
    }

    if (mcells == 1) {
      dk_dp[0] = dKdP_faces[0][f];
    }

    const AmanziGeometry::Point& normal = mesh_->face_normal(f, false, cells[0]);

    ComputeJacobianLocal(mcells, f, method, bc_models, bc_values, pres, dk_dp, Jpp);

    (*Spp_).SumIntoGlobalValues(mcells, cells_GID, Jpp.values());
  }

  Spp_->GlobalAssemble();
}


/* ******************************************************************
* Computation of a local submatrix of 
* Analytical Jacobian (nonlinear part) on a particular face.
****************************************************************** */
void Matrix_TPFA::ComputeJacobianLocal(int mcells,
                                       int face_id,
                                       int Krel_method,
                                       std::vector<int>& bc_models,
                                       std::vector<bc_tuple>& bc_values,
                                       double *pres,
                                       double *dk_dp_cell,
                                       Teuchos::SerialDenseMatrix<int, double>& Jpp)
{
  double dKrel_dp[2];
  double rho_w = *S_->GetScalarData("fluid_density");

  double dpres;
  if (mcells == 2) {
    dpres = pres[0] - pres[1];  // + grn;
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

    Jpp(0, 0) = ((*trans_on_faces_)[face_id]*dpres + (*grav_on_faces_)[face_id])*dKrel_dp[0];
    Jpp(0, 1) = ((*trans_on_faces_)[face_id]*dpres + (*grav_on_faces_)[face_id])*dKrel_dp[1];
    Jpp(1, 0) = -Jpp(0, 0);
    Jpp(1, 1) = -Jpp(0, 1);

  } else if (mcells == 1) {
    if (bc_models[face_id] == FLOW_BC_FACE_PRESSURE) {                   
      pres[1] = bc_values[face_id][0];

      dpres = pres[0] - pres[1];  // + grn;
      Jpp(0,0) = ((*trans_on_faces_)[face_id]*dpres + (*grav_on_faces_)[face_id]) * dk_dp_cell[0];
    } else {
      Jpp(0,0) = 0.0;
    }
  }
}


/* ******************************************************************
* Compute transmissibilities on faces 
****************************************************************** */
void Matrix_TPFA::ComputeTransmissibilities(Epetra_Vector& Trans_faces, Epetra_Vector& grav_faces)
{
  Trans_faces.PutScalar(0.0);

  double rho_ = *S_->GetScalarData("fluid_density");
  double mu_ = *S_->GetScalarData("fluid_viscosity");
  const Epetra_Vector& gravity_ = *S_->GetConstantVectorData("gravity");

  int dim = mesh_->space_dimension();
  AmanziGeometry::Point gravity(dim);
  for (int k = 0; k < dim; k++) gravity[k] = gravity_[k] * rho_;

  AmanziMesh::Entity_ID_List faces;
  AmanziMesh::Entity_ID_List cells;
  AmanziGeometry::Point a_dist;
  double h[2], perm[2], perm_test[2], h_test[2];
  double trans_f;

  std::vector<int> dirs;
  int nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);

  for (int f = 0; f < nfaces_owned; f++) {
    mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
    int ncells = cells.size();
    const AmanziGeometry::Point& normal = mesh_->face_normal(f, false, cells[0]);
    const AmanziGeometry::Point& face_centr = mesh_->face_centroid(f);
    double area = mesh_->face_area(f);

    if (ncells == 2) {
      a_dist = mesh_->cell_centroid(cells[1]) - mesh_->cell_centroid(cells[0]);
    } else if (ncells == 1) {    
      a_dist = face_centr - mesh_->cell_centroid(cells[0]);
    } 

    a_dist *= 1./norm(a_dist);

    for (int i=0; i<ncells; i++) {
      h[i] = norm(face_centr - mesh_->cell_centroid(cells[i]));
      perm[i] = (rho_/mu_) * ((K[cells[i]] * normal) * normal) / area;

      perm_test[i] = (rho_/mu_) * ((K[cells[i]] * normal) * a_dist);
      h_test[i] = pow(-1.0, i)*((face_centr - mesh_->cell_centroid(cells[i]))*normal) / area;
    }

    double factor, grav;
    grav = (gravity * normal) / area;

    if (ncells == 2){
      factor = (perm[0]*perm[1]) / (h[0]*perm[1] + h[1]*perm[0]);
      grav *= (h[0] + h[1]);
    } else if (ncells == 1) {    
      factor = perm[0] / h[0];
      grav *= h[0];
    } 

    trans_f = 0.0;
    for (int i = 0; i < ncells; i++) {
      trans_f += h_test[i] / perm[i];
    }

    trans_f = 1.0 / trans_f;

    Trans_faces[f] = trans_f;
    grav_faces[f] = Trans_faces[f] * grav;
  }

  /* TODO 
  FS->CopyMasterFace2GhostFace(Trans_faces);
  FS->CopyMasterFace2GhostFace(grav_faces);
  */
}
/* ******************************************************************
* Parallel matvec product Spp * Xc.                                              
****************************************************************** */
int Matrix_TPFA::Apply(const CompositeVector& X, CompositeVector& Y) const
{
  const Epetra_MultiVector& Xc = *X.ViewComponent("cell");
  const Epetra_MultiVector& Xf = *X.ViewComponent("face");

  Epetra_MultiVector& Yc = *Y.ViewComponent("cell");
  Epetra_MultiVector& Yf = *Y.ViewComponent("face");

  int ierr = (*Spp_).Multiply(false, Xc, Yc);

  if (ierr) {
    Errors::Message msg("Matrix_TPFA::Apply has failed to calculate y = A*x.");
    Exceptions::amanzi_throw(msg);
  }

  Yf.PutScalar(0.0);

  return 0;
}


/* ******************************************************************
* The OWNED cell-based and face-based d.o.f. are packed together into 
* the X and Y Epetra vectors, with the cell-based in the first part.                                           
****************************************************************** */
int Matrix_TPFA::ApplyInverse(const CompositeVector& X, CompositeVector& Y) const
{
  Teuchos::ParameterList plist;
  Teuchos::ParameterList& slist = plist.sublist("gmres");
  slist.set<string>("iterative method", "gmres");
  slist.set<double>("error tolerance", 1e-8 );
  slist.set<int>("maximum number of iterations", 100);
  Teuchos::ParameterList& vlist = slist.sublist("VerboseObject");
  vlist.set("Verbosity Level", "low");

  AmanziSolvers::LinearOperatorFactory<FlowMatrix, CompositeVector, CompositeVectorSpace> factory;
  Teuchos::RCP<AmanziSolvers::LinearOperator<FlowMatrix, CompositeVector, CompositeVectorSpace> > 
      solver = factory.Create("gmres", plist, preconditioner_, preconditioner_);
   
  solver->ApplyInverse(X, Y);

  return 0;
}


/* ******************************************************************
* 
****************************************************************** */
void Matrix_TPFA::ApplyBoundaryConditions(std::vector<int>& bc_model, 
                                          std::vector<bc_tuple>& bc_values){

  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  rhs_->PutScalar(0.0);
  Epetra_MultiVector& rhs_cells = *rhs_->ViewComponent("cell");

  for (int c = 0; c < ncells; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];
      double value = bc_values[f][0];

      if (bc_model[f] == FLOW_BC_FACE_PRESSURE) {
	rhs_cells[0][c] += value * (*trans_on_faces_)[f] * (*Krel_faces_)[f];
      } else if (bc_model[f] == FLOW_BC_FACE_FLUX) {
        rhs_cells[0][c] -= value * mesh_->face_area(f);
	(*trans_on_faces_)[f] = 0.0;
	(*grav_on_faces_)[f] = 0.0;
      } else if (bc_model[f] == FLOW_BC_FACE_MIXED) {
	Errors::Message msg;
	msg << "Mixed boundary conditions are not supported in TPFA mode\n";
	Exceptions::amanzi_throw(msg);
      }
    }
  }
}


/* ******************************************************************
* Linear algebra operations with matrices: r = A * u - f                                                 
****************************************************************** */
double Matrix_TPFA::ComputeNegativeResidual(const CompositeVector& u, CompositeVector& r)
{
  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  int nfaces_wghost = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::USED);

  AmanziMesh::Entity_ID_List faces;
  AmanziMesh::Entity_ID_List cells;
  std::vector<int> dirs;

  r.PutScalar(0.0);
  Epetra_MultiVector& rc = *r.ViewComponent("cell");

  u.ScatterMasterToGhosted("cell");
  const Epetra_MultiVector& uc = *u.ViewComponent("cell", true);

  for (int f = 0; f < nfaces_wghost; f++) {
    mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
    int ncells = cells.size();
    if (ncells > 1) {
      int sign(1);
      for (int i = 0; i < ncells; i++) {
	int c = cells[i];
	if (c >= ncells_owned) continue;
	rc[0][c] += sign*(*Krel_faces_)[f]*(*trans_on_faces_)[f]*(uc[0][cells[0]] - uc[0][cells[1]]);  
        sign = -sign;
      }
    } else if (ncells == 1) {
      int c = cells[0];
      rc[0][c] += (*Krel_faces_)[f] * (*trans_on_faces_)[f] * uc[0][c];
    }							
  } 
  
  Epetra_MultiVector& rhs_cells = *rhs_->ViewComponent("cell");
  for (int c = 0; c < ncells_owned; c++) {    
    rc[0][c] -= rhs_cells[0][c];
  }
    
  double norm_residual;
  r.Norm2(&norm_residual);
  return norm_residual;
}


/* ******************************************************************
* 
****************************************************************** */
void Matrix_TPFA::DeriveMassFlux(
    const CompositeVector& solution, CompositeVector& darcy_mass_flux,
    std::vector<int>& bc_model, std::vector<bc_tuple>& bc_values)
{
  solution.ScatterMasterToGhosted("cell");
  const Epetra_MultiVector& p = *solution.ViewComponent("cell", true);
  Epetra_MultiVector& flux = *darcy_mass_flux.ViewComponent("face");

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

      if (bc_model[f] == FLOW_BC_FACE_PRESSURE) {
	double value = bc_values[f][0];
	flux[0][f] = dirs[n]*(*Krel_faces_)[f]*((*trans_on_faces_)[f]*(p[0][c] - value) + (*grav_on_faces_)[f]);
      } else if (bc_model[f] == FLOW_BC_FACE_FLUX) {
	double value = bc_values[f][0];
	double area = mesh_->face_area(f);
	flux[0][f] = value*area;
      } else {
	if (f < nfaces_owned && !flag[f]) {
	//if (f < nfaces_owned) {
	  
	  mesh_->face_get_cells(f,  AmanziMesh::USED, &cells);
	  if (cells.size() <= 1 ){
	    Errors::Message msg("Flow PK: Matrix_TPFA. These boundary conditions are not supported by TPFA discratization.");
	    Exceptions::amanzi_throw(msg);
	  }

          if (c == cells[0]){
            flux[0][f] = dirs[n]*(*trans_on_faces_)[f]*(p[0][cells[0]] - p[0][cells[1]]) + (*grav_on_faces_)[f];
          } else {
            flux[0][f] = dirs[n]*(*trans_on_faces_)[f]*(p[0][cells[1]] - p[0][cells[0]]) + (*grav_on_faces_)[f];
          }	    
          flux[0][f] *= (*Krel_faces_)[f];
          flag[f] = 1;
	}
      }
    }
  }
}

}  // namespace AmanziFlow
}  // namespace Amanzi


