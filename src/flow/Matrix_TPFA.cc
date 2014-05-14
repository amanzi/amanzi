/*
  This is the flow component of the Amanzi code. 

  Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
           Daniil Svyatskiy (dasvyat@lanl.gov)
*/

#include <vector>

#include "Epetra_FECrsGraph.h"
#include "mfd3d_diffusion.hh"

#include "FlowDefs.hh"
#include "Matrix_TPFA.hh"

#include "LinearOperatorFactory.hh"
#include "PreconditionerFactory.hh"

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
                         std::vector<WhetStone::Tensor>* K, 
                         Teuchos::RCP<RelativePermeability> rel_perm)
    : Matrix<CompositeVector, CompositeVectorSpace>(S, K, rel_perm)
{
  ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  ncells_wghost = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::USED);

  nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  nfaces_wghost = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::USED);

  Acc_cells_.resize(ncells_owned);
  Fc_cells_.resize(ncells_owned);

  npassed_ = 0;
  nokay_ = 0;
}


/* ******************************************************************
* Constructor.                                           
****************************************************************** */
void Matrix_TPFA::Init() 
{
  const Epetra_BlockMap& fmap_wghost = mesh_->face_map(true);
  transmissibility_ = Teuchos::rcp(new Epetra_Vector(fmap_wghost));
  gravity_term_ = Teuchos::rcp(new Epetra_Vector(fmap_wghost));
  //  faces_dir_ = Teuchos::rcp(new Epetra_Vector(fmap_wghost));
  //face_flag_.assign(nfaces_wghost, 0);

  ComputeTransmissibilities_();
}


/* ******************************************************************
* Initialize Trilinos matrices. It must be called only once. 
* If matrix is non-symmetric, we generate transpose of the matrix 
* block Afc to reuse cf_graph; otherwise, pointer Afc = Acf.   
****************************************************************** */
void Matrix_TPFA::SymbolicAssemble()
{
  if (cvs_.size() == 0) {  // ugly solution (lipnikov@lanl.gov) 
    cvs_.SetMesh(mesh_);
    cvs_.SetGhosted(true);
    cvs_.SetComponent("cell", AmanziMesh::CELL, 1);
  }

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
  Spp_ = Teuchos::rcp(new Epetra_FECrsMatrix(Copy, pp_graph));
  Spp_->GlobalAssemble();

  // create space for right-hand side
  rhs_ = Teuchos::RCP<CompositeVector>(new CompositeVector(cvs_, true));
}


/* ******************************************************************
* Add gravity fluxes to RHS of TPFA approximation                                            
****************************************************************** */
void Matrix_TPFA::AddGravityFluxesRichards(double rho, const AmanziGeometry::Point& gravity, 
                                           std::vector<int>& bc_model)
{
  Epetra_MultiVector& rhs_cells = *rhs_->ViewComponent("cell");
  Epetra_MultiVector& Krel_faces = *rel_perm_->Krel().ViewComponent("face", true);

  AmanziMesh::Entity_ID_List cells, faces;
  std::vector<int> dirs;


  // for (int f = 0; f < nfaces_wghost; f++) {
  //   if (bc_model[f] == FLOW_BC_FACE_FLUX) continue;
  //   mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
  //   int ncells = cells.size();

  //   int sign(1);
  //   for (int i = 0; i < ncells; i++) {
  //     int c = cells[i];

  //     if (c < ncells_owned){
  // 	rhs_cells[0][c] -= sign * (*gravity_term_)[f] * Krel_faces[0][f];  
  //     }
  //     sign = -sign;
  //   }
  // }
  AmanziGeometry::Point face_centr, cell_cntr;
  int sign;
  for (int c=0; c<ncells_owned; c++){
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
 
    int nfaces = faces.size();
    for (int i=0; i<nfaces; i++){
      int f = faces[i];
      if (bc_model[f] == FLOW_BC_FACE_FLUX) continue;
      if (dirs[i] > 0) sign = -1;
      else sign = 1;
      rhs_cells[0][c] += sign * (*gravity_term_)[f] * Krel_faces[0][f];  
    }
   
  }

}


/* ******************************************************************
* 
****************************************************************** */
void Matrix_TPFA::Assemble()
{
  Epetra_MultiVector& Krel_faces = *rel_perm_->Krel().ViewComponent("face", true);

  AmanziMesh::Entity_ID_List faces;

  const Epetra_Map& cmap_wghost = mesh_->cell_map(true);
  AmanziMesh::Entity_ID_List cells;
  //AmanziGeometry::Point face_centr, cell_cntr;

  int cells_LID[2], cells_GID[2];
  Teuchos::SerialDenseMatrix<int, double> Spp_local(2, 2);

  Spp_->PutScalar(0.0);

  for (int f = 0; f < nfaces_owned; f++){
    mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
    int mcells = cells.size();

    for (int n = 0; n < mcells; n++) {
      cells_LID[n] = cells[n];
      cells_GID[n] = cmap_wghost.GID(cells_LID[n]);     
    }

    double tij = (*transmissibility_)[f] * Krel_faces[0][f];

    for (int i=0; i < mcells; i++){
      for (int j=0; j < mcells; j++){
	if (i==j) Spp_local(i,j) = tij;
	else Spp_local(i,j) = -tij;
      }
    }
    (*Spp_).SumIntoGlobalValues(mcells, cells_GID, Spp_local.values());
  }

  Epetra_MultiVector& rhs_cells = *rhs_->ViewComponent("cell", true);

  int mcells = 1;
  for (int c = 0; c < ncells_owned; c++){

    cells_GID[0] = cmap_wghost.GID(c);
    Spp_local(0,0) = Acc_cells_[c];
    (*Spp_).SumIntoGlobalValues(mcells, cells_GID, Spp_local.values());
    rhs_cells[0][c] += Fc_cells_[c];

  }

  Spp_->GlobalAssemble();
}


/* ******************************************************************
* Initialization of the preconditioner                                                 
****************************************************************** */
void Matrix_TPFA::InitPreconditioner(const std::string& prec_name, const Teuchos::ParameterList& prec_list)
{
  AmanziPreconditioners::PreconditionerFactory factory;
  preconditioner_ = factory.Create(prec_name, prec_list);
}


/* ******************************************************************
* Parallel matvec product Spp * Xc.                                              
****************************************************************** */
int Matrix_TPFA::Apply(const CompositeVector& X, CompositeVector& Y) const
{
  const Epetra_MultiVector& Xc = *X.ViewComponent("cell");
  Epetra_MultiVector& Yc = *Y.ViewComponent("cell");

  int ierr = (*Spp_).Multiply(false, Xc, Yc);

  if (ierr) {
    Errors::Message msg("Matrix_TPFA::Apply has failed to calculate y = A*x.");
    Exceptions::amanzi_throw(msg);
  }

  return 0;
}


/* ******************************************************************
* The OWNED cell-based and face-based d.o.f. are packed together into 
* the X and Y Epetra vectors, with the cell-based in the first part.                                           
****************************************************************** */
int Matrix_TPFA::ApplyPreconditioner(const CompositeVector& X, CompositeVector& Y) const
{
  Teuchos::ParameterList plist;
  Teuchos::ParameterList& pre_list = plist.sublist("gmres");
  Teuchos::ParameterList& slist = pre_list.sublist("gmres parameters");

  pre_list.set<string>("iterative method", "gmres");
  slist.set<double>("error tolerance", 1e-7);
  slist.set<int>("maximum number of iterations", 200);
  Teuchos::ParameterList& vlist = slist.sublist("VerboseObject");
  vlist.set("Verbosity Level", "low");

  Teuchos::RCP<const Matrix_TPFA> matrix_tmp = Teuchos::rcp(this, false);

  AmanziSolvers::LinearOperatorFactory<FlowMatrix, CompositeVector, CompositeVectorSpace> factory;
  Teuchos::RCP<AmanziSolvers::LinearOperator<FlowMatrix, CompositeVector, CompositeVectorSpace> > 
      solver = factory.Create("gmres", plist, matrix_tmp, matrix_tmp);

  Y.PutScalar(0.0);
  int ierr = solver->ApplyInverse(X, Y);

  // if (ierr != 1) {
  //   std::cout << "Newton solver (" << solver->name() 
  //        << "): ||r||=" << solver->residual() << " itr=" << solver->num_itrs()
  //        << " code=" << solver->returned_code() << std::endl;
  //   exit(0);
  // }

  return 0;
}


/* ******************************************************************
x * 
****************************************************************** */
int Matrix_TPFA::ApplyInverse(const CompositeVector& X, CompositeVector& Y) const
{
  const Epetra_MultiVector& Xc = *X.ViewComponent("cell");
  Epetra_MultiVector& Yc = *Y.ViewComponent("cell");

  preconditioner_->ApplyInverse(Xc, Yc);
}


/* ******************************************************************
* 
****************************************************************** */
void Matrix_TPFA::ApplyBoundaryConditions(std::vector<int>& bc_model, 
                                          std::vector<bc_tuple>& bc_values)
{
  Epetra_MultiVector& Krel_faces = *rel_perm_->Krel().ViewComponent("face", true);

  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  AmanziMesh::Entity_ID_List faces;

  Epetra_MultiVector& rhs_cells = *rhs_->ViewComponent("cell");

  for (int c = 0; c < ncells; c++) {
    mesh_->cell_get_faces(c, &faces);
    int nfaces = faces.size();

    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];
      double value = bc_values[f][0];


      if (bc_model[f] == FLOW_BC_FACE_PRESSURE) {
	rhs_cells[0][c] += value * (*transmissibility_)[f] * Krel_faces[0][f];
      } else if (bc_model[f] == FLOW_BC_FACE_FLUX) {
        rhs_cells[0][c] -= value * mesh_->face_area(f);
	(*transmissibility_)[f] = 0.0;
	(*gravity_term_)[f] = 0.0;
      } else if (bc_model[f] == FLOW_BC_FACE_MIXED) {
	Errors::Message msg;
	msg << "Mixed boundary conditions are not supported in TPFA mode\n";
	Exceptions::amanzi_throw(msg);
      }
    }
  }
}


/* ******************************************************************
* Adds time derivative to the cell-based part of MFD algebraic system.
****************************************************************** */
void Matrix_TPFA::AddTimeDerivative(
    const Epetra_MultiVector& p, const Epetra_MultiVector& phi, double rho, double dT)
{
  Epetra_MultiVector dSdP(p);
  rel_perm_->DerivedSdP(p, dSdP);

  for (int c = 0; c < ncells_owned; c++) {
    double volume = mesh_->cell_volume(c);
    double factor = rho * phi[0][c] * dSdP[0][c] * volume / dT;
    Acc_cells_[c] += factor;
    Fc_cells_[c] += factor * p[0][c];
  }
}


/* ******************************************************************
* Linear algebra operations with matrices: r = A * u - f                                                 
****************************************************************** */
double Matrix_TPFA::ComputeNegativeResidual(const CompositeVector& u, CompositeVector& r)
{
  Epetra_MultiVector& Krel_faces = *rel_perm_->Krel().ViewComponent("face", true);

  AmanziMesh::Entity_ID_List faces;
  AmanziMesh::Entity_ID_List cells;

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
	if (c < ncells_owned){ 
	  rc[0][c] += sign * Krel_faces[0][f] * (*transmissibility_)[f] * (uc[0][cells[0]] - uc[0][cells[1]]);  
	}
        sign = -sign;
      }
    } else if (ncells == 1) {
      int c = cells[0];
      if (c < ncells_owned) { 
        rc[0][c] += Krel_faces[0][f] * (*transmissibility_)[f] * uc[0][c];
      }
    }							
  } 

  
  
  Epetra_MultiVector& rhs_cells = *rhs_->ViewComponent("cell");

  AmanziGeometry::Point face_centr, cell_cntr;


  for (int c = 0; c < ncells_owned; c++) {    
    rc[0][c] -= rhs_cells[0][c];   
    cell_cntr = mesh_->cell_centroid(c);
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
  Epetra_MultiVector& Krel_faces = *rel_perm_->Krel().ViewComponent("face", true);

  solution.ScatterMasterToGhosted("cell");
  const Epetra_MultiVector& p = *solution.ViewComponent("cell", true);
  Epetra_MultiVector& flux = *darcy_mass_flux.ViewComponent("face", true);

  AmanziMesh::Entity_ID_List faces;
  std::vector<double> dp;
  std::vector<int> dirs;

  AmanziMesh::Entity_ID_List cells;
  std::vector<int> flag(nfaces_wghost, 0);
  //face_flag_.assign(nfaces_wghost, 0);

  for (int c = 0; c < ncells_owned; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];

      if (bc_model[f] == FLOW_BC_FACE_PRESSURE) {

	double value = bc_values[f][0];
	flux[0][f] = dirs[n] * (*transmissibility_)[f] * (p[0][c] - value) + (*gravity_term_)[f];
	flux[0][f] *= Krel_faces[0][f];

      } else if (bc_model[f] == FLOW_BC_FACE_FLUX) {

	double value = bc_values[f][0];
	double area = mesh_->face_area(f);
	flux[0][f] = value*area;

      } else {
	if (f < nfaces_owned && !flag[f]) {
	  mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
	  if (cells.size() <= 1) {
	    Errors::Message msg("Flow PK: These boundary conditions are not supported by TPFA.");
	    Exceptions::amanzi_throw(msg);
	  }
          if (c == cells[0]){
            flux[0][f] = dirs[n] * (*transmissibility_)[f] * (p[0][cells[0]] - p[0][cells[1]]) + (*gravity_term_)[f];
          } else {
            flux[0][f] = dirs[n] * (*transmissibility_)[f] * (p[0][cells[1]] - p[0][cells[0]]) + (*gravity_term_)[f];
          }	    
          flux[0][f] *= Krel_faces[0][f];
          flag[f] = 1;
	}
      }
    }
  }
}


/* ******************************************************************
* Computation of the part of the jacobian which depends on
* analytical derivatives of relative permeabilities
****************************************************************** */
void Matrix_TPFA::AnalyticJacobian_(const CompositeVector& u, 
                                    std::vector<int>& bc_models,
                                    std::vector<bc_tuple>& bc_values)
{
  u.ScatterMasterToGhosted("cell");
  const Epetra_MultiVector& uc = *u.ViewComponent("cell", true);

  const Epetra_Map& cmap_wghost = mesh_->cell_map(true);
  AmanziMesh::Entity_ID_List cells, faces;
  std::vector<int> dirs;

  int cells_LID[2], cells_GID[2];
  double perm_abs_vert[2], perm_abs_horz[2];
  double k_rel[2], dk_dp[2], pres[2];
  AmanziGeometry::Point cntr_cell[2];
  double dist;

  const Epetra_MultiVector& Krel_cells = *rel_perm_->Krel().ViewComponent("cell");
  Epetra_MultiVector& Krel_faces = *rel_perm_->Krel().ViewComponent("face", true);

  const Epetra_MultiVector& dKdP_cells = *rel_perm_->dKdP().ViewComponent("cell");
  Epetra_MultiVector& dKdP_faces = *rel_perm_->dKdP().ViewComponent("face", true);
  int method = rel_perm_->method();

  //face_flag_.assign(nfaces_wghost, 0);
  std::vector<int> flag(nfaces_wghost, 0);
  //  for (int f = 0; f < nfaces_owned; f++){

  AmanziGeometry::Point face_centr;


  for (int c = 0; c < ncells_owned; c++){

    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();
    for (int i=0; i<nfaces; i++){

      int f = faces[i];
      if (f < nfaces_owned && !flag[f]) {
	mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
	int mcells = cells.size();
	Teuchos::SerialDenseMatrix<int, double> Jpp(mcells, mcells);

	

	int face_dir;   	
	// face_dir is equal to 1 if normal direction points from cells[0] ,
	// otherwise face_dir is -1
	if (cells[0] == c)  face_dir = dirs[i];
	else face_dir = -dirs[i];

	for (int n = 0; n < mcells; n++) {
	  cells_LID[n] = cells[n];
	  cells_GID[n] = cmap_wghost.GID(cells_LID[n]);
	  pres[n] = uc[0][cells_LID[n]];
	  dk_dp[n] = dKdP_cells[0][cells_LID[n]];
	}

	if (mcells == 1) {
	  dk_dp[0] = dKdP_faces[0][f];
	}

	const AmanziGeometry::Point& normal = mesh_->face_normal(f, false, cells[0]);
	ComputeJacobianLocal_(mcells, f, face_dir, method, bc_models, bc_values, pres, dk_dp, Jpp);


	(*Spp_).SumIntoGlobalValues(mcells, cells_GID, Jpp.values());

	flag[f] = 1;
      }

    }


  }
    
  Spp_->GlobalAssemble();

  //cout<<*Spp_<<endl;
}


/* ******************************************************************
* Computation of a local submatrix of 
* Analytical Jacobian (nonlinear part) on a particular face.
****************************************************************** */
void Matrix_TPFA::ComputeJacobianLocal_(
    int mcells, int face_id, int face_dir, int Krel_method,
    std::vector<int>& bc_models, std::vector<bc_tuple>& bc_values,
    double *pres, double *dk_dp_cell,
    Teuchos::SerialDenseMatrix<int, double>& Jpp)
{
  double dKrel_dp[2];
  double rho_w = *S_->GetScalarData("fluid_density");

  double dpres;
  if (mcells == 2) {
    dpres = pres[0] - pres[1];  // + grn;
    if (Krel_method == FLOW_RELATIVE_PERM_UPWIND_GRAVITY) {  // Define K and Krel_faces
      double cos_angle = face_dir * (*gravity_term_)[face_id] / (rho_w * mesh_->face_area(face_id));
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
      double flux0to1;
      flux0to1 = (*transmissibility_)[face_id] * dpres + face_dir * (*gravity_term_)[face_id];
      if (flux0to1  > FLOW_RELATIVE_PERM_TOLERANCE) {  // Upwind
        dKrel_dp[0] = dk_dp_cell[0];
        dKrel_dp[1] = 0.0;
      } else if ( flux0to1 < -FLOW_RELATIVE_PERM_TOLERANCE) {  // Upwind
        dKrel_dp[0] = 0.0;
        dKrel_dp[1] = dk_dp_cell[1];
      } else if (fabs(flux0to1) < FLOW_RELATIVE_PERM_TOLERANCE) {  // Upwind
        dKrel_dp[0] = 0.5*dk_dp_cell[0];
        dKrel_dp[1] = 0.5*dk_dp_cell[1];
      }
    } else if (Krel_method == FLOW_RELATIVE_PERM_ARITHMETIC_MEAN) {
      dKrel_dp[0] = 0.5*dk_dp_cell[0];
      dKrel_dp[1] = 0.5*dk_dp_cell[1];
    }

    Jpp(0, 0) = ((*transmissibility_)[face_id] * dpres + face_dir * (*gravity_term_)[face_id]) * dKrel_dp[0];
    Jpp(0, 1) = ((*transmissibility_)[face_id] * dpres + face_dir * (*gravity_term_)[face_id]) * dKrel_dp[1];
    Jpp(1, 0) = -Jpp(0, 0);
    Jpp(1, 1) = -Jpp(0, 1);

  } else if (mcells == 1) {
    if (bc_models[face_id] == FLOW_BC_FACE_PRESSURE) {                   
      pres[1] = bc_values[face_id][0];

      dpres = pres[0] - pres[1];  // + grn;
      Jpp(0,0) = ((*transmissibility_)[face_id] * dpres + face_dir * (*gravity_term_)[face_id]) * dk_dp_cell[0];
    } else {
      Jpp(0,0) = 0.0;
    }
  }
}


/* ******************************************************************
* Compute transmissibilities on faces 
****************************************************************** */
void Matrix_TPFA::ComputeTransmissibilities_()
{
  transmissibility_->PutScalar(0.0);

  double rho_ = *S_->GetScalarData("fluid_density");
  double mu_ = *S_->GetScalarData("fluid_viscosity");
  const Epetra_Vector& gravity_ = *S_->GetConstantVectorData("gravity");

  int dim = mesh_->space_dimension();
  AmanziGeometry::Point gravity(dim);
  for (int k = 0; k < dim; k++) gravity[k] = gravity_[k] * rho_;

  AmanziMesh::Entity_ID_List faces;
  AmanziMesh::Entity_ID_List cells;
  AmanziGeometry::Point a_dist, a[2];
  double h[2], perm[2], perm_test[2], h_test[2], beta[2];
  double trans_f;

  for (int f = 0; f < nfaces_owned; f++) {
    mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
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
      perm[i] = (rho_/mu_) * (((*K_)[cells[i]] * a[i]) * normal) * s;

      perm_test[i] = (rho_/mu_) * (((*K_)[cells[i]] * normal) * a_dist);
      //h_test[i] = pow(-1.0, i)*((face_centr - mesh_->cell_centroid(cells[i]))*normal) / area;  
      double dxn = a[i]*normal;
     
      beta[i] = fabs(perm[i] / dxn);
    }


    double  grav;
    grav = (gravity * normal) / area;

    trans_f = 0.0;

    if (ncells == 2){
      //factor = (perm[0]*perm[1]) / (h[0]*perm[1] + h[1]*perm[0]);
      grav *= (h[0] + h[1]);
      trans_f = (beta[0]*beta[1]) / (beta[0] + beta[1]);
    } else if (ncells == 1) {    
      //factor = perm[0] / h[0];
      grav *= h[0];
      trans_f = beta[0];
    } 



    (*transmissibility_)[f] = trans_f;
    (*gravity_term_)[f] = (*transmissibility_)[f] * grav;

  }

  // parallelization using CV capability
#ifdef HAVE_MPI
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
#endif
}

}  // namespace AmanziFlow
}  // namespace Amanzi


