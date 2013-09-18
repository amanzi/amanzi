/*
This is the flow component of the Amanzi code.
Frequently used bundles of routines are wrapped into computational blocks.

Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided Reconstruction.cppin the top-level COPYRIGHT file.

Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include "Flow_State.hh"
#include "Matrix_MFD.hh"
#include "Matrix_TPFA.hh"
// #include "Matrix_MFD_PLambda.hh"
#include "Richards_PK.hh"


namespace Amanzi {
namespace AmanziFlow {

/* ******************************************************************
*                       LEVEL 1 subroutines                         
****************************************************************** */

/* ******************************************************************
* A wrapper for updating boundary conditions.
****************************************************************** */
void Richards_PK::UpdateSourceBoundaryData(
    double Tp, Epetra_Vector& pressure, Epetra_Vector& lambda)
{
  if (src_sink != NULL) {
    if (src_sink_distribution & Amanzi::Functions::DOMAIN_FUNCTION_ACTION_DISTRIBUTE_PERMEABILITY) {
      src_sink->ComputeDistribute(Tp, Kxy->Values()); 
    } else {
      src_sink->ComputeDistribute(Tp, NULL);
    }
  }

  bc_pressure->Compute(Tp);
  bc_flux->Compute(Tp);
  bc_seepage->Compute(Tp);
  if (shift_water_table_.getRawPtr() == NULL)
    bc_head->Compute(Tp);
  else
    bc_head->ComputeShift(Tp, shift_water_table_->Values());

  ProcessBoundaryConditions(
      bc_pressure, bc_head, bc_flux, bc_seepage,
      pressure, lambda, atm_pressure, rainfall_factor,
      bc_submodel, bc_model, bc_values);
}


/* ******************************************************************
* A wrapper for generating a steady state matrix. 
* Warning: Krel must be initialized before calling this routine. 
****************************************************************** */
void Richards_PK::AssembleSteadyStateMatrix_MFD(Matrix_MFD* matrix)
{ 
  if (experimental_solver_ != FLOW_SOLVER_NEWTON) {
    matrix->CreateMFDstiffnessMatrices(*rel_perm);
    matrix->CreateMFDrhsVectors();
    AddGravityFluxes_MFD(K, matrix, *rel_perm);
    matrix->ApplyBoundaryConditions(bc_model, bc_values);
    matrix->AssembleGlobalMatrices();
  }
  else {
    Matrix_MFD_TPFA* matrix_tpfa = dynamic_cast<Matrix_MFD_TPFA*>(&*matrix_);
    if (matrix_tpfa == 0) {
      Errors::Message msg;
      msg << "Richards PK: cannot cast pointer to class Matrix_MFD_TPFA\n";
      Exceptions::amanzi_throw(msg);
    }
    Epetra_Vector& Krel_faces = rel_perm->Krel_faces();
    //ComputeTransmissibilities(*Transmis_faces, *Grav_term_faces);
    matrix_tpfa -> ApplyBoundaryConditions(bc_model, bc_values, Krel_faces, *Transmis_faces, *Grav_term_faces);
    AddGravityFluxes_TPFA( Krel_faces, *Grav_term_faces, bc_model, matrix_tpfa);
    matrix_tpfa -> AssembleGlobalMatrices(Krel_faces, *Transmis_faces);
  }
}


/* ******************************************************************
* A wrapper for generating a steady state preconditioner. 
* Warning: Krel must be initialized before calling this routine. 
****************************************************************** */
void Richards_PK::AssembleSteadyStatePreconditioner_MFD(Matrix_MFD* preconditioner)
{ 
  if (experimental_solver_ != FLOW_SOLVER_NEWTON) {
    preconditioner->CreateMFDstiffnessMatrices(*rel_perm);
    preconditioner->CreateMFDrhsVectors();
    preconditioner->ApplyBoundaryConditions(bc_model, bc_values);
    preconditioner->AssembleSchurComplement(bc_model, bc_values);
  }
  else {
    Matrix_MFD_TPFA* prec_tpfa = dynamic_cast<Matrix_MFD_TPFA*>(preconditioner);
    if (prec_tpfa == 0) {
      Errors::Message msg;
      msg << "Richards PK: cannot cast pointer to class Matrix_MFD_TPFA\n";
      Exceptions::amanzi_throw(msg);
    }
    Epetra_Vector& Krel_faces = rel_perm->Krel_faces();
    //ComputeTransmissibilities(*Transmis_faces, *Grav_term_faces);
    prec_tpfa -> ApplyBoundaryConditions(bc_model, bc_values, Krel_faces, *Transmis_faces, *Grav_term_faces);
    AddGravityFluxes_TPFA( Krel_faces, *Grav_term_faces, bc_model, prec_tpfa);
    prec_tpfa -> AssembleGlobalMatrices(Krel_faces, *Transmis_faces);  
  }  
}


/* ******************************************************************
*                       LEVEL 2 subroutines                         
****************************************************************** */

/* ******************************************************************
* Gathers together routines to compute MFD matrices.                            
****************************************************************** */
void Richards_PK::AssembleMatrixMFD(const Epetra_Vector& u, double Tp)
{
  Epetra_Vector* u_cells = FS->CreateCellView(u);
  Epetra_Vector* u_faces = FS->CreateFaceView(u);

  rel_perm->Compute(u, bc_model, bc_values);
  UpdateSourceBoundaryData(Tp, *u_cells, *u_faces);
  
  if (experimental_solver_ == FLOW_SOLVER_NEWTON) {
    Matrix_MFD_TPFA* matrix_tpfa = dynamic_cast<Matrix_MFD_TPFA*>(&*matrix_);
    if (matrix_tpfa == 0) {
      Errors::Message msg;
      msg << "Flow PK: cannot cast pointer to class Matrix_MFD_TPFA\n";
      Exceptions::amanzi_throw(msg);
    }

    Teuchos::RCP<Epetra_Vector>& rhs_cells_ = matrix_tpfa -> rhs_cells();
    rhs_cells_->PutScalar(0.0);

    //Epetra_Vector& Trans_faces = matrix_tpfa -> ref_trans_faces();
    //Epetra_Vector& grav_faces = matrix_tpfa -> ref_grav_term_faces();
    Epetra_Vector& Krel_faces = rel_perm->Krel_faces();
    //ComputeTransmissibilities(*Transmis_faces, *Grav_term_faces);
    matrix_tpfa->ApplyBoundaryConditions(
        bc_model, bc_values, Krel_faces, *Transmis_faces, *Grav_term_faces);
    AddGravityFluxes_TPFA( Krel_faces, *Grav_term_faces, bc_model, matrix_tpfa);
    matrix_tpfa->AssembleGlobalMatrices(Krel_faces, *Transmis_faces);
  } else{
    // setup a new algebraic problem
    matrix_->CreateMFDstiffnessMatrices(*rel_perm);
    matrix_->CreateMFDrhsVectors();
    AddGravityFluxes_MFD(K, &*matrix_, *rel_perm);
    matrix_->ApplyBoundaryConditions(bc_model, bc_values);
    matrix_->AssembleGlobalMatrices();
  }

  rhs = matrix_->rhs();
  if (src_sink != NULL) AddSourceTerms(src_sink, *rhs);
}


/* ******************************************************************
* Gathers together routines to compute MFD matrices Axx(u) and
* preconditioner Sff(u) using time step dT.                             
****************************************************************** */
void Richards_PK::AssemblePreconditionerMFD(const Epetra_Vector& u, double Tp, double dTp)
{
  Epetra_Vector* u_cells = FS->CreateCellView(u);
  Epetra_Vector* u_faces = FS->CreateFaceView(u);
  Epetra_Vector& Krel_faces = rel_perm->Krel_faces();
  Epetra_Vector& flux = FS_aux->ref_darcy_flux();

  // update all coefficients, boundary data, and source/sink terms
  //Amanzi::timer_manager.stop("Update precon");
  if (experimental_solver_ == FLOW_SOLVER_NEWTON) {
    Matrix_MFD_TPFA* matrix_tpfa = dynamic_cast<Matrix_MFD_TPFA*>(&*matrix_);
    if (matrix_tpfa == 0) {
      Errors::Message msg;
      msg << "Richards_PK: cannot cast pointer to class Matrix_MFD_TPFA\n";
      Exceptions::amanzi_throw(msg);
    }

    matrix_tpfa->DeriveDarcyMassFlux(*u_cells, Krel_faces, *Transmis_faces, *Grav_term_faces, bc_model, bc_values, flux);

    for (int f = 0; f < nfaces_owned; f++) flux[f] /= rho_;    

    // cout<<"After DeriveDarcyMassFlux\n"<<endl;
    // for (int f = 48; f < nfaces_owned; f++){
    //   cout<<"flux "<<f<<":   "<<flux[f]<<endl;
    // }
  }

  rel_perm->Compute(u, bc_model, bc_values);
  UpdateSourceBoundaryData(Tp, *u_cells, *u_faces);
  //exit(0);
  //Amanzi::timer_manager.start("Update precon");

  // setup a new algebraic problem
  if (experimental_solver_ != FLOW_SOLVER_NEWTON) {
    preconditioner_->CreateMFDstiffnessMatrices(*rel_perm);
    preconditioner_->CreateMFDrhsVectors();
  } else {
    //u_faces -> PutScalar(0.0);
    //ComputeTransmissibilities(*Transmis_faces, *Grav_term_faces);
    std::vector<double>& Acc_cells = preconditioner_->Acc_cells();
    //std::vector<double>& Fc_cells = preconditioner_->Fc_cells();

    double* Ac = Acc_cells.data();
    //double* Fc = Fc_cells.data();
    int nsize = Acc_cells.size();
    memset(Ac, 0., nsize*sizeof(double));
  }

  AddTimeDerivative_MFD(*u_cells, dTp, &*preconditioner_);
  
  // if (experimental_solver_ == FLOW_SOLVER_PICARD_NEWTON) {
  //   Matrix_MFD_PLambda* matrix_plambda = static_cast<Matrix_MFD_PLambda*>(preconditioner_);
  //   Epetra_Vector& flux = FS->ref_darcy_flux();
  //   Epetra_Vector& Krel_faces = rel_perm->Krel_faces();
  //   rhs = preconditioner_->rhs();
  //   AddNewtonFluxes_MFD(*rel_perm, *u_cells, flux, *rhs, matrix_plambda);
  // }

  if (experimental_solver_ != FLOW_SOLVER_NEWTON) {
    preconditioner_->ApplyBoundaryConditions(bc_model, bc_values);
    preconditioner_->AssembleSchurComplement(bc_model, bc_values);
  } else {
    Matrix_MFD_TPFA* matrix_tpfa = dynamic_cast<Matrix_MFD_TPFA*>(&*preconditioner_);
    if (matrix_tpfa == 0) {
      Errors::Message msg;
      msg << "Flow PK: cannot cast pointer to class Matrix_MFD_TPFA\n";
      Exceptions::amanzi_throw(msg);
    }

    //Amanzi::timer_manager.start("AnalyticJacobian");
    Epetra_Vector& Krel_faces = rel_perm->Krel_faces();

    matrix_tpfa->ApplyBoundaryConditions(bc_model, bc_values, Krel_faces, *Transmis_faces, *Grav_term_faces);
    matrix_tpfa->AssembleSchurComplement(Krel_faces, *Transmis_faces);

    matrix_tpfa->AnalyticJacobian(*u_cells, dim, bc_model, bc_values, *Transmis_faces, *Grav_term_faces, *rel_perm);
    //Amanzi::timer_manager.stop("AnalyticJacobian");
  }

  preconditioner_->UpdatePreconditioner();
}


/* ******************************************************************
* Compute transmissibilities on faces 
****************************************************************** */
void Richards_PK::ComputeTransmissibilities(Epetra_Vector& Trans_faces, Epetra_Vector& grav_faces)
{
  Trans_faces.PutScalar(0.0);

  double rho = FS->ref_fluid_density();
  double vis = FS->ref_fluid_viscosity();

  AmanziGeometry::Point gravity(dim);
  for (int k = 0; k < dim; k++) gravity[k] = (*(FS->gravity()))[k] * rho;

  AmanziMesh::Entity_ID_List faces;
  AmanziMesh::Entity_ID_List cells;
  AmanziGeometry::Point a_dist;
  double h[2], perm[2], perm_test[2], h_test[2];
  double trans_f;

  Epetra_Vector& Krel_faces = rel_perm->Krel_faces();

  std::vector<int> dirs;

  for (int f = 0; f < nfaces_owned; f++) {
    mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
    int ncells = cells.size();
    const AmanziGeometry::Point& normal = mesh_->face_normal(f, false, cells[0]);
    const AmanziGeometry::Point& face_centr = mesh_->face_centroid(f);
    double area = mesh_->face_area(f);

    if (ncells == 2){
      a_dist = mesh_->cell_centroid(cells[1]) - mesh_->cell_centroid(cells[0]);
    } else if (ncells == 1) {    
      a_dist = face_centr - mesh_->cell_centroid(cells[0]);
    } 

    a_dist *= 1./norm(a_dist);

    for (int i=0; i<ncells; i++) {
      h[i] = norm(face_centr - mesh_->cell_centroid(cells[i]));
      perm[i] = (rho/vis) * ((K[cells[i]] * normal) * normal) / area;

      perm_test[i] = (rho/vis) * ((K[cells[i]] * normal) * a_dist);
      h_test[i] = pow(-1.0, i)*((face_centr - mesh_->cell_centroid(cells[i]))*normal) / area;
      //perm[i] = ((K[cells[i]] * normal) * normal)/ area;
      //cout<<i<<" "<<h[i]<<" "<<h_test[i]<<endl;
      //cout<<i<<" "<<perm[i]<<" "<<perm_test[i]<<endl;
    }
    double factor, grav;

    grav = (gravity * normal) / area;

    if (ncells == 2){
      factor = (perm[0]*perm[1]) / (h[0]*perm[1] + h[1]*perm[0]);
      grav *= (h[0] + h[1]);
    } else if (ncells == 1) {    
      factor = perm[0]/h[0];
      grav *= h[0];
    } 

    trans_f = 0.;
    for (int i=0; i<ncells; i++) {
      trans_f += h_test[i]/perm[i];
    }

    trans_f = 1./trans_f;

    //cout<<factor<<" "<<trans_f<<endl;

    Trans_faces[f] = trans_f;
    grav_faces[f] = Trans_faces[f] * grav;
  }


  FS->CopyMasterFace2GhostFace(Trans_faces);
  FS->CopyMasterFace2GhostFace(grav_faces);
}


}  // namespace AmanziFlow
}  // namespace Amanzi



