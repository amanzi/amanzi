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
#include "Matrix_MFD_TPFA.hh"
#include "Matrix_MFD_PLambda.hh"
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
    if (src_sink_distribution & Amanzi::DOMAIN_FUNCTION_ACTION_DISTRIBUTE_PERMEABILITY)
      src_sink->ComputeDistribute(Tp, Kxy->Values()); 
    else
      src_sink->ComputeDistribute(Tp, NULL);
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
  matrix->CreateMFDstiffnessMatrices(*rel_perm);
  matrix->CreateMFDrhsVectors();
  AddGravityFluxes_MFD(K, matrix, *rel_perm);
  matrix->ApplyBoundaryConditions(bc_model, bc_values);
  matrix->AssembleGlobalMatrices();
}


/* ******************************************************************
* A wrapper for generating a steady state preconditioner. 
* Warning: Krel must be initialized before calling this routine. 
****************************************************************** */
void Richards_PK::AssembleSteadyStatePreconditioner_MFD(Matrix_MFD* preconditioner)
{ 
  preconditioner->CreateMFDstiffnessMatrices(*rel_perm);
  preconditioner->CreateMFDrhsVectors();
  preconditioner->ApplyBoundaryConditions(bc_model, bc_values);
  preconditioner->AssembleSchurComplement(bc_model, bc_values);
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

  // setup a new algebraic problem
  matrix_->CreateMFDstiffnessMatrices(*rel_perm);
  matrix_->CreateMFDrhsVectors();
  AddGravityFluxes_MFD(K, matrix_, *rel_perm);
  matrix_->ApplyBoundaryConditions(bc_model, bc_values);
  matrix_->AssembleGlobalMatrices();
 
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

  // update all coefficients, boundary data, and source/sink terms
  rel_perm->Compute(u, bc_model, bc_values);
  UpdateSourceBoundaryData(Tp, *u_cells, *u_faces);

  // setup a new algebraic problem
  preconditioner_->CreateMFDstiffnessMatrices(*rel_perm);
  preconditioner_->CreateMFDrhsVectors();
  AddTimeDerivative_MFD(*u_cells, dTp, preconditioner_);

  if (experimental_solver_ == FLOW_SOLVER_PICARD_NEWTON) {
    Matrix_MFD_PLambda* matrix_plambda = static_cast<Matrix_MFD_PLambda*>(preconditioner_);
    Epetra_Vector& flux = FS->ref_darcy_flux();
    Epetra_Vector& Krel_faces = rel_perm->Krel_faces();
    rhs = preconditioner_->rhs();
    AddNewtonFluxes_MFD(*dKdP_faces, Krel_faces, *u_cells, flux, *rhs, matrix_plambda);
  }

  preconditioner_->ApplyBoundaryConditions(bc_model, bc_values);
  preconditioner_->AssembleSchurComplement(bc_model, bc_values);
 
  if (experimental_solver_ == FLOW_SOLVER_NEWTON) {
    Matrix_MFD_TPFA* matrix_tpfa = dynamic_cast<Matrix_MFD_TPFA*>(preconditioner_);
    if (matrix_tpfa == 0) {
      Errors::Message msg;
      msg << "Flow PK: cannot cast pointer to class Matrix_MFD_TPFA\n";
      Exceptions::amanzi_throw(msg);
    }

    Epetra_Vector& Krel_cells = rel_perm->Krel_cells();
    Epetra_Vector& Krel_faces = rel_perm->Krel_faces();
    matrix_tpfa->AnalyticJacobian(*u_cells, dim, Krel_method, bc_model, bc_values,
                                  Krel_cells, *dKdP_cells,
                                  Krel_faces, *dKdP_faces);
  }

  preconditioner_->UpdatePreconditioner();
}

}  // namespace AmanziFlow
}  // namespace Amanzi



