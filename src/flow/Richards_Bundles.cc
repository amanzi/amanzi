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
* A wrapper for updating relative permeabilities.
****************************************************************** */
void Richards_PK::CalculateRelativePermeability(const Epetra_Vector& u)
{
  Epetra_Vector* u_cells = FS->CreateCellView(u);
  Epetra_Vector* u_faces = FS->CreateFaceView(u);

  if (Krel_method == FLOW_RELATIVE_PERM_UPWIND_GRAVITY ||
      Krel_method == FLOW_RELATIVE_PERM_UPWIND_DARCY_FLUX ||
      Krel_method == FLOW_RELATIVE_PERM_ARITHMETIC_MEAN) {
    CalculateRelativePermeabilityFace(*u_cells);
    Krel_cells->PutScalar(1.0);
    if (experimental_solver_ == FLOW_SOLVER_NEWTON || 
        experimental_solver_ == FLOW_SOLVER_PICARD_NEWTON) {
      CalculateDerivativePermeabilityFace(*u_cells);
    }
  } else if (Krel_method == FLOW_RELATIVE_PERM_EXPERIMENTAL) {
    CalculateRelativePermeabilityFace(*u_cells);
  } else {
    CalculateRelativePermeabilityCell(*u_cells);
    Krel_faces->PutScalar(1.0);
  }
}

 
/* ******************************************************************
* A wrapper for updating boundary conditions.
* Warning: routine is marked as obsolete.
****************************************************************** */
void Richards_PK::UpdateBoundaryConditions(double Tp, Epetra_Vector& p_faces)
{
  bc_pressure->Compute(Tp);
  bc_flux->Compute(Tp);
  bc_seepage->Compute(Tp);
  if (shift_water_table_.getRawPtr() == NULL)
    bc_head->Compute(Tp);
  else
    bc_head->ComputeShift(Tp, shift_water_table_->Values());

  ProcessBoundaryConditions(
      bc_pressure, bc_head, bc_flux, bc_seepage,
      p_faces, atm_pressure, rainfall_factor,
      bc_submodel, bc_model, bc_values);
}


/* ******************************************************************
* A wrapper for updating boundary conditions.
****************************************************************** */
void Richards_PK::UpdateSourceBoundaryData(double Tp, Epetra_Vector& p_faces)
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
      p_faces, atm_pressure, rainfall_factor,
      bc_submodel, bc_model, bc_values);
}


/* ******************************************************************
* A wrapper for generating a steady state matrix. 
* Warning: Krel must be initialized before calling this routine. 
****************************************************************** */
void Richards_PK::AssembleSteadyStateMatrix_MFD(Matrix_MFD* matrix)
{ 
  matrix->CreateMFDstiffnessMatrices(*Krel_cells, *Krel_faces, Krel_method);
  matrix->CreateMFDrhsVectors();
  AddGravityFluxes_MFD(K, *Krel_cells, *Krel_faces, Krel_method, matrix);
  matrix->ApplyBoundaryConditions(bc_model, bc_values);
  matrix->AssembleGlobalMatrices();
}


/* ******************************************************************
* A wrapper for generating a steady state preconditioner. 
* Warning: Krel must be initialized before calling this routine. 
****************************************************************** */
void Richards_PK::AssembleSteadyStatePreconditioner_MFD(Matrix_MFD* preconditioner)
{ 
  preconditioner->CreateMFDstiffnessMatrices(*Krel_cells, *Krel_faces, Krel_method);
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

  CalculateRelativePermeability(u);
  UpdateSourceBoundaryData(Tp, *u_faces);

  // setup a new algebraic problem
  matrix_->CreateMFDstiffnessMatrices(*Krel_cells, *Krel_faces, Krel_method);
  matrix_->CreateMFDrhsVectors();
  AddGravityFluxes_MFD(K, *Krel_cells, *Krel_faces, Krel_method, matrix_);
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
  CalculateRelativePermeability(u);
  UpdateSourceBoundaryData(Tp, *u_faces);

  // setup a new algebraic problem
  preconditioner_->CreateMFDstiffnessMatrices(*Krel_cells, *Krel_faces, Krel_method);
  preconditioner_->CreateMFDrhsVectors();
  AddTimeDerivative_MFD(*u_cells, dTp, preconditioner_);

  if (experimental_solver_ == FLOW_SOLVER_PICARD_NEWTON) {
    Matrix_MFD_PLambda* matrix_plambda = static_cast<Matrix_MFD_PLambda*>(preconditioner_);
    Epetra_Vector& flux = FS->ref_darcy_flux();
    rhs = preconditioner_->rhs();
    AddNewtonFluxes_MFD(*dKdP_faces, *Krel_faces, *u_cells, flux, *rhs, matrix_plambda);
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

    matrix_tpfa->AnalyticJacobian(*u_cells, dim, Krel_method, bc_model, bc_values,
                                  *Krel_cells, *dKdP_cells,
                                  *Krel_faces, *dKdP_faces);
  }

  preconditioner_->UpdatePreconditioner();
}

}  // namespace AmanziFlow
}  // namespace Amanzi



