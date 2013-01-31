/*
This is the flow component of the Amanzi code.
Frequently used bundles of routines are wrapped into computational blocks. 

Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided Reconstruction.cppin the top-level COPYRIGHT file.

Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include "Flow_State.hpp"
#include "Matrix_MFD.hpp"
#include "Matrix_MFD_TPFA.hpp"
#include "Matrix_MFD_PLambda.hpp"
#include "Richards_PK.hpp"

namespace Amanzi {
namespace AmanziFlow {

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

  // use the exisiting mass matrices to recover the Darcy flux
  // may exist problems with the first call (lipnikov@lanl.gov)
  /*
  if (experimental_solver_ == FLOW_SOLVER_PICARD_NEWTON) {
    Epetra_Vector& flux = FS->ref_darcy_flux();
    preconditioner_->DeriveDarcyMassFlux(u, *face_importer_, flux);
    AddGravityFluxes_DarcyFlux(K, *Krel_cells, *Krel_faces, Krel_method, flux);
    for (int f = 0; f < nfaces_owned; f++) flux[f] /= rho;
  }
  */

  // update all coefficients, boundary data, and source/sink terms
  CalculateRelativePermeability(u);
  UpdateSourceBoundaryData(Tp, *u_faces);

  // setup a new algebraic problem
  preconditioner_->CreateMFDstiffnessMatrices(*Krel_cells, *Krel_faces, Krel_method);
  preconditioner_->CreateMFDrhsVectors();
  AddGravityFluxes_MFD(K, *Krel_cells, *Krel_faces, Krel_method, preconditioner_);
  AddTimeDerivative_MFD(*u_cells, dTp, preconditioner_);
  if (experimental_solver_ == FLOW_SOLVER_PICARD_NEWTON) {
    Matrix_MFD_PLambda* matrix_plambda = static_cast<Matrix_MFD_PLambda*>(preconditioner_);
    Epetra_Vector& flux = FS->ref_darcy_flux();
    rhs = preconditioner_->rhs();
    AddNewtonFluxes_MFD(*dKdP_faces, *Krel_faces, *u_cells, flux, *rhs, matrix_plambda);
  }
  preconditioner_->ApplyBoundaryConditions(bc_model, bc_values);
  preconditioner_->AssembleGlobalMatrices();
 
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

  preconditioner_->ComputeSchurComplement(bc_model, bc_values);
  preconditioner_->UpdatePreconditioner();
}


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

  // DEBUG
  // if (MyPID == 0 && verbosity >= FLOW_VERBOSITY_EXTREME)
  //     std::printf("Flow PK: updating source/boundary data at T=%14.9e [sec]\n", Tp);
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
* A wrapper for generating a steady state problem. 
* Warning: Krel must be initialized before calling this routine. 
* OBSOLETE routine: It will eventually phase away.
****************************************************************** */
void Richards_PK::AssembleSteadyStateProblem_MFD(Matrix_MFD* matrix_operator, bool add_preconditioner)
{ 
  matrix_operator->CreateMFDstiffnessMatrices(*Krel_cells, *Krel_faces, Krel_method);
  matrix_operator->CreateMFDrhsVectors();
  AddGravityFluxes_MFD(K, *Krel_cells, *Krel_faces, Krel_method, matrix_operator);
  matrix_operator->ApplyBoundaryConditions(bc_model, bc_values);
  matrix_operator->AssembleGlobalMatrices();

  if (add_preconditioner) {
    matrix_operator->ComputeSchurComplement(bc_model, bc_values);
  }

  // DEBUG
  // Matrix_Audit audit(mesh_, matrix);
  // audit.InitAudit();
  // audit.CheckSpectralBounds();
}


/* ******************************************************************
* A wrapper for generating a transient problem. 
* Warning: Krel must be initialized before calling this routine. 
* OBSOLETE routine: It will eventually phase away.
****************************************************************** */
void Richards_PK::AssembleTransientProblem_MFD(Matrix_MFD* matrix_operator, double dTp,
                                               Epetra_Vector& p, bool add_preconditioner)
{ 
  matrix_operator->CreateMFDstiffnessMatrices(*Krel_cells, *Krel_faces, Krel_method);
  matrix_operator->CreateMFDrhsVectors();
  AddGravityFluxes_MFD(K, *Krel_cells, *Krel_faces, Krel_method, matrix_operator);
  AddTimeDerivative_MFD(p, dTp, matrix_operator);
  matrix_operator->ApplyBoundaryConditions(bc_model, bc_values);
  matrix_operator->AssembleGlobalMatrices();

  if (add_preconditioner) {
    matrix_operator->ComputeSchurComplement(bc_model, bc_values);
  }
}

}  // namespace AmanziFlow
}  // namespace Amanzi



