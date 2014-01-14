/*
  This is the flow component of the Amanzi code.
  Frequently used bundles of routines are wrapped into computational blocks.

  Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include "Matrix_MFD.hh"
#include "Matrix_TPFA.hh"
#include "Richards_PK.hh"


namespace Amanzi {
namespace AmanziFlow {

/* ******************************************************************
*                       LEVEL 1 subroutines                         
****************************************************************** */

/* ******************************************************************
* A wrapper for updating boundary conditions.
****************************************************************** */
void Richards_PK::UpdateSourceBoundaryData(double Tp, const CompositeVector& u)
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

  ComputeBCs(u);
}


/* ******************************************************************
* A wrapper for generating a steady state matrix. 
* Warning: Krel must be initialized before calling this routine. 
****************************************************************** */
void Richards_PK::AssembleSteadyStateMatrix(FlowMatrix* matrix)
{ 
  matrix->CreateStiffnessMatricesRichards();
  matrix->CreateRHSVectors();
  matrix->AddGravityFluxesRichards(rho_, gravity_, bc_model);
  matrix->ApplyBoundaryConditions(bc_model, bc_values);
  matrix->Assemble();
  /*
    Epetra_Vector& Krel_faces = rel_perm->Krel_faces();

    matrix ->ApplyBoundaryConditions(bc_model, bc_values);
    matrix_->AddGravityFluxes(Krel_faces, *Grav_term_faces, bc_model, &*matrix);
    matrix ->Assemble();
  */
}


/* ******************************************************************
* A wrapper for generating a steady state preconditioner. 
* Warning: Krel must be initialized before calling this routine. 
****************************************************************** */
void Richards_PK::AssembleSteadyStatePreconditioner(FlowMatrix* preconditioner)
{ 
  CompositeVector& u = *preconditioner->rhs();  // TODO u is dummy

  preconditioner->CreateStiffnessMatricesRichards();
  preconditioner->CreateRHSVectors();
  preconditioner->ApplyBoundaryConditions(bc_model, bc_values);
  preconditioner->AssembleDerivatives(u, bc_model, bc_values);
}


/* ******************************************************************
*                       LEVEL 2 subroutines                         
****************************************************************** */

/* ******************************************************************
* Gathers together routines to compute MFD matrices.                            
****************************************************************** */
void Richards_PK::AssembleMatrixMFD(const CompositeVector& u, double Tp)
{
  rel_perm->Compute(u, bc_model, bc_values);
  UpdateSourceBoundaryData(Tp, u);
  
  matrix_->CreateStiffnessMatricesRichards();
  matrix_->CreateRHSVectors();
  matrix_->AddGravityFluxesRichards(rho_, gravity_, bc_model);
  matrix_->ApplyBoundaryConditions(bc_model, bc_values);
  matrix_->Assemble();

  Teuchos::RCP<CompositeVector> rhs = matrix_->rhs();
  if (src_sink != NULL) AddSourceTerms(*rhs);
}


/* ******************************************************************
* Gathers together routines to compute MFD matrices Axx(u) and
* preconditioner Sff(u) using time step dT.                             
****************************************************************** */
void Richards_PK::AssemblePreconditionerMFD(const CompositeVector& u, double Tp, double dTp)
{
  // update all coefficients, boundary data, and source/sink terms
  const Epetra_MultiVector& u_cells = *u.ViewComponent("cell");

  rel_perm->Compute(u, bc_model, bc_values);
  UpdateSourceBoundaryData(Tp, u);

  preconditioner_->CreateStiffnessMatricesRichards();
  preconditioner_->CreateRHSVectors();

  if (dTp > 0.0) {
    const Epetra_MultiVector& phi = *S_->GetFieldData("porosity")->ViewComponent("cell");
    preconditioner_->AddTimeDerivative(u_cells, phi, rho_, dTp);
  }
  
  preconditioner_->ApplyBoundaryConditions(bc_model, bc_values);
  preconditioner_->AssembleDerivatives(u, bc_model, bc_values);
  preconditioner_->UpdatePreconditioner();
}

}  // namespace AmanziFlow
}  // namespace Amanzi



