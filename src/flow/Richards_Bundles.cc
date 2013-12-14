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
void Richards_PK::UpdateSourceBoundaryData(double Tp, const CompositeVector& pressure)
{
  const Epetra_MultiVector& p_cells = *pressure.ViewComponent("cell");
  const Epetra_MultiVector& p_faces = *pressure.ViewComponent("face"); 

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

  ComputeBCs(pressure);
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
  preconditioner->CreateStiffnessMatricesRichards();
  preconditioner->CreateRHSVectors();
  preconditioner->ApplyBoundaryConditions(bc_model, bc_values);
  preconditioner->AssembleSchurComplement(bc_model, bc_values);

  /*
    Epetra_Vector& Krel_faces = rel_perm->Krel_faces();

    preconditioner -> ApplyBoundaryConditions(bc_model, bc_values);
    AddGravityFluxes_TPFA( Krel_faces, *Grav_term_faces, bc_model, preconditioner);
    preconditioner -> AssembleGlobalMatrices();  
  */ 
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
  
  /* TODO
  if (experimental_solver_ == FLOW_SOLVER_NEWTON) {
    Teuchos::RCP<Epetra_Vector> rhs_cells_ = matrix_->rhs_cells();
    rhs_cells_->PutScalar(0.0);

    Epetra_Vector& Krel_faces = rel_perm->Krel_faces();

    matrix_->ApplyBoundaryConditions(bc_model, bc_values);
    AddGravityFluxes_TPFA(Krel_faces, *Grav_term_faces, bc_model, &*matrix_);
    matrix_->AssembleGlobalMatrices();
  } else{
  */
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
  /* TODO
  if (experimental_solver_ == FLOW_SOLVER_NEWTON) {
    matrix_->DeriveMassFlux(*u_cells, flux, bc_model, bc_values);
    for (int f = 0; f < nfaces_owned; f++) flux[f] /= rho_;    
  }
  */

  const Epetra_MultiVector& u_cells = *u.ViewComponent("cell");

  rel_perm->Compute(u, bc_model, bc_values);
  UpdateSourceBoundaryData(Tp, u);

  preconditioner_->CreateStiffnessMatricesRichards();
  preconditioner_->CreateRHSVectors();
  /* TODO
    std::vector<double>& Acc_cells = preconditioner_->Acc_cells();

    double* Ac = Acc_cells.data();
    int nsize = Acc_cells.size();
    memset(Ac, 0.0, nsize*sizeof(double));
  */

  if (dTp > 0.0) {
    const Epetra_MultiVector& phi = *S_->GetFieldData("porosity")->ViewComponent("cell");
    preconditioner_->AddTimeDerivative(u_cells, phi, rho_, dTp);
  }
  
  preconditioner_->ApplyBoundaryConditions(bc_model, bc_values);
  preconditioner_->AssembleSchurComplement(bc_model, bc_values);
  /* TODO
    preconditioner_->ApplyBoundaryConditions(bc_model, bc_values);
    preconditioner_->AssembleSchurComplement(bc_model, bc_values);

    preconditioner_->AnalyticJacobian(*u_cells, bc_model, bc_values, *rel_perm);
  */

  preconditioner_->UpdatePreconditioner();
}

}  // namespace AmanziFlow
}  // namespace Amanzi



