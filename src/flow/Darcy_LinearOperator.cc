/*
This is the flow component of the Amanzi code.

Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided in the top-level COPYRIGHT file.

Authors: Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
*/

#include <vector>

#include "Darcy_PK.hh"
#include "LinearOperatorFactory.hh"

namespace Amanzi {
namespace AmanziFlow {

/* ******************************************************************
* Gathers together routines to compute steady-state MFD matrices.                            
****************************************************************** */
void Darcy_PK::AssembleMatrixMFD()
{
  matrix_->CreateMFDstiffnessMatrices();
  matrix_->CreateMFDrhsVectors();
  AddGravityFluxes_MFD(&*matrix_);
  matrix_->ApplyBoundaryConditions(bc_model, bc_values);
  matrix_->AssembleGlobalMatrices();
  matrix_->AssembleSchurComplement(bc_model, bc_values);
  matrix_->UpdatePreconditioner();
}


/* ******************************************************************
* Calculates steady-state solution assuming that absolute permeability 
* does not depend on time. The boundary conditions are calculated
* only once, during the initialization step.                                                
****************************************************************** */
void Darcy_PK::SolveFullySaturatedProblem(double Tp, CompositeVector& u)
{
  // calculate and assemble elemental stifness matrices
  matrix_->CreateMFDstiffnessMatrices();
  matrix_->CreateMFDrhsVectors();
  AddGravityFluxes_MFD(&*matrix_);
  matrix_->ApplyBoundaryConditions(bc_model, bc_values);
  matrix_->AssembleGlobalMatrices();
  matrix_->AssembleSchurComplement(bc_model, bc_values);
  matrix_->UpdatePreconditioner();

  // create linear solver
  LinearSolver_Specs& ls_specs = ti_specs->ls_specs;

  AmanziSolvers::LinearOperatorFactory<Matrix_MFD, CompositeVector, CompositeVectorSpace> factory;
  Teuchos::RCP<AmanziSolvers::LinearOperator<Matrix_MFD, CompositeVector, CompositeVectorSpace> >
     solver = factory.Create(ls_specs.solver_name, linear_operator_list_, matrix_);

  CompositeVector& rhs = *matrix_->rhs();
  solver->ApplyInverse(rhs, *solution);

  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
    int num_itrs = solver->num_itrs();
    double residual = solver->residual();

    Teuchos::OSTab tab = vo_->getOSTab();
    *(vo_->os()) << "pressure solver (" << solver->name() 
                 << "): ||r||=" << residual << " itr=" << num_itrs << endl;
  }
}


/* ******************************************************************
* Calculates steady-state solution using a user-given rhs vector. 
* The matrix has to be assembled before this call.
****************************************************************** */
void Darcy_PK::SolveFullySaturatedProblem(double Tp, const CompositeVector& rhs, CompositeVector& u)
{
  LinearSolver_Specs& ls_specs = ti_specs->ls_specs;

  AmanziSolvers::LinearOperatorFactory<Matrix_MFD, CompositeVector, CompositeVectorSpace> factory;
  Teuchos::RCP<AmanziSolvers::LinearOperator<Matrix_MFD, CompositeVector, CompositeVectorSpace> >
     solver = factory.Create(ls_specs.solver_name, linear_operator_list_, matrix_);

  solver->ApplyInverse(rhs, u);

  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
    int num_itrs = solver->num_itrs();
    double residual = solver->residual();

    Teuchos::OSTab tab = vo_->getOSTab();
    *(vo_->os()) << "pressure solver: ||r||=" << residual << " itr=" << num_itrs << endl;
  }
}

}  // namespace AmanziFlow
}  // namespace Amanzi


