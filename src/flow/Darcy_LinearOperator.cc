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
* Calculates steady-state solution assuming that absolute permeability 
* does not depend on time. The boundary conditions are calculated
* only once, during the initialization step.                                                
****************************************************************** */
void Darcy_PK::SolveFullySaturatedProblem(double Tp, CompositeVector& u)
{
  // calculate and assemble elemental stifness matrices
  matrix_->CreateStiffnessMatricesDarcy(mfd3d_method_);
  matrix_->CreateRHSVectors();
  matrix_->AddGravityFluxesDarcy(rho_, gravity_);
  matrix_->ApplyBoundaryConditions(bc_model, bc_values);
  matrix_->Assemble();
  matrix_->AssembleDerivatives(u, bc_model, bc_values);
  matrix_->UpdatePreconditioner();

  // create linear solver
  LinearSolver_Specs& ls_specs = ti_specs->ls_specs;

  AmanziSolvers::LinearOperatorFactory<FlowMatrix, CompositeVector, CompositeVectorSpace> factory;
  Teuchos::RCP<AmanziSolvers::LinearOperator<FlowMatrix, CompositeVector, CompositeVectorSpace> >
     solver = factory.Create(ls_specs.solver_name, linear_operator_list_, matrix_);

  CompositeVector& rhs = *matrix_->rhs();
  int ierr = solver->ApplyInverse(rhs, *solution);

  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
    int num_itrs = solver->num_itrs();
    double residual = solver->residual();
    int code = solver->returned_code();

    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "pressure solver (" << solver->name() 
               << "): ||r||=" << residual << " itr=" << num_itrs
               << " code=" << code << std::endl;
  }
  if (ierr != 0) {
    Errors::Message msg;
    msg << "\nLinear solver returned an unrecoverable error code.\n";
    Exceptions::amanzi_throw(msg);
  }
}

}  // namespace AmanziFlow
}  // namespace Amanzi


