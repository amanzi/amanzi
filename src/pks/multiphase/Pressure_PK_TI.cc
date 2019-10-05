/*
  MultiPhase

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Quan Bui (mquanbui@math.umd.edu)

  This routine implements the interface functions for the implicit
  time integrator.
*/

#include "Pressure_PK.hh"
#include "Epetra_Vector.h"

namespace Amanzi {
namespace Multiphase {

/* ******************************************************************
* Calculate f(u, du/dt) = a d(s(u))/dt + A*u - rhs.
This is basically the residual
****************************************************************** */
void Pressure_PK::Functional(double t_old, double t_new, 
						  Teuchos::RCP<TreeVector> u_old,
			    		Teuchos::RCP<TreeVector> u_new, 
			    		Teuchos::RCP<TreeVector> f)
{
  op_matrix_->ApplyBCs(true, true, true);
  op_matrix_->global_operator()->SymbolicAssembleMatrix();
  op_matrix_->global_operator()->AssembleMatrix();

  // Update source term
  Teuchos::RCP<CompositeVector> rhs = op_matrix_->global_operator()->rhs();
  //if (src_sink_ != NULL) AddSourceTerms(*rhs);
  
  op_matrix_->global_operator()->ComputeNegativeResidual(*u_new->Data(), *f->Data());
}


/* ******************************************************************
* Apply preconditioner inv(B) * X.                                                 
****************************************************************** */
int Pressure_PK::ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, 
                                		Teuchos::RCP<TreeVector> Pu)
{
  AmanziSolvers::LinearOperatorFactory<Operators::Operator, CompositeVector, CompositeVectorSpace> factory;
  Teuchos::RCP<AmanziSolvers::LinearOperator<Operators::Operator, CompositeVector, CompositeVectorSpace> >
     solver = factory.Create(solver_name_, *linear_operator_list_, op_);

  solver->ApplyInverse(*u->Data(), *Pu->Data());

  return 0;
}


/* ******************************************************************
* Update new preconditioner B(p, dT_prec).                                   
****************************************************************** */
void Pressure_PK::UpdatePreconditioner(double Tp, Teuchos::RCP<const TreeVector> u, double dTp)
{
  // Calculate total mobility needed to initialize diffusion operator
  Teuchos::RCP<const CompositeVector> water_saturation = S_->GetFieldData("water_saturation");
  //std::cout << "Pressure_PK::UpdatePreconditioner::water_saturation: " << *water_saturation->ViewComponent("cell") << "\n";
  rel_perm_w_->Compute(*water_saturation);
  rel_perm_n_->Compute(*water_saturation);
  *tot_mobility_ = *rel_perm_w_->Krel();
  tot_mobility_->Update(1.0/mu_[1], *rel_perm_n_->Krel(), 1.0/mu_[0]);
  //std::cout << "Pressure_PK::UpdatePreconditioner::tot_mobility_: " << *tot_mobility_->ViewComponent("cell") << "\n";

  if (include_gravity_) {
    // create scaling factor for gravity term
    gravity_factor_->PutScalar(0.0);
    gravity_factor_->Update(rho_[0]/mu_[0], *rel_perm_w_->Krel(), rho_[1]/mu_[1], *rel_perm_n_->Krel(), 0.0);
    //gravity_factor_->ReciprocalMultiply(1.0, *tot_mobility_, *gravity_factor_, 0.0);
  }
  //std::cout << "gravity factor: " << *gravity_factor_->ViewComponent("cell") << "\n";
  SetAbsolutePermeabilityTensor(gravity_factor_);
  Teuchos::RCP<std::vector<WhetStone::Tensor> > Kptr = Teuchos::rcpFromRef(K_);
  op_->Init();
  op_preconditioner_->Setup(Kptr, Teuchos::null, Teuchos::null);
  op_preconditioner_->UpdateMatrices(Teuchos::null, Teuchos::null);
  op_preconditioner_->ApplyBCs(true, true, true);

  rhs_ = Teuchos::rcp(new CompositeVector(*op_->rhs()));
  //if (src_sink_ != NULL) AddSourceTerms(*rhs_);

  SetAbsolutePermeabilityTensor(tot_mobility_);
  
  op_->Init();
  op_preconditioner_->Setup(Kptr, Teuchos::null, Teuchos::null);
  op_preconditioner_->UpdateMatrices(Teuchos::null, Teuchos::null);
  op_preconditioner_->ApplyBCs(true, true, true);

  op_->SymbolicAssembleMatrix();
  op_->AssembleMatrix();
  op_->InitPreconditioner(preconditioner_name_, *preconditioner_list_);
}

}  // namespace Multiphase
}  // namespace Amanzi
