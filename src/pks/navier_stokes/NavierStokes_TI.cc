/*
  Navier Stokes PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include "NavierStokes_PK.hh"

namespace Amanzi {
namespace NavierStokes {

/* ******************************************************************
* Calculate f(u, du/dt) = a d(s(u))/dt + A*u - rhs.
****************************************************************** */
void NavierStokes_PK::Functional(double t_old, double t_new, 
                                 Teuchos::RCP<TreeVector> u_old,
                                 Teuchos::RCP<TreeVector> u_new, 
                                 Teuchos::RCP<TreeVector> f)
{ 
  double dtp = t_new - t_old;

  Teuchos::RCP<CompositeVector> uu = u_old->SubVector(0)->Data();
  Teuchos::RCP<CompositeVector> fu = f->SubVector(0)->Data();
  Teuchos::RCP<CompositeVector> fp = f->SubVector(1)->Data();

  // refresh data
  UpdateSourceBoundaryData_(t_old, t_new);

  // assemble residual using linear operator
  op_matrix_elas_->global_operator()->Init();
  op_matrix_elas_->UpdateMatrices();
  op_matrix_elas_->ApplyBCs(true, true);

  CompositeVector one(*uu);
  one.PutScalar(1.0);  // FIXME
  op_acc_->AddAccumulationDelta(*uu, one, dtp, "node");
  op_acc_->ApplyBCs(bcv_);

  op_div_->UpdateMatrices(*u_old->SubVector(0)->Data());
  op_div_->ApplyBCs(false, true);

  // Teuchos::RCP<CompositeVector> rhs = op_matrix_->rhs();
  // AddSourceTerms(*rhs);

  // compute negative residual
  op_matrix_->Apply(*u_new, *f);
  fu->Update(-1.0, *op_matrix_elas_->global_operator()->rhs(), 1.0);
  fp->Update(-1.0, *op_div_->global_operator()->rhs(), 1.0);

  // add accumulation term 

  // add convection term
}


/* ******************************************************************
* Apply preconditioner inv(B) * X.                                                 
****************************************************************** */
int NavierStokes_PK::ApplyPreconditioner(Teuchos::RCP<const TreeVector> X, 
                                         Teuchos::RCP<TreeVector> Y)
{
  Y->PutScalar(0.0);
  return op_pc_solver_->ApplyInverse(*X, *Y);
}


/* ******************************************************************
* Update new preconditioner on the interval (tp-dtp, tp].
****************************************************************** */
void NavierStokes_PK::UpdatePreconditioner(double tp, Teuchos::RCP<const TreeVector> u, double dtp)
{
  double t_old = tp - dtp;

  // refresh data
  UpdateSourceBoundaryData_(t_old, tp);

  // populate elastic operator
  Teuchos::RCP<Operators::Operator> global_op = op_preconditioner_elas_->global_operator();
  global_op->Init();
  op_preconditioner_elas_->UpdateMatrices();
  op_preconditioner_elas_->ApplyBCs(true, true);
  global_op->AssembleMatrix();
  global_op->InitPreconditioner(preconditioner_name_, *preconditioner_list_);

  // populate pressure operator
  global_op = op_mass_->global_operator();
  global_op->AssembleMatrix();
  global_op->InitPreconditioner("Diagonal", *preconditioner_list_);

  // add time derivative
  if (dtp > 0.0) {
    // op_acc_->AddAccumulationTerm(*u->Data(), dwc_dp, dtp, "cell");
  }

  // finalize global preconditioner
  op_preconditioner_->InitBlockDiagonalPreconditioner();
}


/* ******************************************************************
* Check difference du between the predicted and converged solutions.
* This is a wrapper for various error control methods. 
****************************************************************** */
double NavierStokes_PK::ErrorNorm(Teuchos::RCP<const TreeVector> u, 
                                  Teuchos::RCP<const TreeVector> du)
{
}
 

/* ******************************************************************
* A wrapper for updating boundary conditions.
****************************************************************** */
void NavierStokes_PK::UpdateSourceBoundaryData_(double t_old, double t_new)
{
  /*
  for (int i = 0; i < srcs.size(); ++i) {
    srcs[i]->Compute(t_old, t_new);
  }
  */

  for (int i = 0; i < bcs_.size(); i++) {
    bcs_[i]->Compute(t_old, t_new);
    bcs_[i]->ComputeSubmodel(mesh_);
  }

  // populate global arrays with data
  // ComputeOperatorBCs();
}

}  // namespace NavierStokes
}  // namespace Amanzi
 
