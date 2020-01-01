/*
  MultiPhase PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Quan Bui (mquanbui@math.umd.edu)
           Konstantin Lipnikov (lipnikov@lanl.gov)

  Multiphase multi-component flow, see Multiphase_PK.cc for more detail.
*/


// TPLs
#include "Teuchos_RCP.hpp"

// Multiphase
#include "Multiphase_PK.hh"

namespace Amanzi {
namespace Multiphase {

/* ******************************************************************
* F(u) implementation
****************************************************************** */
void Multiphase_PK::FunctionalResidual(double t_old, double t_new, 
                                       Teuchos::RCP<TreeVector> u_old,
                                       Teuchos::RCP<TreeVector> u_new,
                                       Teuchos::RCP<TreeVector> f) 
{
  // comp_w_pk_->FunctionalResidual(t_old, t_new, u_old, u_new, f->SubVector(0));
  // comp_h_pk_->FunctionalResidual(t_old, t_new, u_old, u_new, f->SubVector(1));
  // gas_constraint_pk_->FunctionalResidual(t_old, t_new, u_old, u_new, f->SubVector(2));
  // rhs_ = f;
}


/* ******************************************************************
* Update new preconditioner on the interval (tp-dtp, tp].
****************************************************************** */
void Multiphase_PK::UpdatePreconditioner(double tp, Teuchos::RCP<const TreeVector> u, double dtp)
{
  double t_old = tp - dtp;
}


/* ******************************************************************
* Apply preconditioner inv(B) * X.                                                 
****************************************************************** */
int Multiphase_PK::ApplyPreconditioner(Teuchos::RCP<const TreeVector> X, 
                                       Teuchos::RCP<TreeVector> Y)
{
  Y->PutScalar(0.0);
  // return op_pc_solver_->ApplyInverse(*X->Data(), *Y->Data());
  return 0;
}


/* ******************************************************************
* Monitor l2 norm of residual
****************************************************************** */
double Multiphase_PK::ErrorNorm(Teuchos::RCP<const TreeVector> u,
                                Teuchos::RCP<const TreeVector> du) 
{
  double du_l2 = 0.0;
  double resnorm_p, resnorm_s, resnorm_r;

  du->SubVector(0)->Data()->Norm2(&resnorm_p);
  du->SubVector(1)->Data()->Norm2(&resnorm_s);
  du->SubVector(2)->Data()->Norm2(&resnorm_r);
  printf("resnorm_p = %4.6e, resnorm_s = %4.6e, resnorm_r = %4.6e \n", resnorm_p, resnorm_s, resnorm_r);

  du->Norm2(&du_l2);
  return du_l2;
}


/********************************************************************
* Modifies nonlinear update du using .. TBW
****************************************************************** */
AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
Multiphase_PK::ModifyCorrection(double h, Teuchos::RCP<const TreeVector> res,
                               Teuchos::RCP<const TreeVector> u,
                               Teuchos::RCP<TreeVector> du)
{
  Teuchos::RCP<CompositeVector> rho_next = Teuchos::rcp(new CompositeVector(*u->SubVector(2)->Data()));
  rho_next->Update(-1.0, *du->SubVector(2)->Data(), 1.0);
  ClipConcentration_(rho_next);
  du->SubVector(2)->Data()->Update(1.0, *u->SubVector(2)->Data(), -1.0, *rho_next, 0.0);

  return AmanziSolvers::FnBaseDefs::CORRECTION_MODIFIED;
}


/********************************************************************
* Helper function for modifying correction.
****************************************************************** */
void Multiphase_PK::ClipConcentration_(Teuchos::RCP<CompositeVector> rho)
{
  Epetra_MultiVector& rho_c = *rho->ViewComponent("cell");
  for (int c = 0; c < rho_c.MyLength(); c++) {
    rho_c[0][c] = std::max(0.0, rho_c[0][c]);
  }
}

}  // namespace Multiphase
}  // namespace Amanzi

