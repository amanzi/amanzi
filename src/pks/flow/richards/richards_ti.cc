/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
A base two-phase, thermal Richard's equation with water vapor.

License: BSD
Authors: Neil Carlson (version 1) 
         Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
         Ethan Coon (ATS version) (ecoon@lanl.gov)
*/

#include "richards.hh"

namespace Amanzi {
namespace AmanziFlow {

/* ******************************************************************
 * Calculate f(u, du/dt) = d(s u)/dt + A*u - g.                                         
 ****************************************************************** */
void Richards::fun(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
                       Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> g) {
  S_inter_->set_time(t_old);
  S_next_->set_time(t_new);

  // pointer-copy temperature into states and update any auxilary data
  solution_to_state(u_old, S_inter_);
  solution_to_state(u_new, S_next_);
  UpdateSecondaryVariables_(S_inter_);
  UpdateSecondaryVariables_(S_next_);

  // update boundary conditions
  bc_pressure_->Compute(t_new);
  bc_head_->Compute(t_new);
  bc_flux_->Compute(t_new);
  UpdateBoundaryConditions_();

  // zero out residual
  Teuchos::RCP<CompositeVector> res = g->data();
  res->PutScalar(0.0);

  // diffusion term, treated implicitly
  ApplyDiffusion_(S_next_, res);

  // accumulation term
  AddAccumulation_(res);
};

/* ******************************************************************
* .                                                 
****************************************************************** */
// applies preconditioner to u and returns the result in Pu
void TwoPhase::precon(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu) {
  std::cout << "Precon application:" << std::endl;
  std::cout << "  u: " << (*u->data())("cell",0,0) << " " << (*u->data())("face",0,0) << std::endl;
  preconditioner_->ApplyInverse(*u->data(), Pu->data());
  std::cout << "  Pu: " << (*Pu->data())("cell",0,0) << " " << (*Pu->data())("face",0,0) << std::endl;
};

  
/* ******************************************************************
* Compute new preconditioner B(p, dT_prec). For BDF2 method, we need
* a separate memory allocation.                                              
****************************************************************** */
void TwoPhase::update_precon(double t, Teuchos::RCP<const TreeVector> up, double h) {
  S_next_->set_time(t);
  PK::solution_to_state(up, S_next_); // not sure why this isn't getting found? --etc
  UpdateSecondaryVariables_(S_next_);

  // upwind the rel perm into faces?
  
  // div K_e grad u
  Teuchos::RCP<CompositeVector> rel_perm =
    S_next_->GetFieldData("relative_permeability", "flow");

  // update boundary conditions
  bc_pressure_->Compute(S_next_->time());
  bc_head_->Compute(S_next_->time());
  bc_flux_->Compute(S_next_->time());
  UpdateBoundaryConditions_();

  preconditioner_->CreateMFDstiffnessMatrices(K_, *rel_perm);
  preconditioner_->CreateMFDrhsVectors();

  // update with gravity and accumulation terms
}

}  // namespace Flow
}  // namespace Amanzi



