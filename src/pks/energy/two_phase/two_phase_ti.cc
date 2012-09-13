/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon
------------------------------------------------------------------------- */

#include "boost/math/special_functions/fpclassify.hpp"

#include "boundary_function.hh"
#include "field_evaluator.hh"
#include "two_phase.hh"

namespace Amanzi {
namespace Energy {

#define DEBUG_FLAG 1

// TwoPhase is a BDFFnBase
// -----------------------------------------------------------------------------
// computes the non-linear functional g = g(t,u,udot)
// -----------------------------------------------------------------------------
void TwoPhase::fun(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
                       Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> g) {
  // VerboseObject stuff.
  Teuchos::OSTab tab = getOSTab();

  S_inter_->set_time(t_old);
  S_next_->set_time(t_new);

  Teuchos::RCP<CompositeVector> u = u_new->data();
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
    *out_ << "Two-Phase Residual calculation:" << std::endl;
    *out_ << "  T0: " << (*u)("cell",0) << " " << (*u)("face",3) << std::endl;
    *out_ << "  T1: " << (*u)("cell",99) << " " << (*u)("face",497) << std::endl;
  }

  // pointer-copy temperature into states and update any auxilary data
  solution_to_state(u_new, S_next_);

  // update boundary conditions
  bc_temperature_->Compute(t_new);
  bc_flux_->Compute(t_new);
  UpdateBoundaryConditions_();

  // zero out residual
  Teuchos::RCP<CompositeVector> res = g->data();
  res->PutScalar(0.0);

  // diffusion term, implicit
  ApplyDiffusion_(S_next_, res);
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
    *out_ << "  res0 (after diffusion): " << (*res)("cell",0) << " " << (*res)("face",3) << std::endl;
    *out_ << "  res1 (after diffusion): " << (*res)("cell",99) << " " << (*res)("face",497) << std::endl;
  }

  // accumulation term
  AddAccumulation_(res);
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
    *out_ << "  res0 (after accumulation): " << (*res)("cell",0) << " " << (*res)("face",3) << std::endl;
    *out_ << "  res1 (after accumulation): " << (*res)("cell",99) << " " << (*res)("face",497) << std::endl;
  }

  // advection term, implicit
  AddAdvection_(S_next_, res, false);
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
    *out_ << "  res0 (after advection): " << (*res)("cell",0) << " " << (*res)("face",3) << std::endl;
    *out_ << "  res1 (after advection): " << (*res)("cell",99) << " " << (*res)("face",497) << std::endl;
  }
};


// -----------------------------------------------------------------------------
// Apply the preconditioner to u and return the result in Pu.
// -----------------------------------------------------------------------------
void TwoPhase::precon(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu) {
  // VerboseObject stuff.
  Teuchos::OSTab tab = getOSTab();

  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
    *out_ << "Precon application:" << std::endl;
    *out_ << "  T0: " << (*u->data())("cell",0) << " " << (*u->data())("face",3) << std::endl;
    *out_ << "  T1: " << (*u->data())("cell",99) << " " << (*u->data())("face",497) << std::endl;
  }

  preconditioner_->ApplyInverse(*u->data(), Pu->data());

  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
    *out_ << "  PC*T0: " << (*Pu->data())("cell",0) << " " << (*Pu->data())("face",3) << std::endl;
    *out_ << "  PC*T1: " << (*Pu->data())("cell",99) << " " << (*Pu->data())("face",497) << std::endl;
  }
};


// -----------------------------------------------------------------------------
// Update the preconditioner at time t and u = up
// -----------------------------------------------------------------------------
void TwoPhase::update_precon(double t, Teuchos::RCP<const TreeVector> up, double h) {
  // VerboseObject stuff.
  Teuchos::OSTab tab = getOSTab();

  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
    *out_ << "Precon update at t = " << t << std::endl;
  }

  // update state with the solution up.
  S_next_->set_time(t);
  PKDefaultBase::solution_to_state(up, S_next_);

  // update boundary conditions
  bc_temperature_->Compute(S_next_->time());
  bc_flux_->Compute(S_next_->time());
  UpdateBoundaryConditions_();

  // div K_e grad u
  S_next_->GetFieldEvaluator("thermal_conductivity")
    ->HasFieldChanged(S_next_.ptr(), "energy_pk");
  Teuchos::RCP<const CompositeVector> thermal_conductivity =
    S_next_->GetFieldData("thermal_conductivity");

  preconditioner_->CreateMFDstiffnessMatrices(thermal_conductivity);
  preconditioner_->CreateMFDrhsVectors();

  // update with accumulation terms
  // -- update the accumulation derivatives
  S_next_->GetFieldEvaluator("energy")
      ->HasFieldDerivativeChanged(S_next_.ptr(), "energy_pk", "temperature");

  // -- get the accumulation deriv
  Teuchos::RCP<const CompositeVector> de_dT =
      S_next_->GetFieldData("denergy_dtemperature");
  Teuchos::RCP<const CompositeVector> temp =
      S_next_->GetFieldData("temperature");

  // -- get the matrices/rhs that need updating
  std::vector<double>& Acc_cells = preconditioner_->Acc_cells();
  std::vector<double>& Fc_cells = preconditioner_->Fc_cells();

  // -- update the diagonal
  int ncells = temp->size("cell");
  for (int c=0; c!=ncells; ++c) {
    Acc_cells[c] += (*de_dT)("cell",c) / h;
    //    Fc_cells[c] += (*de_dT)("cell",c) / h * (*temp)("cell",c);
  }

  // -- assemble
  preconditioner_->ApplyBoundaryConditions(bc_markers_, bc_values_);
  preconditioner_->AssembleGlobalMatrices();

  // -- form and prep the Schur complement for inversion
  preconditioner_->ComputeSchurComplement(bc_markers_, bc_values_);
  preconditioner_->UpdatePreconditioner();

};

} // namespace Energy
} // namespace Amanzi
