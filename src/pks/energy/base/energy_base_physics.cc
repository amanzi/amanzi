/* -*-  mode++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Solves:

de/dt + q dot grad h = div Ke grad T + S?
------------------------------------------------------------------------- */

#include "advection.hh"
#include "FieldEvaluator.hh"
#include "energy_base.hh"
#include "Op.hh"

namespace Amanzi {
namespace Energy {

// -------------------------------------------------------------
// Accumulation of energy term de/dt
// -------------------------------------------------------------
void EnergyBase::AddAccumulation_(const Teuchos::Ptr<CompositeVector>& g) {
  double dt = S_next_->time() - S_inter_->time();

  // update the energy at both the old and new times.
  S_next_->GetFieldEvaluator(energy_key_)->HasFieldChanged(S_next_.ptr(), name_);
  S_inter_->GetFieldEvaluator(energy_key_)->HasFieldChanged(S_inter_.ptr(), name_);

  // get the energy at each time
  Teuchos::RCP<const CompositeVector> e1 = S_next_->GetFieldData(energy_key_);
  Teuchos::RCP<const CompositeVector> e0 = S_inter_->GetFieldData(energy_key_);

  // Update the residual with the accumulation of energy over the
  // timestep, on cells.
  g->ViewComponent("cell", false)
    ->Update(1.0/dt, *e1->ViewComponent("cell", false),
          -1.0/dt, *e0->ViewComponent("cell", false), 1.0);
};


// -------------------------------------------------------------
// Advective term for transport of enthalpy, q dot grad h.
// -------------------------------------------------------------
void EnergyBase::AddAdvection_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& g, bool negate) {

  // set up the operator

  // NOTE: fluxes are a MOLAR flux by choice of the flow pk, i.e.
  // [flux] =  mol/s

  // NOTE: this will be the eventual way to ensure it is up to date,
  // but there is no FieldEvaluator for darcy flux yet.  When there
  // is, we can take the evaluation out of Flow::commit_state(),
  // but for now we'll leave it there and assume it has been updated. --etc
  //  S->GetFieldEvaluator(flux_key_)->HasFieldChanged(S.ptr(), name_);
  Teuchos::RCP<const CompositeVector> flux = S->GetFieldData(flux_key_);
  db_->WriteVector(" adv flux", flux.ptr(), true);
  matrix_adv_->global_operator()->Init();
  matrix_adv_->Setup(*flux);
  matrix_adv_->UpdateMatrices(*flux);

  // apply to enthalpy
  S->GetFieldEvaluator(enthalpy_key_)->HasFieldChanged(S.ptr(), name_);
  Teuchos::RCP<const CompositeVector> enth = S->GetFieldData(enthalpy_key_);;
  ApplyDirichletBCsToEnthalpy_(S.ptr());
  matrix_adv_->ApplyBCs(bc_adv_, true);

  // apply
  matrix_adv_->global_operator()->ComputeNegativeResidual(*enth, *g, false);
}

// -------------------------------------------------------------
// Diffusion term, div K grad T
// -------------------------------------------------------------
void EnergyBase::ApplyDiffusion_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& g) {
  // update the thermal conductivity
  UpdateConductivityData_(S_next_.ptr());
  Teuchos::RCP<const CompositeVector> conductivity =
      S_next_->GetFieldData(uw_conductivity_key_);

  // update the stiffness matrix
  matrix_diff_->global_operator()->Init();
  matrix_diff_->SetScalarCoefficient(conductivity, Teuchos::null);
  matrix_diff_->UpdateMatrices(Teuchos::null, Teuchos::null);
  Teuchos::RCP<const CompositeVector> temp = S->GetFieldData(key_);

  // update the flux if needed
  Teuchos::RCP<CompositeVector> flux = S->GetFieldData(energy_flux_key_, name_);
  matrix_diff_->UpdateFlux(*temp, *flux);

  // finish assembly of the stiffness matrix
  matrix_diff_->ApplyBCs(true, true);

  // calculate the residual
  matrix_diff_->global_operator()->ComputeNegativeResidual(*temp, *g);
};


// ---------------------------------------------------------------------
// Add in energy source, which are accumulated by a single evaluator.
// Note that that evaluator applies the factor of cell volume.
// ---------------------------------------------------------------------
void EnergyBase::AddSources_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& g) {
  Teuchos::OSTab tab = vo_->getOSTab();

  // external sources of energy
  if (is_source_term_) {
    Epetra_MultiVector& g_c = *g->ViewComponent("cell",false);

    // Update the source term
    S->GetFieldEvaluator(source_key_)->HasFieldChanged(S, name_);
    const Epetra_MultiVector& source1 =
        *S->GetFieldData(source_key_)->ViewComponent("cell",false);

    // Add into residual
    unsigned int ncells = g_c.MyLength();
    for (unsigned int c=0; c!=ncells; ++c) {
      g_c[0][c] -= source1[0][c];
    }

    if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
      *vo_->os() << "Adding external source term" << std::endl;
      db_->WriteVector("  Q_ext", S->GetFieldData(source_key_).ptr(), false);
      db_->WriteVector("res (src)", g, false);
    }
  }
}


void EnergyBase::AddSourcesToPrecon_(const Teuchos::Ptr<State>& S, double h) {
  // external sources of energy (temperature dependent source)
  if (is_source_term_ && S->GetFieldEvaluator(source_key_)->IsDependency(S, key_)) {
    std::vector<double>& Acc_cells = preconditioner_acc_->local_matrices()->vals;

    S->GetFieldEvaluator(source_key_)->HasFieldDerivativeChanged(S, name_, key_);
    const Epetra_MultiVector& dsource_dT =
        *S->GetFieldData(dsource_dT_key_)->ViewComponent("cell",false);
    unsigned int ncells = dsource_dT.MyLength();
    for (unsigned int c=0; c!=ncells; ++c) {
      Acc_cells[c] -= dsource_dT[0][c];
    }
  }
}

} //namespace Energy
} //namespace Amanzi
