/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Derived MPC for flow and energy.  This couples using a block-diagonal coupler.
------------------------------------------------------------------------- */

#include "field_evaluator.hh"
#include "mpc_diagonal_flow_energy.hh"


namespace Amanzi {

RegisteredPKFactory<MPCDiagonalFlowEnergy> MPCDiagonalFlowEnergy::reg_("energy-flow diagonally coupled");


MPCDiagonalFlowEnergy::MPCDiagonalFlowEnergy(Teuchos::ParameterList& mpc_plist,
        const Teuchos::RCP<State>& S, const Teuchos::RCP<TreeVector>& soln) :
    StrongMPC(mpc_plist, S, soln) {};


void MPCDiagonalFlowEnergy::initialize(const Teuchos::RCP<State>& S) {
  StrongMPC::initialize(S);

  Teuchos::RCP<const CompositeVector> pres = S->GetFieldData("pressure");
  Teuchos::RCP<const CompositeVector> temp = S->GetFieldData("temperature");

  ASSERT(temp->size("cell") == pres->size("cell"));
  D_pT_ = Teuchos::rcp(new Epetra_MultiVector(*pres->ViewComponent("cell", false)));
  //  D_Tp_ = Teuchos::rcp(new Epetra_MultiVector(*D_pT_));
};


void MPCDiagonalFlowEnergy::update_precon(double t, Teuchos::RCP<const TreeVector> up,
        double h) {
  StrongMPC::update_precon(t, up, h);

  // d pressure_residual / d T
  // -- update the accumulation derivative
  S_next_->GetFieldEvaluator("water_content")
      ->HasFieldDerivativeChanged(S_next_.ptr(), "richards_pk", "temperature");

  // -- get the accumulation deriv
  Teuchos::RCP<const CompositeVector> dwc_dT =
      S_next_->GetFieldData("dwater_content_dtemperature");
  for (int c=0; c!=dwc_dT->size("cell"); ++c) {
    (*D_pT_)[0][c] = (*dwc_dT)("cell",c) / h;
  }

  // // d temperature_residual / d p
  // // -- update the accumulation derivative
  // S_next_->GetFieldEvaluator("energy")
  //     ->HasFieldDerivativeChanged(S_next_.ptr(), "energy_pk", "pressure");

  // // -- get the accumulation deriv
  // Teuchos::RCP<const CompositeVector> de_dp =
  //     S_next_->GetFieldData("denergy_dpressure");
  // for (int c=0; c!=de_dp->size("cell"); ++c) {
  //   (*D_Tp_)[0][c] = (*de_dp)("cell",c) / h;
  // }

};


// -----------------------------------------------------------------------------
// Apply the preconditioner to u and return the result in Pu.
// -----------------------------------------------------------------------------
void MPCDiagonalFlowEnergy::precon(Teuchos::RCP<const TreeVector> u,
        Teuchos::RCP<TreeVector> Pu) {

  // pull out the temperature sub-vector
  Teuchos::RCP<const TreeVector> uT = u->SubVector("energy");
  Teuchos::RCP<TreeVector> PuT = Pu->SubVector("energy");

  // apply temperature inverse
  sub_pks_[1]->precon(uT, PuT);

  // subtract off the d pressure residual dT part.
  Teuchos::RCP<TreeVector> up_minusT = Teuchos::rcp(new TreeVector(*u->SubVector("flow")));
  *up_minusT = *u->SubVector("flow");
  up_minusT->data()->ViewComponent("cell",false)->Multiply(-1.0, *PuT->data()->ViewComponent("cell", false), *D_pT_, 1.0);

  // Apply the pressure inverse
  Teuchos::RCP<TreeVector> Pup = Pu->SubVector("flow");
  sub_pks_[0]->precon(up_minusT, Pup);
};

} //namespace
