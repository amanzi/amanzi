/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Implementation for the derived StrongMPC class.  Is both a PK and a Model
Evalulator, providing needed methods for BDF time integration of the coupled
system.

Completely automated and generic to any sub PKs, this uses a block diagonal
preconditioner.

See additional documentation in the base class src/pks/mpc/MPC.hh
------------------------------------------------------------------------- */

#include "bdf1_time_integrator.hh"

#include "strong_mpc.hh"

namespace Amanzi {

RegisteredPKFactory<StrongMPC> StrongMPC::reg_("strong MPC");


// -----------------------------------------------------------------------------
// Setup
// -----------------------------------------------------------------------------
void StrongMPC::setup(const Teuchos::Ptr<State>& S) {
  MPC<PKBDFBase>::setup(S);
  PKBDFBase::setup(S);
};


// -----------------------------------------------------------------------------
// Required unique initailize(), just calls both of its base class
// initialize() methods.
// -----------------------------------------------------------------------------
void StrongMPC::initialize(const Teuchos::Ptr<State>& S) {
  // Initialize all sub PKs.
  MPC<PKBDFBase>::initialize(S);

  // Initialize my timestepper.
  PKBDFBase::initialize(S);
};


// -----------------------------------------------------------------------------
// Compute the non-linear functional g = g(t,u,udot).
// -----------------------------------------------------------------------------
void StrongMPC::fun(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
                    Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> g) {
  // loop over sub-PKs
  for (MPC<PKBDFBase>::SubPKList::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {

    // pull out the old solution sub-vector
    Teuchos::RCP<TreeVector> pk_u_old = u_old->SubVector((*pk)->name());
    if (pk_u_old == Teuchos::null) {
      Errors::Message message("MPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    }

    // pull out the new solution sub-vector
    Teuchos::RCP<TreeVector> pk_u_new = u_new->SubVector((*pk)->name());
    if (pk_u_new == Teuchos::null) {
      Errors::Message message("MPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    }

    // pull out the residual sub-vector
    Teuchos::RCP<TreeVector> pk_g = g->SubVector((*pk)->name());
    if (pk_g == Teuchos::null) {
      Errors::Message message("MPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    }

    // fill the nonlinear function with each sub-PKs contribution
    (*pk)->fun(t_old, t_new, pk_u_old, pk_u_new, pk_g);
  }
};


// -----------------------------------------------------------------------------
// Applies preconditioner to u and returns the result in Pu.
// -----------------------------------------------------------------------------
void StrongMPC::precon(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu) {
  // loop over sub-PKs
  for (MPC<PKBDFBase>::SubPKList::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {

    // pull out the u sub-vector
    Teuchos::RCP<const TreeVector> pk_u = u->SubVector((*pk)->name());
    if (pk_u == Teuchos::null) {
      Errors::Message message("MPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    }

    // pull out the preconditioned u sub-vector
    Teuchos::RCP<TreeVector> pk_Pu = Pu->SubVector((*pk)->name());
    if (pk_Pu == Teuchos::null) {
      Errors::Message message("MPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    }

    // Fill the preconditioned u as the block-diagonal product using each sub-PK.
    (*pk)->precon(pk_u, pk_Pu);
  }
};


// -----------------------------------------------------------------------------
// Compute a norm on u-du and returns the result.
// For a Strong MPC, the enorm is just the max of the sub PKs enorms.
// -----------------------------------------------------------------------------
double StrongMPC::enorm(Teuchos::RCP<const TreeVector> u,
                        Teuchos::RCP<const TreeVector> du){
  double norm = 0.0;

  // loop over sub-PKs
  for (MPC<PKBDFBase>::SubPKList::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {

    // pull out the u sub-vector
    Teuchos::RCP<const TreeVector> pk_u = u->SubVector((*pk)->name());
    if (pk_u == Teuchos::null) {
      Errors::Message message("MPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    }

    // pull out the du sub-vector
    Teuchos::RCP<const TreeVector> pk_du = du->SubVector((*pk)->name());
    if (pk_du == Teuchos::null) {
      Errors::Message message("MPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    }

    // norm is the max of the sub-PK norms
    double tmp_norm = (*pk)->enorm(pk_u, pk_du);
    norm = std::max(norm, tmp_norm);
  }
  return norm;
};


// -----------------------------------------------------------------------------
// Update the preconditioner.
// -----------------------------------------------------------------------------
void StrongMPC::update_precon(double t, Teuchos::RCP<const TreeVector> up, double h) {
  PKDefaultBase::solution_to_state(up, S_next_);

  // loop over sub-PKs
  for (MPC<PKBDFBase>::SubPKList::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {

    // pull out the up sub-vector
    Teuchos::RCP<const TreeVector> pk_up = up->SubVector((*pk)->name());
    if (pk_up == Teuchos::null) {
      Errors::Message message("MPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    }

    // update precons of each of the sub-PKs
    (*pk)->update_precon(t, pk_up, h);
  };
};


// -----------------------------------------------------------------------------
// Experimental approach -- calling this indicates that the time integration
// scheme is changing the value of the solution in state.
// -----------------------------------------------------------------------------
void StrongMPC::changed_solution() {
  // loop over sub-PKs
  for (MPC<PKBDFBase>::SubPKList::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    (*pk)->changed_solution();
  }
};

} // namespace
