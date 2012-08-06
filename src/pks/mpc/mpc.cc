/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
   ATS

   License: see $ATS_DIR/COPYRIGHT
   Author: Ethan Coon

   Implementation for the Base MPC class.  A multi process coordinator is a PK
   (process kernel) which coordinates several PKs.  Each of these coordinated PKs
   may be MPCs themselves, or physical PKs.  Note this does NOT provide a full
   implementation of PK -- it does not supply the advance() method.  Therefore
   this class cannot be instantiated, but must be inherited by derived classes
   which finish supplying the functionality.  Instead, this provides the data
   structures and methods (which may be overridden by derived classes) for
   managing multiple PKs.

   Most of these methods simply loop through the coordinated PKs, calling their
   respective methods.
   ------------------------------------------------------------------------- */

#include <string>

#include "errors.hh"

#include "dbc.hh"
#include "mpc.hh"

namespace Amanzi {

// loop over sub-PKs, calling their initialization methods
void MPC::initialize(const Teuchos::RCP<State>& S) {
  for (std::vector< Teuchos::RCP<PK> >::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    (*pk)->initialize(S);
  }
};

// loop over sub-PKs, calling their state_to_solution method
void MPC::state_to_solution(const Teuchos::RCP<State>& S,
                            const Teuchos::RCP<TreeVector>& soln) {
  for (std::vector< Teuchos::RCP<PK> >::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    Teuchos::RCP<TreeVector> pk_soln = soln->SubVector((*pk)->name());
    if (pk_soln == Teuchos::null) {
      Errors::Message message("MPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    }
    (*pk)->state_to_solution(S, pk_soln);
  }
};

// loop over sub-PKs, calling their solution_to_state method
void MPC::solution_to_state(const Teuchos::RCP<TreeVector>& soln,
                            const Teuchos::RCP<State>& S) {
  for (std::vector< Teuchos::RCP<PK> >::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    Teuchos::RCP<TreeVector> pk_soln = soln->SubVector((*pk)->name());
    if (pk_soln == Teuchos::null) {
      Errors::Message message("MPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    }
    (*pk)->solution_to_state(pk_soln, S);
  }
};

// min(pk.get_dt()) for each sub-PK
double MPC::get_dt() {
  double dt = 1.e99;
  for (std::vector< Teuchos::RCP<PK> >::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    dt = std::min(dt, (*pk)->get_dt());
  }
  return dt;
}

// loop over sub-PKs, calling their commit_state method
void MPC::commit_state(double dt, const Teuchos::RCP<State>& S) {
  for (std::vector< Teuchos::RCP<PK> >::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    (*pk)->commit_state(dt, S);
  }
};

// loop over sub-PKs, calling their calculate_diagnostics method
void MPC::calculate_diagnostics(const Teuchos::RCP<State>& S) {
  for (std::vector< Teuchos::RCP<PK> >::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    (*pk)->calculate_diagnostics(S);
  }
};

// loop over sub-PKs, calling their set_states() methods
void MPC::set_states(const Teuchos::RCP<const State>& S,
                     const Teuchos::RCP<State>& S_inter,
                     const Teuchos::RCP<State>& S_next) {
  S_ = S;
  S_inter_ = S_inter;
  S_next_ = S_next;

  // do the loop
  for (std::vector< Teuchos::RCP<PK> >::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    (*pk)->set_states(S, S_inter, S_next);
  }
};
} // namespace
