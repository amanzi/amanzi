/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
   ATS

   License: see $ATS_DIR/COPYRIGHT
   Author: Ethan Coon

   Implementation for a basic Richards PK.

   Example usage:

   <ParameterList name="flow">
   <Parameter name="PK model" type="string" value="Richards"/>
   </ParameterList>

   ------------------------------------------------------------------------- */

#ifndef PKS_FLOW_RICHARDS_RICHARDSPK_HH_
#define PKS_FLOW_RICHARDS_RICHARDSPK_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "PK.hh"
#include "State.hh"
#include "RichardsProblem.hh"
#include "BDF2_Dae.hpp"

namespace Amanzi {

class RichardsPK : public PK, public ImplicitTIBDF2fnBase {

public:
  RichardsPK(Teuchos::ParameterList& flow_plist, Teuchos::RCP<State>& S,
             Teuchos::RCP<TreeVector>& solution);

  void initialize(Teuchos::RCP<State>& S);


  // -- transfer operators
  void state_to_solution(Teuchos::RCP<State>& S, Teuchos::RCP<TreeVector>& soln);
  void state_to_solution(Teuchos::RCP<State>& S, Teuchos::RCP<TreeVector>& soln,
                         Teuchos::RCP<TreeVector>& soln_dot);
  void solution_to_state(Teuchos::RCP<TreeVector>& soln, Teuchos::RCP<State>& S);
  void solution_to_state(Teuchos::RCP<TreeVector>& soln, Teuchos::RCP<TreeVector>& soln_dot,
                         Teuchos::RCP<State>& S);

  // -- Choose a time step compatible with physics.
  double get_dt() { return h_; }

  // -- Advance from state S to state S_next at time S0.time + dt.
  bool advance(double dt);

  // -- Commit any secondary (dependent) variables.
  void commit_state(double dt, Teuchos::RCP<State>& S) {}

  // -- Update diagnostics for vis.
  void calculate_diagnostics(Teuchos::RCP<State>& S);

  // BDF2 interface
  // computes the non-linear functional f = f(t,u,udot)
  void fun(double t, Teuchos::RCP<const TreeVector> soln,
           Teuchos::RCP<const TreeVector> soln_dot, Teuchos::RCP<TreeVector> f);

  // applies preconditioner to u and returns the result in Pu
  void precon(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu);

  // computes a norm on u-du and returns the result
  double enorm(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<const TreeVector> du);

  // updates the preconditioner
  void update_precon(double t, Teuchos::RCP<const TreeVector> up, double h, int* errc);

  // check the admissibility of a solution
  // override with the actual admissibility check
  bool is_admissible(Teuchos::RCP<const TreeVector> up) { return true; }

private:
  // TODO: these should be scoped pointers
  Teuchos::RCP<RichardsProblem> problem_;
  Teuchos::RCP<Amanzi::ImplicitTIBDF2> time_stepper_;
  double atol_;
  double rtol_;

  Teuchos::ParameterList flow_plist_;

  double h0_; // initial timestep size
  double h_, hnext_; // current, next step sizes
};

} // close namespace Amanzi

#endif
