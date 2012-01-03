/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Interface for the NullEnergy PK.  This PK simply provides a constant
temperature, and is provided for testing with other PKs that depend upon an
energy equation.  This could easily be provided by the state as an independent
variable, but this is nice for testing the full hierarchy with a simple PK.

Example usage:

  <ParameterList name="energy">
    <Parameter name="PK model" type="string" value="Constant Temperature"/>
    <Parameter name="Constant Temperature" type="double" value="290.0"/>
  </ParameterList>

------------------------------------------------------------------------- */

#ifndef PKS_ENERGY_NULLENERGYPK_HH_
#define PKS_ENERGY_NULLENERGYPK_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_ParameterList.hpp"

#include "State.hh"
#include "Vector.hh"
#include "TreeVector.hh"

#include "PK.hh"
#include "BDF2_fnBase.hpp"

namespace Amanzi {

class NullEnergyPK : public PK {
  //class NullEnergyPK : public PK, public BDF2::fnBase {

public:
  NullEnergyPK(Teuchos::ParameterList& energy_plist, Teuchos::RCP<State>& S,
               Teuchos::RCP<TreeVector>& soln);
  ~NullEnergyPK() {}

  // PK interface
  // -- Initialize owned (dependent) variables.
  void initialize(Teuchos::RCP<State>& S);

  // -- transfer operators
  void state_to_solution(State& S, Teuchos::RCP<TreeVector>& soln);
  void state_to_solution(State& S, Teuchos::RCP<TreeVector>& soln,
                         Teuchos::RCP<TreeVector>& soln_dot);
  void solution_to_state(TreeVector& soln, Teuchos::RCP<State>& S);
  void solution_to_state(TreeVector& soln, TreeVector& soln_dot,
                         Teuchos::RCP<State>& S);

  // -- Choose a time step compatible with physics.
  double get_dt() { return 1.e99; }

  // -- Advance from state S to state S_next at time S0.time + dt.
  bool advance(double dt);

  // -- Commit any secondary (dependent) variables.
  void commit_state(double dt, Teuchos::RCP<State>& S) {}

  // -- Update diagnostics for vis.
  void calculate_diagnostics(Teuchos::RCP<State>& S) {}

  // -- Set states for use in current and next timestep, and update solution from S_next
  void set_states(Teuchos::RCP<const State>& S, Teuchos::RCP<State>& S_next);

  // BDF2 interface
  // computes the non-linear functional f = f(t,u,udot)
  void fun(const double t, const TreeVector& soln, const TreeVector& soln_dot,
           TreeVector& f);

  // applies preconditioner to u and returns the result in Pu
  void precon(const TreeVector& u, TreeVector& Pu);

  // computes a norm on u-du and returns the result
  double enorm(const TreeVector& u, const TreeVector& du);

  // updates the preconditioner
  void update_precon(const double t, const TreeVector& up, const double h, int& errc) {};

  // check the admissibility of a solution
  // override with the actual admissibility check
  bool is_admissible(const TreeVector& up) { return true; }

private:
  // A few options for advance
  bool advance_analytic(double dt);
  bool advance_bdf(double dt);

  // states
  double T_;

  // misc setup information
  Teuchos::ParameterList energy_plist_;

  // for the model evaluator component
  double atol_;
  double rtol_;
};
} // namespace

#endif
