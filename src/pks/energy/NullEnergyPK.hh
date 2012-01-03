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

namespace Amanzi {

class NullEnergyPK : public PK {

public:
  NullEnergyPK(Teuchos::ParameterList& energy_plist, Teuchos::RCP<State>& S,
               Teuchos::RCP<TreeVector>& soln);
  ~NullEnergyPK() {}

  // PK interface
  // -- Initialize owned (dependent) variables.
  void initialize(Teuchos::RCP<State>& S, Teuchos::RCP<TreeVector>& soln);

  // -- transfer operators
  void state_to_solution(const State& S,
                         Teuchos::RCP<TreeVector>& soln);
  void state_to_solution(const State& S,
                         Teuchos::RCP<TreeVector>& soln,
                         Teuchos::RCP<TreeVector>& soln_dot);
  void solution_to_state(const TreeVector& soln,
                         Teuchos::RCP<State>& S);
  void solution_to_state(const TreeVector& soln,
                         const TreeVector& soln_dot,
                         Teuchos::RCP<State>& S);

  // -- Choose a time step compatible with physics.
  double get_dt() { return 1.e99; }

  // -- Advance from state S0 to state S1 at time S0.time + dt.
  bool advance(double dt, Teuchos::RCP<TreeVector>& solution);

  // -- Commit any secondary (dependent) variables.
  void commit_state(double dt, Teuchos::RCP<State>& S) {}

  // -- Update diagnostics for vis.
  void calculate_diagnostics(Teuchos::RCP<State>& S) {}

  // extras, will eventually be put into a model evaluator.
  void compute_f(const double t, const Vector& u,
                 const Vector& udot, Vector& f);

private:
  // states
  double T_;

  // misc setup information
  Teuchos::ParameterList energy_plist_;
};
} // namespace

#endif
