/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

#ifndef _NULLENERGYPK_HH_
#define _NULLENERGYPK_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_ParameterList.hpp"

#include "State.hh"
#include "Vector.hh"
#include "TreeVector.hh"
#include "PK.hh"

namespace Amanzi {

// Energy equation that simply maintains constant temperature
// Here for testing PK coupling.
class NullEnergyPK : public PK {

public:
  NullEnergyPK(Teuchos::ParameterList& energy_plist, Teuchos::RCP<State>& S,
               Teuchos::RCP<TreeVector>& soln);
  ~NullEnergyPK() {}

  // Initialize owned (dependent) variables.
  void initialize(Teuchos::RCP<State>& S, Teuchos::RCP<TreeVector>& soln);

  // transfer operators
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

  // Choose a time step compatible with physics.
  double get_dT() { return 1.e99; }

  // Advance from state S0 to state S1 at time S0.time + dt.
  bool advance(double dt, Teuchos::RCP<const State>& S0,
          Teuchos::RCP<State>& S1, Teuchos::RCP<TreeVector>& solution);

  void compute_f(const double t, const Vector& u,
                 const Vector& udot, Vector& f);

  // Commit any secondary (dependent) variables.
  void commit_state(double dt, Teuchos::RCP<State>& S) {}

  // update diagnostics for vis
  void calculate_diagnostics(Teuchos::RCP<State>& S) {}

private:
  // states
  Teuchos::RCP<State> S_;
  double T_;

  // misc setup information
  Teuchos::ParameterList energy_plist_;

};
} // namespace

#endif
