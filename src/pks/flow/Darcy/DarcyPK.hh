/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

#ifndef PKS_FLOW_DARCY_DARCYPK_HH_
#define PKS_FLOW_DARCY_DARCYPK_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Epetra_Vector.h"
#include "AztecOO.h"

#include "PK.hh"
#include "State.hh"
#include "DarcyProblem.hh"

namespace Amanzi {

class DarcyPK : public PK {

public:
  // Populate state with pressure, porosity, perm
  DarcyPK(Teuchos::ParameterList& plist, Teuchos::RCP<State>& S,
           Teuchos::RCP<TreeVector>& soln);

  // Initialize will do all the work, running to steady state.
  void initialize(Teuchos::RCP<State>& S, Teuchos::RCP<TreeVector>& soln);

  // transfer operators
  virtual void state_to_solution(const State& S, Teuchos::RCP<TreeVector>& soln);
  virtual void state_to_solution(const State& S, Teuchos::RCP<TreeVector>& soln,
          Teuchos::RCP<TreeVector>& soln_dot);
  virtual void solution_to_state(const TreeVector& soln, Teuchos::RCP<State>& S);
  virtual void solution_to_state(const TreeVector& soln, const TreeVector& soln_dot,
          Teuchos::RCP<State>& S);

  // Steady state Darcy -- no timestep limitations
  double get_dt() { return 1.e99; }

  // Steady state Darcy -- transient does nothing.
  bool advance(double dt, Teuchos::RCP<TreeVector>& solution) { return 0; }

  // Transient does nothing, so commit not needed.
  void commit_state(double dt, Teuchos::RCP<State>& S) {}

  // Calculate flow on cells for diagnostic purposes.
  void calculate_diagnostics(Teuchos::RCP<State>& S);

private:
  void calculate_velocity(Teuchos::RCP<Epetra_MultiVector>& q) const;
  void advance_to_steady_state(Teuchos::RCP<TreeVector>& soln);

  Teuchos::RCP<DarcyProblem> problem_;
  Teuchos::RCP<AztecOO> solver_;

  int max_itr_;      // max number of linear solver iterations
  double err_tol_;   // linear solver convergence error tolerance

  Teuchos::ParameterList flow_plist_;
};

} //namespace

#endif
