/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Virtual interface for Process Kernels.  All physical kernels and MPCs
must implement this interface for use within weak and strongly coupled
hierarchies.
------------------------------------------------------------------------- */


#ifndef PKS_PK_HH_
#define PKS_PK_HH_

#include "Teuchos_RCP.hpp"
#include "PK.hh"
#include "State.hh"
#include "TreeVector.hh"

namespace Amanzi {

class PK_ATS : public PK{
  //class PK_ATS {

public:
  // Virtual destructor
  virtual ~PK_ATS() {}

  // -- Setup
  virtual void setup(const Teuchos::Ptr<State>& S) = 0;

  // -- Initialize owned (dependent) variables.
  virtual void initialize(const Teuchos::Ptr<State>& S) = 0;

  // -- transfer operators
  virtual void state_to_solution(const Teuchos::RCP<State>& S,
                                 TreeVector& soln) = 0;
  virtual void solution_to_state(TreeVector& soln,
                                 const Teuchos::RCP<State>& S) = 0;

  // -- Choose a time step compatible with physics.
  virtual double get_dt() = 0;

  // -- Advance from state S0 to state S1 at time S0.time + dt.
  virtual bool advance(double dt) = 0;

  // -- Commit any secondary (dependent) variables.
  virtual void commit_state(double dt, const Teuchos::RCP<State>& S) = 0;

  // -- Calculate any diagnostics prior to doing vis
  virtual void calculate_diagnostics(const Teuchos::RCP<State>& S) = 0;

  // -- set pointers to State, and point the solution vector to the data in S_next
  virtual void set_states(const Teuchos::RCP<const State>& S,
                          const Teuchos::RCP<State>& S_inter,
                          const Teuchos::RCP<State>& S_next) = 0;

  virtual std::string name() = 0;

  virtual void set_dt(double dt) {};

  virtual void Setup(const Teuchos::Ptr<State>& S){
    setup(S);
  };
  virtual void Initialize(const Teuchos::Ptr<State>& S){
    initialize(S);
  };
  virtual void State_to_Solution(const Teuchos::RCP<State>& S,
                                 TreeVector& soln){
    state_to_solution(S, soln);
  };
  virtual void Solution_to_State(TreeVector& soln,
                                 const Teuchos::RCP<State>& S){
    solution_to_state(soln, S);
  };
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit){
    advance(t_new - t_old);
  };
  virtual void CommitStep(double t_old, double t_new, 
                          const Teuchos::RCP<State>& S) {
    commit_state(t_new - t_old, S);
  };

  virtual void CalculateDiagnostics(const Teuchos::RCP<State>& S) {
    calculate_diagnostics(S);
  };


};

} // namespace

#endif
