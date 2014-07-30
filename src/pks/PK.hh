/* 
  ATS & ASCEM

  Licenses: see $ATS_DIR/COPYRIGHT, $ASCEM_DIR/COPYRIGHT
  Author: Ethan Coon

  Virtual interface for Process Kernels. 
*/

#ifndef PK_HH_
#define PK_HH_

#include "Teuchos_RCP.hpp"

#include "State.hh"
#include "TreeVector.hh"

namespace Amanzi {

class PK {
 public:
  // Virtual destructor
  virtual ~PK() {}

  // Setup
  virtual void Setup(const Teuchos::Ptr<State>& S) = 0;

  // Initialize owned (dependent) variables.
  virtual void Initialize(const Teuchos::Ptr<State>& S) = 0;

  // Transfer operators
  virtual void StateToSolution(const Teuchos::RCP<State>& S,
                               const Teuchos::RCP<TreeVector>& soln) = 0;
  virtual void SolutionToState(const Teuchos::RCP<TreeVector>& soln,
                               const Teuchos::RCP<State>& S) = 0;

  // Choose a time step compatible with physics.
  virtual double GetDt() = 0;

  // Advance from state S0 to state S1 at time S0.time + dt.
  // Due to Flow PK / MPC conflict (FIXME when MPC will be upgraded)
  virtual bool Advance(double dt, double& dt_actual) = 0;
  // virtual bool Advance(double dt) = 0;

  // Commit any secondary (dependent) variables.
  virtual void CommitState(double dt, const Teuchos::RCP<State>& S) = 0;

  // Calculate any diagnostics prior to doing vis
  virtual void CalculateDiagnostics(const Teuchos::RCP<State>& S) = 0;

  // -- set pointers to State, and point the solution vector to the data in S_next
  virtual void SetStates(const Teuchos::RCP<const State>& S,
                         const Teuchos::RCP<State>& S_inter,
                         const Teuchos::RCP<State>& S_next) = 0;

  virtual std::string name() = 0;
};

}  // namespace Amanzi

#endif
