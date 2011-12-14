/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

#ifndef __DARCY_PK_HH__
#define __DARCY_PK_HH__

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Epetra_Vector.h"
#include "AztecOO.h"

#include "PK.hh"
#include "State.hh"
#include "DarcyProblem.hh"

namespace Amanzi {

class Darcy_PK : public PK {

public:
  // Populate state with pressure, porosity, perm
  Darcy_PK(Teuchos::ParameterList &plist, Teuchos::RCP<State> &S);

  // Initialize will do all the work, running to steady state.
  void initialize(Teuchos::RCP<State> &S);

  // Steady state Darcy -- no timestep limitations
  double get_dT() { return 1.e99; }

  // Steady state Darcy -- transient does nothing.
  bool advance_transient(double dt, const Teuchos::RCP<State> &S0,
                         Teuchos::RCP<State> &S1) { return 0; }

  // Transient does nothing, so commit not needed.
  void commit_state(double dt, Teuchos::RCP<State> &S) {}

  // Access to solution (are these necessary?)
  const Epetra_MultiVector& get_pressure() const { return *pressure; }
  const Epetra_MultiVector& get_darcy_flux() const { return *darcy_flux; }

  // Compute the components of the Darcy velocity on cells (can this be private?)
  void calculate_velocity(Epetra_MultiVector &q) const {
    problem->DeriveDarcyVelocity(*solution, q); }

private:
  void advance_to_steady_state();

  Teuchos::RCP<DarcyProblem> problem;
  Teuchos::RCP<AztecOO> solver;

  Teuchos::RCP<Epetra_Vector> solution;   // full cell/face solution
  Teuchos::RCP<Epetra_Vector> pressure;   // cell pressures
  Teuchos::RCP<Epetra_Vector> darcy_flux; // Darcy face fluxes

  int max_itr;      // max number of linear solver iterations
  double err_tol;   // linear solver convergence error tolerance

  Teuchos::ParameterList plist;
};

} //namespace

#endif
