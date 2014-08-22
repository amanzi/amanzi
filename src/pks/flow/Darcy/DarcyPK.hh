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
    void state_to_solution(Teuchos::RCP<State>& S, Teuchos::RCP<TreeVector>& soln);
    void state_to_solution(Teuchos::RCP<State>& S, Teuchos::RCP<TreeVector>& soln,
                           Teuchos::RCP<TreeVector>& soln_dot);
    void solution_to_state(Teuchos::RCP<TreeVector>& soln, Teuchos::RCP<State>& S);
    void solution_to_state(Teuchos::RCP<TreeVector>& soln,
                           Teuchos::RCP<TreeVector>& soln_dot,Teuchos::RCP<State>& S);

    // Steady state Darcy -- no timestep limitations
    double get_dt() {
      return 1.e99;
    }

    // Steady state Darcy -- transient does nothing.
    bool advance(double dt, Teuchos::RCP<TreeVector>& solution) {
      return 0;
    }

    // Transient does nothing, so commit not needed.
    void commit_state(double dt, Teuchos::RCP<State>& S) {}

    // Calculate flow on cells for diagnostic purposes.
    void calculate_diagnostics(Teuchos::RCP<State>& S);

  private:
    void calculate_velocity(Teuchos::RCP<Epetra_MultiVector>& q) const;
    void advance_to_steady_state(Teuchos::RCP<TreeVector>& soln);

    // TODO(etc): Make these scoped pointers instead of RCPs -- they should NOT get copied
    Teuchos::RCP<DarcyProblem> problem_;
    Teuchos::RCP<AztecOO> solver_;
    Teuchos::RCP<Matrix_MFD> matrix_;
    Teuchos::RCP<Matrix_MFD> preconditioner_;
    

    int num_itrs_sss_;      // count of linear solver iterations
    int max_itrs_sss_;      // max number of linear solver iterations
    double convergence_tol_sss_;
    double residual_sss_;

    Teuchos::ParameterList flow_plist_;
  };

} //namespace

#endif
