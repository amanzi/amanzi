/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Interface for the Base MPC class.  A multi process coordinator is a PK
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

#ifndef PKS_MPC_MPC_HH_
#define PKS_MPC_MPC_HH_

#include <vector>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Epetra_MpiComm.h"

#include "state.hh"
#include "tree_vector.hh"

#include "PK.hh"
#include "pk_factory.hh"

namespace Amanzi {

class MPC : public PK {

public:
  MPC(Teuchos::ParameterList& mpc_plist, const Teuchos::RCP<State>& S,
        const Teuchos::RCP<TreeVector>& soln) :
      mpc_plist_(mpc_plist),
      n_iter_(0) {};

  // PK methods
  // -- calls all sub-PK initialize() methods
  virtual void initialize(const Teuchos::RCP<State>& S);

  // transfer operators
  virtual void state_to_solution(const Teuchos::RCP<State>& S,
          const Teuchos::RCP<TreeVector>& soln);
  virtual void solution_to_state(const Teuchos::RCP<TreeVector>& soln,
          const Teuchos::RCP<State>& S);

  // -- min(get_dt()) for each sub-PK
  virtual double get_dt();

  // -- loops over sub-PKs
  virtual bool advance(double dt) = 0;
  virtual void commit_state(double dt, const Teuchos::RCP<State>& S);
  virtual void calculate_diagnostics(const Teuchos::RCP<State>& S);

  // set States
  virtual void set_states(const Teuchos::RCP<const State>& S,
                          const Teuchos::RCP<State>& S_inter,
                          const Teuchos::RCP<State>& S_next);

protected:
  std::vector< Teuchos::RCP<PK> > sub_pks_;
  Teuchos::ParameterList mpc_plist_;
  int n_iter_;
};
} // close namespace Amanzi

#endif
