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

#include "State.hh"
#include "Vector.hh"
#include "TreeVector.hh"

#include "PK.hh"
#include "PK_Factory.hh"

namespace Amanzi
{

class MPC : public PK {

public:
  // all PK constructors must look like this
  MPC(Teuchos::ParameterList& mpc_plist,
      Teuchos::RCP<State>& S, Teuchos::RCP<TreeVector>& solution);

  ~MPC() {};

  // PK methods
  // -- calls all sub-PK initialize() methods
  virtual void initialize(Teuchos::RCP<State>& S);

  // transfer operators
  virtual void state_to_solution(State& S, Teuchos::RCP<TreeVector>& soln);
  virtual void state_to_solution(State& S, Teuchos::RCP<TreeVector>& soln,
                                 Teuchos::RCP<TreeVector>& soln_dot);
  virtual void solution_to_state(TreeVector& soln, Teuchos::RCP<State>& S);
  virtual void solution_to_state(TreeVector& soln, TreeVector& soln_dot,
          Teuchos::RCP<State>& S);

  // -- min(get_dt()) for each sub-PK
  virtual double get_dt();

  // -- loops over sub-PKs
  virtual bool advance(double dt) = 0;
  virtual void commit_state(double dt, Teuchos::RCP<State>& S);
  virtual void calculate_diagnostics(Teuchos::RCP<State>& S);

  // set States
  virtual void set_states(Teuchos::RCP<const State>& S, Teuchos::RCP<State>& S_next);

protected:
  // PK container and factory
  PK_Factory pk_factory_;
  std::vector< Teuchos::RCP<PK> > sub_pks_;

  // misc setup information
  Teuchos::ParameterList mpc_plist_;
};
} // close namespace Amanzi

#endif
