#ifndef _MPC_HH_
#define _MPC_HH_

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

class MPC : public Teuchos::VerboseObject<MPC>, public PK {

public:
  MPC(Teuchos::ParameterList& mpc_plist,
      Teuchos::RCP<State>& S, Teuchos::RCP<TreeVector>& soln);

  ~MPC() {};

  // PK methods
  virtual void initialize(Teuchos::RCP<State>& S,
                          Teuchos::RCP<TreeVector>& soln);

  // transfer operators
  virtual void state_to_solution(const State& S,
                                 Teuchos::RCP<TreeVector>& soln);
  virtual void state_to_solution(const State& S,
                                 Teuchos::RCP<TreeVector>& soln,
                                 Teuchos::RCP<TreeVector>& soln_dot);
  virtual void solution_to_state(const TreeVector& soln,
                                 Teuchos::RCP<State>& S);
  virtual void solution_to_state(const TreeVector& soln,
                                 const TreeVector& soln_dot,
                                 Teuchos::RCP<State>& S);

  virtual double get_dT();

  virtual bool advance(double dt, Teuchos::RCP<const State>& S0,
             Teuchos::RCP<State>& S1, Teuchos::RCP<TreeVector>& solution) = 0;

  virtual void commit_state(double dt, Teuchos::RCP<State>& S);

  virtual void calculate_diagnostics(Teuchos::RCP<State>& S);

protected:
  // PK container and factory
  PK_Factory pk_factory_;
  std::vector< Teuchos::RCP<PK> > sub_pks_;

  // states
  Teuchos::RCP<State> S_;

  // misc setup information
  Teuchos::ParameterList mpc_plist_;
};

} // close namespace Amanzi

#endif
