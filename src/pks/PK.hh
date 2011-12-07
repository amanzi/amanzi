/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

#ifndef _PK_HH_
#define _PK_HH_


#include "Teuchos_RCP.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_ParameterList.hpp"
#include "State.hh"

namespace Amanzi {

// purely virtual interface class for PKs
class PK : public Teuchos::VerboseObject<PK> {

public:
  // Constructor should populate state with independent and dependent variables,
  // and have a signature of:
  // PK(Teuchos::ParameterList&, Teuchos::RCP<State> S);

  //virtual ~PK() = 0;

  // Initialize owned (dependent) variables.
  virtual void initialize_state(Teuchos::RCP<State> S) = 0;

  // Choose a time step compatible with physics.
  virtual double get_dT() = 0;

  // Advance from state S0 to state S1 at time S0.time + dt.
  virtual bool advance_transient(double dt, const Teuchos::RCP<State> S0,
                                 Teuchos::RCP<State> S1) = 0;

  // Commit any secondary (dependent) variables.
  virtual void commit_state(double dt, Teuchos::RCP<State> S) = 0;
};

} // namespace

#endif
