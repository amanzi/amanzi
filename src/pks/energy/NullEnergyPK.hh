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
  NullEnergyPK(Teuchos::ParameterList&, Teuchos::RCP<State> S);
  ~PK() {}

  // Initialize owned (dependent) variables.
  void initialize();

  // Choose a time step compatible with physics.
  double get_dT() { return 1.e99; }

  // Advance from state S0 to state S1 at time S0.time + dt.
  bool advance(double dt, const Teuchos::RCP<State> &S0,
          Teuchos::RCP<State> &S1, Teuchos::RCP<TreeVector> &solution);

  void compute_f(const double t, const Vector& u,
                         const Vector& udot, Vector& f);

  // Take solution for u,udot and put the various components into
  // the state.
  void solution_to_state(const TreeVector& u, const TreeVector& udot);
  void solution_to_state(const TreeVector& u, Teuchos::RCP<State> &S);

  // Commit any secondary (dependent) variables.
  virtual void commit_state(double dt, Teuchos::RCP<State> &S) {}

  // get and set name
  std::string name() { return name_; }
  void set_name(std::string) { name_ = name; }

private: 
  std::string name_;

  // states
  Teuchos::RCP<State> S_;
  double T_;

  // misc setup information
  Teuchos::ParameterList energy_plist_;
  
};
} // namespace

#endif
