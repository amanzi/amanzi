/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

#ifndef _PK_HH_
#define _PK_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_ParameterList.hpp"

#include "State.hh"
#include "Vector.hh"
#include "TreeVector.hh"

namespace Amanzi {

// purely virtual interface class for PKs
class PK : public Teuchos::VerboseObject<PK> {

public:
  // Constructor should populate state with independent and dependent variables,
  // and have a signature of:
  // PK(Teuchos::ParameterList&, Teuchos::RCP<State> S);

  //virtual ~PK() = 0;

  // Initialize owned (dependent) variables.
  virtual void initialize() = 0;

  // Choose a time step compatible with physics.
  virtual double get_dT() = 0;

  // Advance from state S0 to state S1 at time S0.time + dt.
  // NOTE: advance gets called from above, and therefore can 
  // know that the solution is in fact a TreeVector.
  virtual bool advance(double dt, const Teuchos::RCP<State> &S0,
          Teuchos::RCP<State> &S1, Teuchos::RCP<TreeVector> &solution) = 0;

  // compute the residual function, for use as a sub-pk
  // NOTE CONTINUED: On the other hand, compute_f gets called by
  // the BDF/integration code, which means that it must accept
  // vectors and dynamic_cast them to TreeVectors.
  virtual void compute_f(const double t, const Vector& u,
                         const Vector& udot, Vector& f) = 0;

  // Take solution for u,udot and put the various components into
  // the state.
  virtual void solution_to_state(const TreeVector& u,
                                 const TreeVector& udot) = 0;

  // Commit any secondary (dependent) variables.
  virtual void commit_state(double dt, Teuchos::RCP<State> &S) = 0;

  // get and set name
  std::string name() { return name_; }
  void set_name(std::string) { name_ = name; }

private: 
  std::string name_;
  
};
} // namespace

#endif
