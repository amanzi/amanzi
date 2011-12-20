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
  // PK(Teuchos::ParameterList& plist, Teuchos::RCP<State>& S, Teuchos::RCP<TreeVector>& soln);

  //virtual ~PK() = 0;

  // Initialize owned (dependent) variables.
  virtual void initialize(Teuchos::RCP<State>& S, Teuchos::RCP<TreeVector>& soln) = 0;

  // transfer operators
  virtual void state_to_solution(Teuchos::RCP<const State>& S,
                                 Teuchos::RCP<TreeVector>& soln) = 0;
  virtual void state_to_solution(Teuchos::RCP<const State>& S,
                                 Teuchos::RCP<TreeVector>& soln,
                                 Teuchos::RCP<TreeVector>& soln_dot) = 0;
  virtual void solution_to_state(Teuchos::RCP<const TreeVector>& soln,
                                 Teuchos::RCP<State>& S) = 0;
  virtual void solution_to_state(Teuchos::RCP<const TreeVector>& soln,
                                 Teuchos::RCP<const TreeVector>& soln_dot,
                                 Teuchos::RCP<State>& S) = 0;

  // Choose a time step compatible with physics.
  virtual double get_dT() = 0;

  // Advance from state S0 to state S1 at time S0.time + dt.
  virtual bool advance(double dt, Teuchos::RCP<const State>& S0,
          Teuchos::RCP<State>& S1, Teuchos::RCP<TreeVector>& solution) = 0;

  // compute the residual function, for use as a sub-pk 
  // NOTE: compute_f gets called by the BDF/integration code, which means that
  // it must accept vectors and dynamic_cast them to TreeVectors.
  virtual void compute_f(const double t, const Vector& soln,
                         const Vector& soln_dot, Vector& f) = 0;

  // Commit any secondary (dependent) variables.
  virtual void commit_state(double dt, Teuchos::RCP<State>& S) = 0;

  // Calculate any diagnostics prior to doing vis
  virtual void calculate_diagnostics(Teuchos::RCP<State>& S) = 0;

  // get and set name
  std::string name() { return name_; }
  void set_name(std::string name) { name_ = name; }

private: 
  std::string name_;
  
};
} // namespace

#endif
