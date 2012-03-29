/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

simple passive tracer transport

-- designed to test the advection operators
------------------------------------------------------------------------- */

#ifndef PKS_TRANSPORT_PASSIVE_TRACER_HH_
#define PKS_TRANSPORT_PASSIVE_TRACER_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_ParameterList.hpp"

#include "composite_vector.hh"
#include "tree_vector.hh"
#include "state.hh"
#include "advection.hh"
#include "boundary-function.hh"

#include "PK.hh"
#include "bdf_fn_base.hh"
#include "bdf_time_integrator.hh"


namespace Amanzi {
namespace Transport {

class PassiveTracer : public PK, public BDFFnBase {

public:

  PassiveTracer(Teuchos::ParameterList& transport_plist, Teuchos::RCP<State>& S,
               Teuchos::RCP<TreeVector>& solution);

  // ConstantTemperature is a PK
  // -- Initialize owned (dependent) variables.
  virtual void initialize(const Teuchos::RCP<State>& S);

  // -- transfer operators -- pointer copy
  virtual void state_to_solution(const Teuchos::RCP<State>& S,
                                 const Teuchos::RCP<TreeVector>& soln);
  virtual void solution_to_state(const Teuchos::RCP<TreeVector>& soln,
                                 const Teuchos::RCP<State>& S);

  // -- Choose a time step compatible with physics.
  virtual double get_dt() { return dt_; }

  // -- Advance from state S to state S_next at time S0.time + dt.
  virtual bool advance(double dt);

  // -- Commit any secondary (dependent) variables.
  virtual void commit_state(double dt, const Teuchos::RCP<State>& S) {}

  // -- Update diagnostics for vis.
  virtual void calculate_diagnostics(const Teuchos::RCP<State>& S) {}

  // ConstantTemperature is a BDFFnBase
  // computes the non-linear functional f = f(t,u,udot)
  virtual void fun(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
                   Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> f);

  // applies preconditioner to u and returns the result in Pu
  virtual void precon(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu);

  // computes a norm on u-du and returns the result
  virtual double enorm(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<const TreeVector> du);

  // updates the preconditioner
  virtual void update_precon(double t, Teuchos::RCP<const TreeVector> up, double h) {}

private:
  void process_parameter_list(const Teuchos::RCP<State>& S);

  void AddAccumulation_(Teuchos::RCP<CompositeVector> f);
  void AddAdvection_(Teuchos::RCP<CompositeVector> f, bool negate=false);

  // states
  double C_;
  double dt_;

  // misc setup information
  Teuchos::ParameterList transport_plist_;

  // time integration
  Teuchos::RCP<Amanzi::BDFTimeIntegrator> time_stepper_;
  double atol_;
  double rtol_;

  // operators
  Teuchos::RCP<Operators::Advection> advection_;
  double cfl_;
  Teuchos::RCP< std::vector< Teuchos::RCP<BoundaryFunction> > > bcs_;
  Teuchos::RCP< std::vector<int> > bcs_dof_;
};

} // namespace
} // namespace

#endif
