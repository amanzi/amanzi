/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Default base with default implementations of methods for a PK integrated using
BDF.
------------------------------------------------------------------------- */

#ifndef AMANZI_PK_BDF_BASE_HH_
#define AMANZI_PK_BDF_BASE_HH_

#include "bdf_fn_base.hh"
#include "bdf_time_integrator.hh"
#include "pk_default_base.hh"

namespace Amanzi {

class PKBDFBase : public virtual PKDefaultBase, public BDFFnBase {

 public:

  PKBDFBase(Teuchos::ParameterList& plist,
            const Teuchos::RCP<TreeVector>& solution) :
      PKDefaultBase(plist,solution) {}

  // Default implementations of PK methods.
  // -- setup
  virtual void setup(const Teuchos::Ptr<State>& S);

  // -- initialize
  virtual void initialize(const Teuchos::Ptr<State>& S);

  // -- Choose a time step compatible with physics.
  virtual double get_dt();

  // -- Advance from state S0 to state S1 at time S0.time + dt.
  virtual bool advance(double dt);

  // -- Check the admissibility of a solution.
  virtual bool is_admissible(Teuchos::RCP<const TreeVector> up);

  // -- Possibly modify the predictor that is going to be used as a
  //    starting value for the nonlinear solve in the time integrator.
  virtual bool modify_predictor(double h, Teuchos::RCP<TreeVector> up);

 protected: // data
  // timestep control
  double dt_;
  double dt0_;
  double time_step_reduction_factor_;
  Teuchos::RCP<BDFTimeIntegrator> time_stepper_;

  // error criteria
  double atol_, rtol_;

};

} // namespace

#endif
