/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Default base with default implementations of methods for a PK integrated using
Explicit.
------------------------------------------------------------------------- */

#ifndef AMANZI_PK_EXPLICIT_BASE_HH_
#define AMANZI_PK_EXPLICIT_BASE_HH_

#include "Teuchos_TimeMonitor.hpp"

#include "Explicit_TI_RK.hh"
#include "Explicit_TI_FnBase.hh"
#include "pk_default_base.hh"

namespace Amanzi {

class PKExplicitBase : public virtual PKDefaultBase,
                       public Explicit_TI::fnBase<TreeVector> {

 public:

  PKExplicitBase(const Teuchos::RCP<Teuchos::ParameterList>& plist,
            Teuchos::ParameterList& FElist,
            const Teuchos::RCP<TreeVector>& solution) :
      PKDefaultBase(plist, FElist, solution) {}

  // Virtual destructor
  virtual ~PKExplicitBase() {}

  // Default implementations of PK methods.
  // -- setup
  virtual void setup(const Teuchos::Ptr<State>& S);

  // -- initialize
  virtual void initialize(const Teuchos::Ptr<State>& S);

  // -- Choose a time step compatible with physics.
  virtual double get_dt();

  // -- Advance from state S0 to state S1 at time S0.time + dt.
  virtual bool advance(double dt);

 protected: // data
  // timestep control
  double dt_;
  Teuchos::RCP<Explicit_TI::RK<TreeVector> > time_stepper_;

  // timing
  Teuchos::RCP<Teuchos::Time> step_walltime_;

  // solution at the old timestep
  Teuchos::RCP<TreeVector> solution_old_;

};

} // namespace

#endif
