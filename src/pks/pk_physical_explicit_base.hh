/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Standard base for most diffusion-dominated PKs, this combines both
domains/meshes of PKPhysicalBase and Explicit methods of PKExplicitBase.
------------------------------------------------------------------------- */

#ifndef ATS_PK_PHYSICAL_EXPLICIT_BASE_HH_
#define ATS_PK_PHYSICAL_EXPLICIT_BASE_HH_

#include "errors.hh"
#include "pk_default_base.hh"
#include "pk_explicit_base.hh"
#include "pk_physical_base.hh"

namespace Amanzi {

class PKPhysicalExplicitBase : public PKExplicitBase, public PKPhysicalBase {

 public:
  PKPhysicalExplicitBase(const Teuchos::RCP<Teuchos::ParameterList>& plist,
                    Teuchos::ParameterList& FElist,
                    const Teuchos::RCP<TreeVector>& solution) :
      PKDefaultBase(plist, FElist, solution),
      PKPhysicalBase(plist, FElist, solution),
      PKExplicitBase(plist, FElist, solution) {}

  virtual void setup(const Teuchos::Ptr<State>& S) {
    PKPhysicalBase::setup(S);
    PKExplicitBase::setup(S);
  }

  // initialize.  Note both ExplicitBase and PhysicalBase have initialize()
  // methods, so we need a unique overrider.
  virtual void initialize(const Teuchos::Ptr<State>& S) {
    PKPhysicalBase::initialize(S);
    PKExplicitBase::initialize(S);
  }

  // -- Advance from state S0 to state S1 at time S0.time + dt.
  virtual bool advance(double dt) {
    PKExplicitBase::advance(dt);
    ChangedSolution();
    return false;
  }


  // -- Experimental approach -- calling this indicates that the time
  //    integration scheme is changing the value of the solution in
  //    state.
  virtual void ChangedSolution() {
    solution_evaluator_->SetFieldAsChanged(S_next_.ptr());
  }

};


} // namespace

#endif
