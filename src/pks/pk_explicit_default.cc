/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Default base with default implementations of methods for a PK integrated using
Explicit.
------------------------------------------------------------------------- */

#include "Teuchos_TimeMonitor.hpp"
#include "PK.hh"
#include "State.hh"
#include "boost/algorithm/string.hpp"
#include "pk_explicit_default.hh"

namespace Amanzi {
  
// -----------------------------------------------------------------------------
// Setup
// -----------------------------------------------------------------------------
void PK_Explicit_Default::Setup(const Teuchos::Ptr<State>& S) {


  // initial timestep
  dt_ = plist_->get<double>("initial time step", 1.);

};


// -----------------------------------------------------------------------------
// Initialization of timestepper.
// -----------------------------------------------------------------------------
void PK_Explicit_Default::Initialize(const Teuchos::Ptr<State>& S) {
  // set up the timestepping algorithm
  if (!plist_->get<bool>("strongly coupled PK", false)) {
    // -- instantiate time stepper
    Teuchos::ParameterList& ti_plist = plist_->sublist("time integrator");
    ti_plist.set("initial time", S->time());
    time_stepper_ = Teuchos::rcp(new Explicit_TI::RK<TreeVector>(*this, ti_plist, *solution_));

    solution_old_ = Teuchos::rcp(new TreeVector(*solution_));
  }

};


// -----------------------------------------------------------------------------
// Initialization of timestepper.
// -----------------------------------------------------------------------------
double PK_Explicit_Default::get_dt() { return dt_; }

void PK_Explicit_Default::set_dt(double dt) { dt_ = dt; }


// -----------------------------------------------------------------------------
// Advance from state S to state S_next at time S.time + dt.
// -----------------------------------------------------------------------------
bool PK_Explicit_Default::AdvanceStep(double t_old, double t_new, bool reinit) {

  double dt = t_new - t_old;  
  
  Teuchos::OSTab out = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_HIGH))
    *vo_->os() << "----------------------------------------------------------------" << std::endl
               << "Advancing: t0 = " << S_inter_->time()
               << " t1 = " << S_next_->time() << " h = " << dt << std::endl
               << "----------------------------------------------------------------" << std::endl;

  State_to_Solution(S_inter_, *solution_old_);
  State_to_Solution(S_next_, *solution_);

  // take a timestep
  if (true) { // this is here simply to create a context for timer,
              // which stops the clock when it is destroyed at the
              // closing brace.
    time_stepper_->TimeStep(S_inter_->time(), dt, *solution_old_, *solution_);
  }
  return false;
};


} // namespace
