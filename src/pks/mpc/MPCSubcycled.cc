/* -------------------------------------------------------------------------
Amanzi

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Class for subcycling a slave step within a master step.
Assumes that intermediate_time() can be used (i.e. this is not nestable?)

See additional documentation in the base class src/pks/mpc/MPC.hh
------------------------------------------------------------------------- */

#include "MPCSubcycled.hh"

namespace Amanzi {


// -----------------------------------------------------------------------------
// Constructor
// -----------------------------------------------------------------------------
MPCSubcycled::MPCSubcycled(const Teuchos::RCP<Teuchos::ParameterList>& plist,
             Teuchos::ParameterList& FElist,
             const Teuchos::RCP<TreeVector>& soln) :
    MPCTmp<PK>(plist, FElist, soln) {

  // Master PK is the PK whose time step size sets the size, the slave is subcycled.
  master_ = plist->get<int>("master PK index", 0);
  slave_ = master_ == 1 ? 0 : 1;

  if (sub_pks_.size() != 2 || master_ > 1) {
    Errors::Message message("MPCSubcycled: only MPCs with two sub-PKs can currently be subcycled.");
    Exceptions::amanzi_throw(message);
  }

  // min dt allowed in subcycling
  min_dt_ = plist->get<double>("mininum subcycled relative dt", 1.e-5);

}
  

// -----------------------------------------------------------------------------
// Calculate the min of sub PKs timestep sizes.
// -----------------------------------------------------------------------------
double MPCSubcycled::get_dt() {
  master_dt_ = sub_pks_[master_]->get_dt();
  slave_dt_ = sub_pks_[slave_]->get_dt();
  if (slave_dt_ > master_dt_) slave_dt_ = master_dt_;
  
  return master_dt_;
};


// -----------------------------------------------------------------------------
// Advance each sub-PK individually, returning a failure as soon as possible.
// -----------------------------------------------------------------------------
bool MPCSubcycled::Advance(double dt) {
  bool fail = false;

  // advance the master PK
  fail = sub_pks_[master_]->Advance(dt);
  if (fail) return fail;

  // --etc: unclear if state should be commited?
  sub_pks_[master_]->CommitState(dt, S_.ptr());

  // advance the slave, subcycling if needed
  S_->set_intermediate_time(S_->last_time());
  bool done = false;

  double dt_next = slave_dt_;
  double dt_done = 0.;
  while (!done) {
    // do not overstep
    if (dt_done + dt_next > dt) {
      dt_next = dt - dt_done;
    }

    // set the intermediate time
    S_->set_intermediate_time(S_->last_time() + dt_done + dt_next);

    // take the step
    fail = sub_pks_[slave_]->Advance(dt_next);

    if (fail) {
      // if fail, cut the step and try again
      dt_next /= 2;
    } else {
      // if success, commit the state and increment to next intermediate
      // -- etc: unclear if state should be commited or not?
      sub_pks_[slave_]->CommitState(dt_next, S_.ptr());
      dt_done += dt_next;
    }

    // check for done condition
    done = (std::abs(dt_done - dt) / dt < 0.1*min_dt_) || // finished the step
      (dt_next / dt < min_dt_); // failed
  }

  if (std::abs(dt_done - dt) / dt < 0.1*min_dt_) {
    // done, success
    // --etc: unclear if state should be commited or not?
    CommitState(dt, S_.ptr());
    return true;
  } else {
    return false;
  }
}


} // namespace
