/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*
  MPC PKs

  Class for subcycling a slave step within a master step.
  Assumes that intermediate_time() can be used (i.e. this is not nestable?)

  NOTe: if master passes, slave will try to subcycle with adaptive time step
  or simulation fails.

  See additional documentation in the base class src/pks/mpc_pk/PK_MPC.hh
*/

#include "PK_MPCSubcycled.hh"

namespace Amanzi {

// -----------------------------------------------------------------------------
// Constructor
// -----------------------------------------------------------------------------
PK_MPCSubcycled::PK_MPCSubcycled(Teuchos::ParameterList& pk_tree,
                                 const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                                 const Teuchos::RCP<State>& S,
                                 const Teuchos::RCP<TreeVector>& soln)
  : PK_MPC<PK>(pk_tree, global_list, S, soln)
{
  // Master PK is the PK whose timestep size sets the size, the slave is subcycled.
  master_ = my_list_->get<int>("master PK index", 0);
  slave_ = (master_ == 1) ? 0 : 1;

  if (sub_pks_.size() != 2) {
    Errors::Message message(
      "PK_MPCSubcycled: only MPCs with two sub-PKs can currently be subcycled.");
    Exceptions::amanzi_throw(message);
  }

  // min dt allowed in subcycling
  min_dt_ = my_list_->get<double>("minimum subcycled relative dt", 1.e-5);
}


// -----------------------------------------------------------------------------
// Calculate the min of sub PKs timestep sizes.
// -----------------------------------------------------------------------------
double
PK_MPCSubcycled::get_dt()
{
  master_dt_ = sub_pks_[master_]->get_dt();
  slave_dt_ = sub_pks_[slave_]->get_dt();
  if (slave_dt_ > master_dt_) slave_dt_ = master_dt_;

  return master_dt_;
}


// -----------------------------------------------------------------------------
// Advance each sub-PK individually, returning a failure as soon as possible.
// -----------------------------------------------------------------------------
bool
PK_MPCSubcycled::AdvanceStep(double t_old, double t_new, bool reinit)
{
  bool fail = false;

  // advance the master PK using the full step size
  fail = sub_pks_[master_]->AdvanceStep(t_old, t_new, reinit);
  if (fail) return fail;

  master_dt_ = t_new - t_old;
  sub_pks_[master_]->CommitStep(t_old, t_new, Tags::DEFAULT);

  slave_dt_ = sub_pks_[slave_]->get_dt();
  if (slave_dt_ > master_dt_) slave_dt_ = master_dt_;

  // advance the slave, subcycling if needed
  S_->set_intermediate_time(t_old);
  bool done = false;

  double dt_next = slave_dt_;
  double dt_done = 0.;
  while (!done) {
    // do not overstep
    if (t_old + dt_done + dt_next > t_new) {
      dt_next = t_new - t_old - dt_done;
    }

    // take the step
    fail = sub_pks_[slave_]->AdvanceStep(t_old + dt_done, t_old + dt_done + dt_next, reinit);

    if (fail) {
      // if fail, cut the step and try again
      dt_next /= 2;
    } else {
      // set the intermediate time
      S_->set_intermediate_time(t_old + dt_done + dt_next);

      // -- etc: unclear if state should be commited or not?
      sub_pks_[slave_]->CommitStep(t_old + dt_done, t_old + dt_done + dt_next, Tags::DEFAULT);
      dt_done += dt_next;

      // allow dt to grow only when success
      dt_next = sub_pks_[slave_]->get_dt();
    }

    // check for done condition
    done =
      (std::abs(t_old + dt_done - t_new) / (t_new - t_old) < 0.1 * min_dt_) || (dt_next < min_dt_);
  }

  return false;
}


// -----------------------------------------------------------------------------
// Make necessary operatios by the end of the timesteps.
// -----------------------------------------------------------------------------
void
PK_MPCSubcycled::CommitStep(double t_old, double t_new, const Tag& tag)
{
  sub_pks_[slave_]->CommitStep(t_old, t_new, tag);
}

} // namespace Amanzi
