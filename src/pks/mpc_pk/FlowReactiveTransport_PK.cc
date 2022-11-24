/*
  This is the mpc_pk component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatskiy

  PK for coupling of Flow_PK with Transport_PK and Chemistry_PK
*/

#include "FlowReactiveTransport_PK.hh"

namespace Amanzi {

// -----------------------------------------------------------------------------
// Constructor
// -----------------------------------------------------------------------------
FlowReactiveTransport_PK::FlowReactiveTransport_PK(
  Teuchos::ParameterList& pk_tree,
  const Teuchos::RCP<Teuchos::ParameterList>& global_list,
  const Teuchos::RCP<State>& S,
  const Teuchos::RCP<TreeVector>& soln)
  : Amanzi::PK_MPCSubcycled(pk_tree, global_list, S, soln)
{}


// -----------------------------------------------------------------------------
// Calculate the min of sub PKs timestep sizes.
// -----------------------------------------------------------------------------
double
FlowReactiveTransport_PK::get_dt()
{
  //double dt = Amanzi::PK_MPCSubcycled::get_dt();
  double dt = sub_pks_[master_]->get_dt();
  return dt;
}


// -----------------------------------------------------------------------------
// Set master dt
// -----------------------------------------------------------------------------
void
FlowReactiveTransport_PK::set_dt(double dt)
{
  master_dt_ = dt;
  sub_pks_[master_]->set_dt(dt);
}


// -----------------------------------------------------------------------------
// Make necessary operatios by the end of the time steps.
// -----------------------------------------------------------------------------
void
FlowReactiveTransport_PK::CommitStep(double t_old, double t_new, const Tag& tag)
{
  sub_pks_[slave_]->CommitStep(t_old, t_new, tag);
}


// -----------------------------------------------------------------------------
// Advance each sub-PK individually, returning a failure as soon as possible.
// -----------------------------------------------------------------------------
bool
FlowReactiveTransport_PK::AdvanceStep(double t_old, double t_new, bool reinit)
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
    if (t_old + dt_done + dt_next > t_new) { dt_next = t_new - t_old - dt_done; }

    // take the step
    fail = sub_pks_[slave_]->AdvanceStep(t_old + dt_done, t_old + dt_done + dt_next, reinit);

    if (fail) {
      // if fail, cut the step and try again
      dt_next /= 2;
    } else {
      // if success, commit the state and increment to next intermediate
      // -- etc: unclear if state should be commited or not?
      // set the intermediate time
      S_->set_intermediate_time(t_old + dt_done + dt_next);
      sub_pks_[slave_]->CommitStep(t_old + dt_done, t_old + dt_done + dt_next, Tags::DEFAULT);
      dt_done += dt_next;
      // allow dt to grow only when success
      dt_next = sub_pks_[slave_]->get_dt();
    }

    // dt_next = sub_pks_[slave_]->get_dt();
    // no state recovery (e.g. pressure) is made, so the only option is to fail.
    if (dt_next < min_dt_)
      Exceptions::amanzi_throw("Failure in ReactiveTransport_PK: small time step.");

    // check for subcycling condition
    done = std::abs(t_old + dt_done - t_new) / (t_new - t_old) < 0.1 * min_dt_;
  }

  // we reach this point when subcycling has been completed
  return false;
}

} // namespace Amanzi
