/*
  This is the mpc_pk component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatskiy

  PK for coupling of Flow_PK with Transport_PK and Chemestry_PK
*/

#include "mpc_flowreactivetransport_pk.hh"

namespace Amanzi {

// -----------------------------------------------------------------------------
// Constructor
// -----------------------------------------------------------------------------
FlowReactiveTransport_PK_ATS::FlowReactiveTransport_PK_ATS(
    Teuchos::ParameterList& pk_tree,
    const Teuchos::RCP<Teuchos::ParameterList>& global_list,
    const Teuchos::RCP<State>& S,
    const Teuchos::RCP<TreeVector>& soln) :
    PK_MPCSubcycled_ATS(pk_tree, global_list, S, soln) { 

    Teuchos::ParameterList vlist;
    vlist.sublist("verbose object") = global_list -> sublist("verbose object");
    vo_ =  Teuchos::rcp(new VerboseObject("FlowandTransportPK", vlist)); 

}


// -----------------------------------------------------------------------------
// Calculate the min of sub PKs timestep sizes.
// -----------------------------------------------------------------------------
double FlowReactiveTransport_PK_ATS::get_dt() {
  double dt = Amanzi::PK_MPCSubcycled_ATS::get_dt();
  //double dt = sub_pks_[master_]->get_dt();
  return dt;
}


// -----------------------------------------------------------------------------
// Set master dt
// -----------------------------------------------------------------------------
void FlowReactiveTransport_PK_ATS::set_dt(double dt) {
  master_dt_ = dt;
  sub_pks_[master_]->set_dt(dt);
}


// -----------------------------------------------------------------------------
// Make necessary operatios by the end of the time steps.
// -----------------------------------------------------------------------------
void FlowReactiveTransport_PK_ATS::CommitStep(double t_old, double t_new, const Teuchos::RCP<State>& S) {
  sub_pks_[slave_]->CommitStep(t_old, t_new, S);
}


// -----------------------------------------------------------------------------
// Advance each sub-PK individually, returning a failure as soon as possible.
// -----------------------------------------------------------------------------
bool FlowReactiveTransport_PK_ATS::AdvanceStep(double t_old, double t_new, bool reinit) {
  bool fail = false;

  // advance the master PK using the full step size
  
  fail = sub_pks_[master_]->AdvanceStep(t_old, t_new, reinit);
  fail |= !sub_pks_[master_]->ValidStep();
  
  if (fail) return fail;

  //return fail;

  master_dt_ = t_new - t_old;

  sub_pks_[master_]->CommitStep(t_old, t_new, S_next_);

  //  S_next_->WriteStatistics(vo_);  

  slave_dt_ = sub_pks_[slave_]->get_dt();

  if (slave_dt_ > master_dt_) slave_dt_ = master_dt_;

  // advance the slave, subcycling if needed
  S_->set_intermediate_time(t_old);
  S_next_->set_intermediate_time(t_old);
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
      // if success, commit the state and increment to next intermediate
      // -- etc: unclear if state should be commited or not?
      // set the intermediate time
      S_ -> set_intermediate_time(t_old + dt_done + dt_next);
      S_next_ -> set_intermediate_time(t_old + dt_done + dt_next);
      //sub_pks_[slave_]->CommitStep(t_old + dt_done, t_old + dt_done + dt_next, S_);
      sub_pks_[slave_]->CommitStep(t_old + dt_done, t_old + dt_done + dt_next, S_next_);
      dt_done += dt_next;
    }

    // check for done condition
    done = (std::abs(t_old + dt_done - t_new) / (t_new - t_old) < 0.1*min_dt_) || // finished the step
        (dt_next  < min_dt_); // failed
  }


  if (std::abs(t_old + dt_done - t_new) / (t_new - t_old) < 0.1*min_dt_) {
    // done, success
    return false;
  } else {
    return true;
  }  
}

}  // namespace Amanzi

