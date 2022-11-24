/*
  License: see $AMANZI_DIR/COPYRIGHT
  Authors: Daniil Svyatskiy

  PK for coupling of Flow_PK with Transport_PK and Chemestry_PK

*/


#include "PressureSaturation_PK.hh"

namespace Amanzi {

// -----------------------------------------------------------------------------
// Constructor
// -----------------------------------------------------------------------------
PressureSaturation_PK::PressureSaturation_PK(
  Teuchos::ParameterList& pk_tree,
  const Teuchos::RCP<Teuchos::ParameterList>& global_list,
  const Teuchos::RCP<State>& S,
  const Teuchos::RCP<TreeVector>& soln)
  : Amanzi::MPCSubcycled(pk_tree, global_list, S, soln)
{}


// -----------------------------------------------------------------------------
// Calculate the min of sub PKs timestep sizes.
// -----------------------------------------------------------------------------
double
PressureSaturation_PK::get_dt()
{
  double dt = Amanzi::MPCSubcycled::get_dt();
  return dt;
};

// -----------------------------------------------------------------------------
// Set master dt
// -----------------------------------------------------------------------------

void
PressureSaturation_PK::set_dt(double dt)
{
  master_dt_ = dt;
  slave_dt_ = sub_pks_[slave_]->get_dt();
  //  std::cout<<"master_dt_ "<<master_dt_<<" slave_dt_ "<<slave_dt_<<"\n";
  if (slave_dt_ > master_dt_) slave_dt_ = master_dt_;
}


// -----------------------------------------------------------------------------
// Advance each sub-PK individually, returning a failure as soon as possible.
// -----------------------------------------------------------------------------
bool
PressureSaturation_PK::AdvanceStep(double t_old, double t_new, bool reinit)
{
  //   bool fail = Amanzi::MPCSubcycled::AdvanceStep(t_old, t_new);
  bool fail = false;

  // advance the master PK using the full step size
  //std::cout << "Master PK: " << sub_pks_[master_]->name() << "\n";
  //std::cout<<"Advance Master PK: "<<t_old<<" "<<t_new<<"\n";
  fail = sub_pks_[master_]->AdvanceStep(t_old, t_new, reinit);
  if (fail) return fail;

  master_dt_ = t_new - t_old;
  if (slave_dt_ > master_dt_) slave_dt_ = master_dt_;

  // --etc: unclear if state should be commited?
  sub_pks_[master_]->CommitStep(t_old, t_new);

  // advance the slave, subcycling if needed
  S_->set_intermediate_time(t_old);
  bool done = false;

  double dt_next = slave_dt_;
  double dt_done = 0.;
  while (!done) {
    // do not overstep
    if (t_old + dt_done + dt_next > t_new) { dt_next = t_new - t_old - dt_done; }

    // set the intermediate time
    S_->set_intermediate_time(t_old + dt_done + dt_next);

    // take the step
    //std::cout << "Slave PK: " << sub_pks_[slave_]->name() << "\n";
    //std::cout<<"Advance Slave PK: "<<t_old<<" "<<t_new<<"\n";
    fail = sub_pks_[slave_]->AdvanceStep(t_old + dt_done, t_old + dt_done + dt_next, reinit);

    if (fail) {
      // if fail, cut the step and try again
      dt_next /= 2;
    } else {
      // if success, commit the state and increment to next intermediate
      // -- etc: unclear if state should be commited or not?
      sub_pks_[slave_]->CommitStep(t_old + dt_done, t_old + dt_done + dt_next);
      dt_done += dt_next;
    }

    //std::cout << "t_old: " << t_old << "; t_new: " << t_new << "; dt_done: " << dt_done << "; dt_next: " << dt_next << "\n";

    // check for done condition
    done =
      (std::abs(t_old + dt_done - t_new) / (t_new - t_old) < 0.1 * min_dt_) || // finished the step
      (dt_next < min_dt_);                                                     // failed
  }

  if (std::abs(t_old + dt_done - t_new) / (t_new - t_old) < 0.1 * min_dt_) {
    // done, success
    // --etc: unclear if state should be commited or not?
    CommitStep(t_old, t_new);
    return false;
  } else {
    return true;
  }
};


} // namespace Amanzi
