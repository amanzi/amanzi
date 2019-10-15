/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon (coonet@ornl.gov)
*/

//! A mixin class with default implementations of methods for a subcycled BDF
//! integrated PK.

/*!

``PK_MixinImplicitSubcycled`` is a mixin class providing functionality for
implicitly integrated PKs which cannot fail over the requested timestep.  The
implicit step is subcycled until the entire integrated step is met.

Note this inherits everything in ``PK_MixinImplicit_``

* `"maximum time step`" ``[double]`` The maximum step size of the outer step.

*/

#ifndef AMANZI_PK_MIXIN_IMPLICIT_SUBCYCLED_HH_
#define AMANZI_PK_MIXIN_IMPLICIT_SUBCYCLED_HH_

#include "Teuchos_ParameterList.hpp"

#include "Key.hh"
#include "PK.hh"
#include "PK_MixinImplicit.hh"
#include "TimeStepManager.hh"

namespace Amanzi {

template <class Base_t>
class PK_MixinImplicitSubcycled : public PK_MixinImplicit<Base_t> {
 public:
  PK_MixinImplicitSubcycled(
    const Teuchos::RCP<Teuchos::ParameterList>& pk_tree,
    const Teuchos::RCP<Teuchos::ParameterList>& global_plist,
    const Teuchos::RCP<State>& S);

  void Setup();
  bool AdvanceStep(const Key& tag_old, const Key& tag_new);
  void CommitStep(const Key& tag_old, const Key& tag_new);
  double get_dt() { return dt_max_; }

 protected:
  using PK_MixinImplicit<Base_t>::tag_old_;
  using PK_MixinImplicit<Base_t>::tag_new_;
  using Base_t::S_;
  using Base_t::plist_;
  using Base_t::vo_;

  double dt_max_;
};

template <class Base_t>
PK_MixinImplicitSubcycled<Base_t>::PK_MixinImplicitSubcycled(
  const Teuchos::RCP<Teuchos::ParameterList>& pk_tree,
  const Teuchos::RCP<Teuchos::ParameterList>& global_plist,
  const Teuchos::RCP<State>& S)
  : PK_MixinImplicit<Base_t>(pk_tree, global_plist, S)
{
  // initial timestep
  dt_max_ = plist_->template get<double>("maximum time step", 1.e10);
};

template <class Base_t>
void
PK_MixinImplicitSubcycled<Base_t>::Setup()
{
  PK_MixinImplicit<Base_t>::Setup();

  // reserve space for inner step
  tag_old_ = this->name() + " implicit inner";
  S_->template Require<double>("time", tag_old_, "time");
  S_->template Require<int>("cycle", tag_old_, "cycle");
  this->SolutionToState(tag_old_, "");

  tag_new_ = this->name() + " implicit inner next";
  S_->template Require<double>("time", tag_new_, "time");
  S_->template Require<int>("cycle", tag_new_, "cycle");
  this->SolutionToState(tag_new_, "");
}

template <class Base_t>
bool
PK_MixinImplicitSubcycled<Base_t>::AdvanceStep(const Key& tag_old,
                                               const Key& tag_new)
{
  // times associated with inner and outer steps
  double t_final = S_->time(tag_new);
  double t_inner = S_->time(tag_old);
  double dt = t_final - t_inner;
  double dt_inner = PK_MixinImplicit<Base_t>::get_dt();

  // logging
  Teuchos::OSTab out = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_HIGH))
    *vo_->os()
      << "----------------------------------------------------------------"
      << std::endl
      << "Advancing: t0 = " << t_inner << " t1 = " << t_final << " h = " << dt
      << std::endl
      << "----------------------------------------------------------------"
      << std::endl;

  // set up the inner state -- times, cycle, and values
  S_->set_time(tag_old_, t_inner);
  S_->set_time(tag_new_, t_inner);
  S_->set_cycle(tag_old_, 0);
  S_->set_cycle(tag_new_, 0);
  this->StateToState(tag_old, tag_old_);
  this->StateToState(tag_old, tag_new_);
  //  this->ChangedSolutionPK(); Not needed? --etc

  TimeStepManager tsm;
  tsm.RegisterTimeEvent(t_final);

  // advance the outer
  bool done = false;
  bool fail = false;
  while (!done) {
    // -- get a valid dt
    double dt_ok = tsm.TimeStep(t_inner, dt_inner, fail);

    // -- advance the inner state
    S_->set_time(tag_new_, t_inner + dt_ok);
    S_->set_cycle(tag_new_, S_->cycle(tag_old_) + 1);

    // -- take the inner step
    fail = PK_MixinImplicit<Base_t>::AdvanceStep(tag_old_, tag_new_);
    if (!fail) {
      this->CommitStep(tag_old_, tag_new_);
      S_->set_time(tag_old_, S_->time(tag_new_));
      S_->set_cycle(tag_old_, S_->cycle(tag_new_));
      t_inner += dt_ok;
      done = t_inner >= t_final - 1.e-10;
    } else {
      this->FailStep(tag_old_, tag_new_);
    }
    dt_inner = PK_MixinImplicit<Base_t>::get_dt();
  }

  // copy from inner result to outer result
  this->StateToState(tag_new_, tag_new);

  // subcycling is always successful to the outer step
  return false;
}

// -- Commit any secondary (dependent) variables.
template <class Base_t>
void
PK_MixinImplicitSubcycled<Base_t>::CommitStep(const Key& tag_old,
                                              const Key& tag_new)
{
  // only update the timestepper if this is an internal call
  if (tag_new == tag_new_) {
    PK_MixinImplicit<Base_t>::CommitStep(tag_old, tag_new);
  } else {
    // otherwise sidestep PK_MixinImplicit
    Base_t::CommitStep(tag_old, tag_new);
  }
}

} // namespace Amanzi

#endif
