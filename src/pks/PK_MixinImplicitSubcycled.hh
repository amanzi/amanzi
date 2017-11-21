/* -*-  mode: c++; indent-tabs-mode: nil -*- */
//! A mixin class with default implementations of methods for an BDF integrated PK.

/*
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/


/*!

*/

#ifndef AMANZI_PK_MIXIN_IMPLICIT_SUBCYCLED_HH_
#define AMANZI_PK_MIXIN_IMPLICIT_SUBCYCLED_HH_

#include "Teuchos_ParameterList.hpp"

#include "Key.hh"
#include "PK_MixinImplicit.hh"
#include "PK.hh"
#include "TimeStepManager.hh"

namespace Amanzi {

template<class Base_t>
class PK_MixinImplicitSubcycled
    : public PK_MixinImplicit<Base_t> {
 public:
  PK_MixinImplicitSubcycled(const Teuchos::RCP<Teuchos::ParameterList>& pk_tree,
                            const Teuchos::RCP<Teuchos::ParameterList>& global_plist,
                            const Teuchos::RCP<State>& S,
                            const Teuchos::RCP<TreeVector>& solution) :
      PK_MixinImplicit<Base_t>(pk_tree, global_plist, S, solution) {}

  void Setup(const TreeVector& soln);
  bool AdvanceStep(const Key& tag_old, const Teuchos::RCP<TreeVector>& soln_old,
                   const Key& tag_new, const Teuchos::RCP<TreeVector>& soln_new);
  double get_dt() { return 1.e10; }

 protected:
  Teuchos::RCP<TreeVector> inner_soln_old_;
  Teuchos::RCP<TreeVector> inner_soln_new_;
  double inner_dt_;

  using PK_MixinImplicit<Base_t>::tag_old_;
  using PK_MixinImplicit<Base_t>::tag_new_;
  using Base_t::S_;
};


template<class Base_t>
void
PK_MixinImplicitSubcycled<Base_t>::Setup(const TreeVector& soln)
{
  PK_MixinImplicit<Base_t>::Setup(soln);

  // reserve space for inner step
  tag_old_ = this->name() + " implicit inner";
  S_->template Require<double>("time", tag_old_, "time");
  S_->template Require<int>("cycle", tag_old_, "cycle");
  inner_soln_old_ = Teuchos::rcp(new TreeVector(soln, INIT_MODE_NOALLOC));
  this->SolutionToState(*inner_soln_old_, tag_old_, "");

  tag_new_ = this->name() + " implicit inner next";
  S_->template Require<double>("time", tag_new_, "time");
  S_->template Require<int>("cycle", tag_new_, "cycle");
  inner_soln_new_ = Teuchos::rcp(new TreeVector(soln, INIT_MODE_NOALLOC));
  this->SolutionToState(*inner_soln_new_, tag_new_, "");
}

template<class Base_t>
bool 
PK_MixinImplicitSubcycled<Base_t>::AdvanceStep(const Key& tag_old,
        const Teuchos::RCP<TreeVector>& soln_old,
        const Key& tag_new,
        const Teuchos::RCP<TreeVector>& soln_new)
{
  bool fail = false;
  bool done = false;

  double t_final = S_->time(tag_new);
  double t_inner = S_->time(tag_old);
  double dt_inner = t_final - t_inner;

  TimeStepManager tsm;
  tsm.RegisterTimeEvent(t_final);
  
  this->CommitStep(tag_old_, inner_soln_old_, tag_old, soln_old);
  this->CommitStep(tag_new_, inner_soln_new_, tag_old, soln_old);
  S_->set_time(tag_old_, t_inner);
  S_->set_time(tag_new_, t_inner);
  S_->set_cycle(tag_old_, 0);
  S_->set_cycle(tag_new_, 0);
  
  //  this->ChangedSolutionPK(inner_soln_new_); ?? ETC
  
  while (!done) {
    double dt_ok = tsm.TimeStep(t_inner, dt_inner, fail);
    S_->set_time(tag_new_, t_inner + dt_ok);
    S_->set_cycle(tag_new_, S_->cycle(tag_old_)+1);
    fail = PK_MixinImplicit<Base_t>::AdvanceStep(tag_old_, inner_soln_old_,
            tag_new_, inner_soln_new_);
    if (!fail) {
      this->CommitStep(tag_old_, inner_soln_old_, tag_new_, inner_soln_new_);
      S_->set_time(tag_old_, S_->time(tag_new_));
      S_->set_cycle(tag_old_, S_->cycle(tag_new_));
      t_inner += dt_ok;
      done = t_inner >= t_final - 1.e-10;

    } else {
      this->FailStep(tag_old_, inner_soln_old_, tag_new_, inner_soln_new_);
    }
    dt_inner = PK_MixinImplicit<Base_t>::get_dt();
  }

  this->CommitStep(tag_new, soln_new, tag_new_, inner_soln_new_);
  return false;
}

} // namespace Amanzi

#endif
