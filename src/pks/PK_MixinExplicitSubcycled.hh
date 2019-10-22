/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon (coonet@ornl.gov)
*/

//! A mixin class with default implementations of methods for a subcycled
//! explicit time integrator.

/*!

``PK_MixinExplicitSubcycled`` is a mixin class providing functionality for
explicitly integrated PKs.  Manages the creation of intermediate data and
AdvanceStep().  This uses an extremely simple strategy of taking 1/Nth of the
requested timestep, where N is provided by the user.

TODO: This isn't very useful, as the timestep should likely be set by a CFL or
some other quantity of interest in a physical way.  Might think of ways of
extending the (protected) interface to allow this.

Note this inherits everything in ``PK_MixinExplicit_``

* `"subcycling substeps per outer step`" ``[int]``

  N above.  Number of steps per outer step requested.

*/

#ifndef AMANZI_PK_MIXIN_EXPLICIT_SUBCYCLED_HH_
#define AMANZI_PK_MIXIN_EXPLICIT_SUBCYCLED_HH_

#include "Teuchos_ParameterList.hpp"

#include "Key.hh"
#include "PK.hh"
#include "PK_MixinExplicit.hh"

namespace Amanzi {

template <class Base_t>
class PK_MixinExplicitSubcycled : public PK_MixinExplicit<Base_t> {
 public:
  PK_MixinExplicitSubcycled(
    const Teuchos::RCP<Teuchos::ParameterList>& pk_tree,
    const Teuchos::RCP<Teuchos::ParameterList>& global_plist,
    const Teuchos::RCP<State>& S);

  void Setup();
  bool AdvanceStep(const Key& tag_old, const Key& tag_new);

 protected:
  using PK_MixinExplicit<Base_t>::tag_old_;
  using PK_MixinExplicit<Base_t>::tag_new_;
  using PK_MixinExplicit<Base_t>::tag_inter_;

  using PK_MixinExplicit<Base_t>::vo_;
  using PK_MixinExplicit<Base_t>::plist_;
  using PK_MixinExplicit<Base_t>::S_;
  int subcycled_count_;
};

template <class Base_t>
PK_MixinExplicitSubcycled<Base_t>::PK_MixinExplicitSubcycled(
  const Teuchos::RCP<Teuchos::ParameterList>& pk_tree,
  const Teuchos::RCP<Teuchos::ParameterList>& global_plist,
  const Teuchos::RCP<State>& S)
  : PK_MixinExplicit<Base_t>(pk_tree, global_plist, S)
{
  // this could be generalized, for now just take 1/Nth step size
  subcycled_count_ =
    plist_->template get<int>("subcycling substeps per outer step");
}

template <class Base_t>
void
PK_MixinExplicitSubcycled<Base_t>::Setup()
{
  PK_MixinExplicit<Base_t>::Setup();

  // reserve space for inner step
  tag_old_ = this->name() + " explicit subcycled inner";
  S_->template Require<double>("time", tag_old_, "time");
  S_->template Require<int>("cycle", tag_old_, "cycle");
  this->SolutionToState(tag_old_, "");

  tag_new_ = this->name() + " explicit subcycled inner next";
  S_->template Require<double>("time", tag_new_, "time");
  S_->template Require<int>("cycle", tag_new_, "cycle");
  this->SolutionToState(tag_new_, "");
}

template <class Base_t>
bool
PK_MixinExplicitSubcycled<Base_t>::AdvanceStep(const Key& tag_old,
                                               const Key& tag_new)
{
  // times associated with inner and outer steps
  double t_start = S_->time(tag_old);
  double t_final = S_->time(tag_new);
  double dt_inner = (t_final - t_start) / subcycled_count_;

  // logging
  Teuchos::OSTab out = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_HIGH))
    *vo_->os()
      << "----------------------------------------------------------------"
      << std::endl
      << "Subcycling, Outer step: t0 = " << t_start << " t1 = " << t_final
      << " h = " << t_final - t_start << std::endl
      << "----------------------------------------------------------------"
      << std::endl;

  // copy initial condition to inner
  this->StateToState(tag_old, tag_old_);

  // create solution vectors, old and intermediate, which are pointers into
  // state data
  auto soln_space = this->SolutionSpace();
  auto soln_old = Teuchos::rcp(new TreeVector(soln_space, InitMode::NOALLOC));
  this->StateToSolution(*soln_old, tag_old_, "");

  auto soln_inter = Teuchos::rcp(new TreeVector(soln_space, InitMode::NOALLOC));
  this->StateToSolution(*soln_inter, tag_inter_, "");

  // set the initial inner cycle
  S_->set_cycle(tag_old_, 0);

  // loop through the inner steps til done
  for (int k = 0; k != subcycled_count_; ++k) {
    double t_inner_start = t_start + dt_inner * k;
    double t_inner_end = t_start + dt_inner * (k + 1);

    S_->set_time(tag_old_, t_inner_start);
    S_->set_time(tag_new_, t_inner_end);
    S_->set_cycle(tag_old_, k);
    S_->set_cycle(tag_new_, k + 1);

    bool fail = PK_MixinExplicit<Base_t>::AdvanceStep(tag_old_, tag_new_);
    AMANZI_ASSERT(!fail);
    this->CommitStep(tag_old_, tag_new_);
  }

  // copy from inner result to outer result
  this->StateToState(tag_new_, tag_new);
  return false;
};

} // namespace Amanzi

#endif
