/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon (coonet@ornl.gov)
*/

//! A mixin class with default implementations of methods for an explicotly
//! integrated PK.

/*!

``PK_MixinExplicit`` is a mixin class providing functionality for ixplicitly
integrated PKs.  Manages the creation of intermediate data and AdvanceStep().

* `"initial time step`" ``[double]``

  The initial time step size, in seconds.

* `"time integrator`" ``[time-integrator-spec]``  See
``TimeIntegratorExplicit_``.

*/

#ifndef AMANZI_PK_MIXIN_EXPLICIT_HH_
#define AMANZI_PK_MIXIN_EXPLICIT_HH_

#include "Teuchos_ParameterList.hpp"

#include "Explicit_TI_RK.hh"
#include "Key.hh"
#include "TreeVector.hh"
#include "PK.hh"

namespace Amanzi {

template <class Base_t>
class PK_MixinExplicit : public Base_t {
 public:
  PK_MixinExplicit(const Teuchos::RCP<Teuchos::ParameterList>& pk_tree,
                   const Teuchos::RCP<Teuchos::ParameterList>& global_plist,
                   const Teuchos::RCP<State>& S);

  virtual ~PK_MixinExplicit() =
    default; // here to make this polymorphic and therefore castable

  void Setup();
  bool AdvanceStep(const Key& tag_old, const Key& tag_new);

  double get_dt() { return dt_; }

 protected:
  // timestep size
  double dt_;

  Key tag_old_, tag_new_;

  // timestep algorithm
  Teuchos::RCP<Explicit_TI::RK<TreeVector>> time_stepper_;

  // tag at which evaluators are needed
  Key tag_inter_;

  using Base_t::S_;
  using Base_t::plist_;
  using Base_t::vo_;
};

template <class Base_t>
PK_MixinExplicit<Base_t>::PK_MixinExplicit(
  const Teuchos::RCP<Teuchos::ParameterList>& pk_tree,
  const Teuchos::RCP<Teuchos::ParameterList>& global_plist,
  const Teuchos::RCP<State>& S)
  : Base_t(pk_tree, global_plist, S)
{
  // initial timestep
  dt_ = plist_->template get<double>("initial time step", 1.);
};

template <class Base_t>
void
PK_MixinExplicit<Base_t>::Setup()
{
  Base_t::Setup();

  // create an intermediate tag for derivative evaluation, potentially at
  // multiple stages
  tag_inter_ = this->name() + " explicit ti intermediate";
  this->SolutionToState(tag_inter_, "");
  S_->template Require<double>("time", tag_inter_, "time");
  S_->template Require<int>("cycle", tag_inter_, "cycle");

  // set up tags
  tag_old_ = "";
  tag_new_ = "next";

  // require data at the new and old times
  this->SolutionToState(tag_new_, "");
  this->SolutionToState(tag_old_, "");
}

template <class Base_t>
bool
PK_MixinExplicit<Base_t>::AdvanceStep(const Key& tag_old, const Key& tag_new)
{
  // my local tags, used in physics PKs?  Can we get rid of these? --etc
  tag_old_ = tag_old;
  tag_new_ = tag_new;

  // times associated with those tags
  double t_new = S_->time(tag_new);
  double t_old = S_->time(tag_old);
  double dt = t_new - t_old;

  // logging
  Teuchos::OSTab out = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_HIGH))
    *vo_->os()
      << "----------------------------------------------------------------"
      << std::endl
      << "Advancing: t0 = " << t_old << " t1 = " << t_new << " h = " << dt
      << std::endl
      << "----------------------------------------------------------------"
      << std::endl;

  // create solution vectors, old and intermediate, which are pointers into
  // state data
  auto soln_space = this->SolutionSpace();

  auto soln_old = Teuchos::rcp(new TreeVector(soln_space, InitMode::NOALLOC));
  this->StateToSolution(*soln_old, tag_old, "");

  auto soln_inter = Teuchos::rcp(new TreeVector(soln_space, InitMode::NOALLOC));
  this->StateToSolution(*soln_inter, tag_inter_, "");

  // -- instantiate time stepper if needed
  if (!time_stepper_.get()) {
    Teuchos::ParameterList& ti_plist = plist_->sublist("time integrator");
    ti_plist.set("initial time", S_->time());
    auto this_as_explicit_p =
      dynamic_cast<Explicit_TI::fnBase<TreeVector>*>(this);
    AMANZI_ASSERT(this_as_explicit_p);
    time_stepper_ = Teuchos::rcp(new Explicit_TI::RK<TreeVector>(
      *this_as_explicit_p, ti_plist, *soln_inter));
  }

  // take a timestep
  time_stepper_->TimeStep(t_old, dt, *soln_old, *soln_inter);
  this->StateToState(tag_inter_, tag_new_);
  return false;
};

} // namespace Amanzi

#endif
