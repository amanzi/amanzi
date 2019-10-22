/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon (coonet@ornl.gov)
*/

//! A mixin class with default implementations of methods for an BDF integrated
//! PK.

/*!

``PK_MixinImplicit`` is a mixin class providing functionality for implicitly
integrated PKs.  Manages the creation of intermediate data and AdvanceStep().

* `"initial time step`" ``[double]``

  The initial time step size, in seconds.

* `"time integrator`" ``[time-integrator-spec]``  See ``TimeIntegratorBDF_``.

*/

#ifndef AMANZI_PK_MIXIN_IMPLICIT_HH_
#define AMANZI_PK_MIXIN_IMPLICIT_HH_

#include "Teuchos_ParameterList.hpp"

#include "BDF1_TI.hh"
#include "Key.hh"
#include "PK.hh"

namespace Amanzi {

template <class Base_t>
class PK_MixinImplicit : public Base_t {
 public:
  PK_MixinImplicit(const Teuchos::RCP<Teuchos::ParameterList>& pk_tree,
                   const Teuchos::RCP<Teuchos::ParameterList>& global_plist,
                   const Teuchos::RCP<State>& S);

  virtual ~PK_MixinImplicit() =
    default; // here to make this polymorphic and therefore castable

  void Setup();
  bool AdvanceStep(const Key& tag_old, const Key& tag_new);
  void CommitStep(const Key& tag_old, const Key& tag_new);

  double get_dt() { return dt_; }

 protected:
  // timestep size
  double dt_;
  bool assemble_preconditioner_;

  // tags for start and end of step Note these are here to limit how tightly
  // integrated time integration is to a PK/dag.  If time integrators worked
  // with tags we would not need these.
  Key tag_old_, tag_new_;

  // timestep algorithm
  Teuchos::RCP<BDF1_TI<TreeVector, TreeVectorSpace>> time_stepper_;

  using Base_t::S_;
  using Base_t::plist_;
  using Base_t::vo_;
};

template <class Base_t>
PK_MixinImplicit<Base_t>::PK_MixinImplicit(
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
PK_MixinImplicit<Base_t>::Setup()
{
  Base_t::Setup();

  // preconditioner assembly
  assemble_preconditioner_ =
    plist_->template get<bool>("assemble preconditioner", true);

  // set up tags
  tag_old_ = "";
  tag_new_ = "next";

  // require data at the new and old times
  this->SolutionToState(tag_new_, "");
  this->SolutionToState(tag_old_, "");
}

template <class Base_t>
bool
PK_MixinImplicit<Base_t>::AdvanceStep(const Key& tag_old, const Key& tag_new)
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

  // create solution vectors, old and new, which are pointers into state data
  auto soln_space = this->SolutionSpace();

  auto soln_old = Teuchos::rcp(new TreeVector(soln_space, InitMode::NOALLOC));
  this->StateToSolution(*soln_old, tag_old, "");

  auto soln_new = Teuchos::rcp(new TreeVector(soln_space, InitMode::NOALLOC));
  this->StateToSolution(*soln_new, tag_new, "");

  // create the time integrator if first call
  if (!time_stepper_.get()) {
    // -- instantiate time stepper
    Teuchos::ParameterList& bdf_plist = plist_->sublist("time integrator");
    bdf_plist.set("initial time", S_->time(tag_old));

    auto this_as_bdf_p = dynamic_cast<BDFFnBase<TreeVector>*>(this);
    AMANZI_ASSERT(this_as_bdf_p);
    time_stepper_ = Teuchos::rcp(new BDF1_TI<TreeVector, TreeVectorSpace>(
      *this_as_bdf_p, bdf_plist, soln_new));

    // -- initialize time derivative at 0, used in higher order BDF methods?
    // Can we remove this? --etc
    auto solution_dot = Teuchos::rcp(new TreeVector(soln_space));

    // -- set initial state
    this->StateToState(tag_old,
                       tag_new); // assures new and old time are both the IC
    time_stepper_->SetInitialState(S_->time(tag_old), soln_old, solution_dot);
  }

  // take a bdf timestep
  double dt_solver;
  bool fail = time_stepper_->TimeStep(dt, soln_old, soln_new, dt_solver);

  if (!fail) {
    // check step validity
    bool valid = this->ValidStep(tag_old, tag_new);
    if (valid) {
      // update the timestep size
      if (dt_solver < dt_ && dt_solver >= dt) {
        // We took a smaller step than we recommended, and it worked fine (not
        // suprisingly).  Likely this was due to constraints from other PKs or
        // vis.  Do not reduce our recommendation.
      } else {
        dt_ = dt_solver;
      }
    } else {
      time_stepper_->CommitSolution(dt_, soln_new, valid);
      dt_ = 0.5 * dt_;
    }
  } else {
    // take the decreased timestep size
    dt_ = dt_solver;
  }

  return fail;
};

// -- Commit any secondary (dependent) variables.
template <class Base_t>
void
PK_MixinImplicit<Base_t>::CommitStep(const Key& tag_old, const Key& tag_new)
{
  double dt = S_->time(tag_new) - S_->time(tag_old);

  // this should eventually be removed -- it is manipulating internal state of
  // the time integrator, and that state should go away --etc
  if (dt > 0. && time_stepper_ != Teuchos::null) {
    auto soln_new = Teuchos::rcp(new TreeVector(this->SolutionSpace(), InitMode::NOALLOC));
    this->StateToSolution(*soln_new, tag_new, "");
    time_stepper_->CommitSolution(dt, soln_new, true);
  }

  Base_t::CommitStep(tag_old, tag_new);
}

} // namespace Amanzi

#endif
