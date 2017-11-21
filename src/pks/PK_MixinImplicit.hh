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

#ifndef AMANZI_PK_MIXIN_IMPLICIT_HH_
#define AMANZI_PK_MIXIN_IMPLICIT_HH_

#include "Teuchos_ParameterList.hpp"

#include "Key.hh"
#include "BDF1_TI.hh"
#include "PK.hh"

namespace Amanzi {

template<class Base_t>
class PK_MixinImplicit
    : public Base_t,
      public BDFFnBase<TreeVector> {
 public:
  PK_MixinImplicit(const Teuchos::RCP<Teuchos::ParameterList>& pk_tree,
                   const Teuchos::RCP<Teuchos::ParameterList>& global_plist,
                   const Teuchos::RCP<State>& S,
                   const Teuchos::RCP<TreeVector>& solution);

  void Setup(const TreeVector& soln);
  bool AdvanceStep(const Key& tag_old, const Teuchos::RCP<TreeVector>& soln_old,
                   const Key& tag_new, const Teuchos::RCP<TreeVector>& soln_new);
  void CommitStep(const Key& tag_old, const Teuchos::RCP<TreeVector>& soln_old,
                  const Key& tag_new, const Teuchos::RCP<TreeVector>& soln_new);

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
  Teuchos::RCP<BDF1_TI<TreeVector,TreeVectorSpace> > time_stepper_;

  using Base_t::plist_;
  using Base_t::S_;
  using Base_t::vo_;
};

template<class Base_t>
PK_MixinImplicit<Base_t>::PK_MixinImplicit(const Teuchos::RCP<Teuchos::ParameterList>& pk_tree,
                         const Teuchos::RCP<Teuchos::ParameterList>& global_plist,
                         const Teuchos::RCP<State>& S,
                         const Teuchos::RCP<TreeVector>& solution)
    : Base_t(pk_tree, global_plist, S, solution)
{
  // initial timestep
  dt_ = plist_->template get<double>("initial time step", 1.);
};


template<class Base_t>
void
PK_MixinImplicit<Base_t>::Setup(const TreeVector& soln)
{
  Base_t::Setup(soln);

  // preconditioner assembly
  assemble_preconditioner_ = plist_->template get<bool>("assemble preconditioner", true);

  tag_old_ = "";
  tag_new_ = "next";
}


template<class Base_t>
bool 
PK_MixinImplicit<Base_t>::AdvanceStep(const Key& tag_old, const Teuchos::RCP<TreeVector>& soln_old,
        const Key& tag_new, const Teuchos::RCP<TreeVector>& soln_new)
{
  tag_old_ = tag_old;
  tag_new_ = tag_new;

  // -- ensure state vectors are pushed into solution vectors
  this->StateToSolution(*soln_old, tag_old, "");
  this->StateToSolution(*soln_new, tag_new, "");

  if (!time_stepper_.get()) {
    // -- instantiate time stepper
    Teuchos::ParameterList& bdf_plist = plist_->sublist("time integrator");
    bdf_plist.set("initial time", S_->time(tag_old));
    time_stepper_ = Teuchos::rcp(new BDF1_TI<TreeVector,TreeVectorSpace>(*this, bdf_plist, soln_new));

    // -- initialize time derivative
    auto solution_dot = Teuchos::rcp(new TreeVector(*soln_old, INIT_MODE_ZERO));

    // -- set initial state
    *soln_old = *soln_new;
    time_stepper_->SetInitialState(S_->time(tag_old), soln_old, solution_dot);
  }

  double t_new = S_->time(tag_new);
  double t_old = S_->time(tag_old);
  double dt = t_new - t_old;  
  
  Teuchos::OSTab out = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_HIGH))
    *vo_->os() << "----------------------------------------------------------------" << std::endl
               << "Advancing: t0 = " << t_old
               << " t1 = " << t_new << " h = " << dt << std::endl
               << "----------------------------------------------------------------" << std::endl;


  // take a bdf timestep
  double dt_solver;
  bool fail = time_stepper_->TimeStep(dt, soln_old, soln_new, dt_solver);

  if (!fail) {
    // check step validity
    bool valid = this->ValidStep(tag_old, soln_old, tag_new, soln_new);
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
      dt_ = 0.5*dt_;
    }
  } else {
    // take the decreased timestep size
    dt_ = dt_solver;
  }

  return fail;
};

  
// -- Commit any secondary (dependent) variables.
template<class Base_t>
void
PK_MixinImplicit<Base_t>::CommitStep(const Key& tag_old, const Teuchos::RCP<TreeVector>& soln_old,
        const Key& tag_new, const Teuchos::RCP<TreeVector>& soln_new)
{
  Base_t::CommitStep(tag_old, soln_old, tag_new, soln_new);  

  double dt = S_->time(tag_new) - S_->time(tag_old);
  if (dt > 0. && time_stepper_ != Teuchos::null)
    time_stepper_->CommitSolution(dt, soln_new, true);
}



} // namespace Amanzi

#endif
