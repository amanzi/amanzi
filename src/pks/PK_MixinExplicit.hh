/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/*
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! A mixin class with default implementations of methods for an explicotly integrated PK.

/*!

``PK_MixinExplicit`` is a mixin class providing functionality for ixplicitly
integrated PKs.  Manages the creation of intermediate data and AdvanceStep().

* `"initial time step`" ``[double]``

  The initial time step size, in seconds.

* `"time integrator`" ``[time-integrator-spec]``  See ``TimeIntegratorExplicit_``.

*/

#ifndef AMANZI_PK_MIXIN_EXPLICIT_HH_
#define AMANZI_PK_MIXIN_EXPLICIT_HH_

#include "Teuchos_ParameterList.hpp"

#include "Key.hh"
#include "Explicit_TI_RK.hh"
#include "PK.hh"

namespace Amanzi {

template<class Base_t>
class PK_MixinExplicit : public Base_t, public Explicit_TI::fnBase<TreeVector> {
 public:
  PK_MixinExplicit(const Teuchos::RCP<Teuchos::ParameterList>& pk_tree,
                   const Teuchos::RCP<Teuchos::ParameterList>& global_plist,
                   const Teuchos::RCP<State>& S,
                   const Teuchos::RCP<TreeVector>& solution);

  void Setup(const TreeVector& soln);
  bool AdvanceStep(const Key& tag_old, const Teuchos::RCP<TreeVector>& soln_old,
                   const Key& tag_new, const Teuchos::RCP<TreeVector>& soln_new);

  double get_dt() { return dt_; }


 protected:
  // timestep size
  double dt_;

  Key tag_old_, tag_new_;
  
  // timestep algorithm
  Teuchos::RCP<Explicit_TI::RK<TreeVector> > time_stepper_;

  // previous solution
  Teuchos::RCP<TreeVector> solution_intermediate_;

  // tag at which evaluators are needed
  Key dudt_tag_;
  bool method_requires_intermediate_;
  
  using Base_t::plist_;
  using Base_t::S_;
  using Base_t::vo_;
};

template<class Base_t>
PK_MixinExplicit<Base_t>::PK_MixinExplicit(const Teuchos::RCP<Teuchos::ParameterList>& pk_tree,
                         const Teuchos::RCP<Teuchos::ParameterList>& global_plist,
                         const Teuchos::RCP<State>& S,
                         const Teuchos::RCP<TreeVector>& solution)
    : Base_t(pk_tree, global_plist, S, solution),
      method_requires_intermediate_(false)
{
  // initial timestep
  dt_ = plist_->template get<double>("initial time step", 1.);
};


template<class Base_t>
void
PK_MixinExplicit<Base_t>::Setup(const TreeVector& soln)
{
  Base_t::Setup(soln);

  std::string methodname = plist_->sublist("time integrator")
                           .template get<std::string>("RK method");
  
  // potentially create an intermediate tag for multistage algorithms
  method_requires_intermediate_ = methodname != "forward euler";
  if (method_requires_intermediate_) {
    dudt_tag_ = this->name()+" explicit ti intermediate";
    solution_intermediate_ = Teuchos::rcp(new TreeVector(soln, INIT_MODE_NOALLOC));
    this->SolutionToState(*solution_intermediate_, dudt_tag_, "");
    S_->template Require<double>("time", dudt_tag_, "time");
  }

  tag_old_ = "";
  tag_new_ = "next";
}


template<class Base_t>
bool 
PK_MixinExplicit<Base_t>::AdvanceStep(const Key& tag_old, const Teuchos::RCP<TreeVector>& soln_old,
        const Key& tag_new, const Teuchos::RCP<TreeVector>& soln_new)
{
  tag_old_ = tag_old;
  tag_new_ = tag_new;

  if (solution_intermediate_ == Teuchos::null) {
    solution_intermediate_ = soln_old;
    dudt_tag_ = tag_old;
  }

  // -- ensure state vectors are pushed into solution vectors
  this->StateToSolution(*soln_old, tag_old, "");
  this->StateToSolution(*solution_intermediate_, dudt_tag_, "");
  this->StateToSolution(*soln_new, tag_new, "");
  
  
  // -- instantiate time stepper if needed
  if (!time_stepper_.get()) {
    Teuchos::ParameterList& ti_plist = plist_->sublist("time integrator");
    ti_plist.set("initial time", S_->time());
    time_stepper_ = Teuchos::rcp(new Explicit_TI::RK<TreeVector>(*this, ti_plist, *soln_new));
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


  // take a timestep
  time_stepper_->TimeStep(t_old, dt, *soln_old, *solution_intermediate_);
  *soln_new = *solution_intermediate_;
  return false;
};
  


} // namespace Amanzi

#endif
