/* -*-  mode: c++; indent-tabs-mode: nil -*- */
//! A mixin class with default implementations of methods for an explicotly integrated PK.

/*
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/


/*!

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

  void Setup();
  bool AdvanceStep(const Key& tag_old, const Key& tag_new);

  double get_dt() { return dt_; }


 protected:
  // timestep size
  double dt_;

  // timestep algorithm
  Teuchos::RCP<Explicit_TI::RK<TreeVector> > time_stepper_;

  // previous solution
  Teuchos::RCP<TreeVector> solution_old_;
  Teuchos::RCP<TreeVector> solution_intermediate_;

  Key dudt_tag_;
  bool method_requires_intermediate_;
  
  using Base_t::plist_;
  using Base_t::S_;
  using Base_t::solution_;
  using Base_t::vo_;
};

template<class Base_t>
PK_MixinExplicit<Base_t>::PK_MixinExplicit(const Teuchos::RCP<Teuchos::ParameterList>& pk_tree,
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
PK_MixinExplicit<Base_t>::Setup()
{
  Base_t::Setup();

  std::string methodname = plist_->sublist("time integrator")
                           .template get<std::string>("RK method");
  method_requires_intermediate_ = methodname != "forward euler";
  
  // create the old solution
  solution_old_ = Teuchos::rcp(new TreeVector(*solution_, INIT_MODE_NOALLOC));
  this->SolutionToState(*solution_old_, "", "");

  // potentially create an intermediate tag for multistage algorithms
  if (method_requires_intermediate_) {
    dudt_tag_ = this->name()+" explicit ti intermediate";
    solution_intermediate_ = Teuchos::rcp(new TreeVector(*solution_, INIT_MODE_NOALLOC));
    this->SolutionToState(*solution_intermediate_, dudt_tag_, "");
    S_->template Require<double>("time", dudt_tag_, "time");
  } else {
    solution_intermediate_ = solution_old_;
  }
  
  // create the new solution
  this->SolutionToState(*solution_, "next", "");
}


template<class Base_t>
bool 
PK_MixinExplicit<Base_t>::AdvanceStep(const Key& tag_old, const Key& tag_new)
{
  if (!time_stepper_.get()) {
    // -- ensure state vectors are pushed into solution vectors
    this->StateToSolution(*solution_old_, "", "");
    this->StateToSolution(*solution_intermediate_, dudt_tag_, "");
    this->StateToSolution(*solution_, "next", "");

    // -- instantiate time stepper
    Teuchos::ParameterList& ti_plist = plist_->sublist("time integrator");
    ti_plist.set("initial time", S_->time());
    time_stepper_ = Teuchos::rcp(new Explicit_TI::RK<TreeVector>(*this, ti_plist, *solution_));
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
  time_stepper_->TimeStep(t_old, dt, *solution_old_, *solution_intermediate_);
  *solution_ = *solution_intermediate_;
  return false;
};
  


} // namespace Amanzi

#endif
