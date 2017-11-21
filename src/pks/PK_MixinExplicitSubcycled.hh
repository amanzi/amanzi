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

#ifndef AMANZI_PK_MIXIN_EXPLICIT_SUBCYCLED_HH_
#define AMANZI_PK_MIXIN_EXPLICIT_SUBCYCLED_HH_

#include "Teuchos_ParameterList.hpp"

#include "Key.hh"
#include "PK_MixinExplicit.hh"
#include "PK.hh"

namespace Amanzi {

template<class Base_t>
class PK_MixinExplicitSubcycled
    : public PK_MixinExplicit<Base_t> {
 public:
  PK_MixinExplicitSubcycled(const Teuchos::RCP<Teuchos::ParameterList>& pk_tree,
                            const Teuchos::RCP<Teuchos::ParameterList>& global_plist,
                            const Teuchos::RCP<State>& S,
                            const Teuchos::RCP<TreeVector>& solution);

  void Setup();
  bool AdvanceStep(const Key& tag_old, const Key& tag_new);

  double get_dt() { return 1.e10; }

 protected:
  using PK_MixinExplicit<Base_t>::plist_;
  using PK_MixinExplicit<Base_t>::S_;
  int subcycled_count_;
}


template<class Base_t>
PK_MixinExplicitSubcycled(const Teuchos::RCP<Teuchos::ParameterList>& pk_tree,
                          const Teuchos::RCP<Teuchos::ParameterList>& global_plist,
                          const Teuchos::RCP<State>& S,
                          const Teuchos::RCP<TreeVector>& solution) :
    PK_MixinExplicit<Base_t>(pk_tree, global_plist, S, solution)
{
  subcycled_count_ = plist_->template get<double>("subcycling substeps count");
}


template<class Base_t>
void 
PK_MixinExplicit<Base_t>::Setup()
{
  S_->template Require<double>("time", this->name()+"_explicit");
}

template<class Base_t>
bool 
PK_MixinExplicitSubcycled<Base_t>::AdvanceStep(const Key& tag_old, const Key& tag_new)
{
  // take a bdf timestep
  double dt_solver;
  bool fail = time_stepper_->TimeStep(dt, dt_solver, solution_);

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
      time_stepper_->CommitSolution(dt_, solution_, valid);
      dt_ = 0.5*dt_;
    }
  } else {
    // take the decreased timestep size
    dt_ = dt_solver;
  }

  return fail;
};

} // namespace Amanzi

#endif
