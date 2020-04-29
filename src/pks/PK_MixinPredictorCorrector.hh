/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon (coonet@ornl.gov)
*/

//! A mixin class with default implementations of methods for a
//! predictor-corrector scheme.

/*!

``PK_MixinPredictorCorrector`` is a mixin class providing functionality for
predictor-corrector integrated PKs.  Manages the creation of intermediate data
and AdvanceStep().

* `"initial time step`" ``[double]``

  The initial time step size, in seconds.

* `"time integrator`" ``[time-integrator-spec]``  See
``TimeIntegratorExplicit_``.
* `"predictor`" ``[time-integrator-spec]``  See ``TimeIntegratorExplicit_``.

The predictor scheme must be an explicit time integration spec, and the
corrector an implciit time integration spec.

*/

#ifndef AMANZI_PK_MIXIN_PREDICTOR_CORRECTOR_HH_
#define AMANZI_PK_MIXIN_PREDICTOR_CORRECTOR_HH_

#include "Teuchos_ParameterList.hpp"

#include "Key.hh"
#include "PK.hh"

namespace Amanzi {

template <class Base_t>
class PK_MixinPredictorCorrector
  : public PK_MixinImplicit<PK_MixinExplicit<Base_t>> {
 public:
  using PK_MixinImplicit<PK_MixinExplicit<Base_t>>::PK_MixinImplicit;

  bool ModifyPredictor(double h, Teuchos::RCP<const TreeVector> u0,
                       Teuchos::RCP<TreeVector> u);

  bool is_implicit() const { return true; }
  bool is_explicit() const { return true; }

 protected:
  using Base_t::plist_;
};

template <class Base_t>
bool
PK_MixinPredictorCorrector<Base_t>::ModifyPredictor(
  double h, Teuchos::RCP<const TreeVector> u0, Teuchos::RCP<TreeVector> u)
{
  Teuchos::RCP<Teuchos::ParameterList> ti_plist;
  if (!PK_MixinExplicit<Base_t>::time_stepper_.get()) {
    ti_plist = Teuchos::sublist(plist_, "time integrator");
    auto pre_plist = Teuchos::sublist(plist_, "predictor");
    plist_->set("time integrator", plist_->sublist("predictor"));
  }

  bool fail = PK_MixinExplicit<Base_t>::AdvanceStep(
    PK_MixinExplicit<Base_t>::tag_old_, PK_MixinExplicit<Base_t>::tag_new_);

  if (ti_plist.get()) {
    plist_->set("predictor", plist_->sublist("time integrator"));
    plist_->set("time integrator", *ti_plist);
  }
  return true;
}

} // namespace Amanzi

#endif
