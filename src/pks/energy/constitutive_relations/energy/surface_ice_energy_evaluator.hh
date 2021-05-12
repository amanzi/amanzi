/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)
*/
//! Energy content for a surface water, partially frozen system.
/*!

The energy associated with ponded water, in [KJ], given by:

.. math::
  E = V * ( \eta h u_l n_l + (1 - \eta) h u_i n_i )

Specified with evaluator type: `"surface ice energy`"

.. _field-evaluator-type-surface-ice-energy-spec:
.. admonition:: field-evaluator-type-surface-ice-energy-spec

   DEPENDENCIES:

   - `"ponded depth`"  Height of water above the land surface [m]
   - `"unfrozen fraction`"  The fraction of unfrozen water ranges from 0 to 1. [-]
   - `"molar density liquid`" [mol m^-3]
   - `"internal energy liquid`" [KJ mol^-1]
   - `"molar density ice`" [mol m^-3]
   - `"internal energy ice`" [KJ mol^-1]
   - `"cell volume`" [m^2]

*/

#pragma once

#include "Factory.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace Energy {
namespace Relations {

class SurfaceIceEnergyModel;

class SurfaceIceEnergyEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  explicit
  SurfaceIceEnergyEvaluator(Teuchos::ParameterList& plist);
  SurfaceIceEnergyEvaluator(const SurfaceIceEnergyEvaluator& other);

  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

  Teuchos::RCP<SurfaceIceEnergyModel> get_model() { return model_; }

 protected:
  void InitializeFromPlist_();

  Key h_key_;
  Key eta_key_;
  Key nl_key_;
  Key ul_key_;
  Key ni_key_;
  Key ui_key_;
  Key cv_key_;

  Teuchos::RCP<SurfaceIceEnergyModel> model_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,SurfaceIceEnergyEvaluator> reg_;

};

} //namespace
} //namespace
} //namespace

