/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)
*/
//! Water content for liquid + water vapor.
/*!

.. math::
  \Theta = (n_l s_l + n_g s_g \omega) \phi V


Specified with evaluator type: `"liquid+gas water content`"

.. _field-evaluator-type-liquid-gas-water-content-spec:
.. admonition:: field-evaluator-type-liquid-gas-water-content-spec

   DEPENDENCIES:

   - `"porosity`"
   - `"molar density liquid`"
   - `"molar density gas`"
   - `"saturation liquid`"
   - `"saturation gas`"
   - `"mol frac gas`"
   - `"cell volume`"

*/

#ifndef AMANZI_FLOW_LIQUID_GAS_WATER_CONTENT_EVALUATOR_HH_
#define AMANZI_FLOW_LIQUID_GAS_WATER_CONTENT_EVALUATOR_HH_

#include "Factory.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {

class LiquidGasWaterContentModel;

class LiquidGasWaterContentEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  explicit
  LiquidGasWaterContentEvaluator(Teuchos::ParameterList& plist);
  LiquidGasWaterContentEvaluator(const LiquidGasWaterContentEvaluator& other);

  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

  Teuchos::RCP<LiquidGasWaterContentModel> get_model() { return model_; }

 protected:
  void InitializeFromPlist_();

  Key phi_key_;
  Key sl_key_;
  Key nl_key_;
  Key sg_key_;
  Key ng_key_;
  Key omega_key_;
  Key cv_key_;

  Teuchos::RCP<LiquidGasWaterContentModel> model_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,LiquidGasWaterContentEvaluator> reg_;

};

} //namespace
} //namespace
} //namespace

#endif
