/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)
*/
//! Water content for liquid + water vapor.
/*!

.. math::
  \Theta = (n_l s_l + n_i s_i) \phi V


Specified with evaluator type: `"liquid+ice water content`"

.. _field-evaluator-type-liquid-ice-water-content-spec:
.. admonition:: field-evaluator-type-liquid-ice-water-content-spec

   DEPENDENCIES:

   - `"porosity`"
   - `"molar density liquid`"
   - `"molar density ice`"
   - `"saturation liquid`"
   - `"saturation ice`"
   - `"cell volume`"

*/

#ifndef AMANZI_FLOW_LIQUID_ICE_WATER_CONTENT_EVALUATOR_HH_
#define AMANZI_FLOW_LIQUID_ICE_WATER_CONTENT_EVALUATOR_HH_

#include "Factory.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {

class LiquidIceWaterContentModel;

class LiquidIceWaterContentEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  explicit
  LiquidIceWaterContentEvaluator(Teuchos::ParameterList& plist);
  LiquidIceWaterContentEvaluator(const LiquidIceWaterContentEvaluator& other);

  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

  Teuchos::RCP<LiquidIceWaterContentModel> get_model() { return model_; }

 protected:
  void InitializeFromPlist_();

  Key phi_key_;
  Key sl_key_;
  Key nl_key_;
  Key si_key_;
  Key ni_key_;
  Key cv_key_;

  Teuchos::RCP<LiquidIceWaterContentModel> model_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,LiquidIceWaterContentEvaluator> reg_;

};

} //namespace
} //namespace
} //namespace

#endif
