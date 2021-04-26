/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)
*/
//! Richards water content evaluator: the standard form as a function of liquid saturation.
/*!

.. math::
  \Theta = n s \phi V

Specified with evaluator type: `"richards water content`"

.. _field-evaluator-type-richards-water-content-spec:
.. admonition:: field-evaluator-type-richards-water-content-spec

   DEPENDENCIES:

   - `"porosity`"
   - `"molar density liquid`"
   - `"saturation liquid`"
   - `"cell volume`"

*/

#ifndef AMANZI_FLOW_RICHARDS_WATER_CONTENT_EVALUATOR_HH_
#define AMANZI_FLOW_RICHARDS_WATER_CONTENT_EVALUATOR_HH_

#include "Factory.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {

class RichardsWaterContentModel;

class RichardsWaterContentEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  explicit
  RichardsWaterContentEvaluator(Teuchos::ParameterList& plist);
  RichardsWaterContentEvaluator(const RichardsWaterContentEvaluator& other);

  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

  Teuchos::RCP<RichardsWaterContentModel> get_model() { return model_; }

 protected:
  void InitializeFromPlist_();

  Key phi_key_;
  Key sl_key_;
  Key nl_key_;
  Key cv_key_;

  Teuchos::RCP<RichardsWaterContentModel> model_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,RichardsWaterContentEvaluator> reg_;

};

} //namespace
} //namespace
} //namespace

#endif
