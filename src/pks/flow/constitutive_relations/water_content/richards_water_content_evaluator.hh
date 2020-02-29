/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/*
  ATS is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)
*/
//! The Richards water content evaluator is an algebraic evaluator for liquid only water content

/*!
  Generated via evaluator_generator with:
Richards water content evaluator: the standard form as a function of liquid saturation.

.. math::
  Theta = n * s * phi * cell volume

``[field-evaluator-type-richards-water-content-spec]``

* `"porosity key`" ``[string]`` **DOMAIN-porosity** 
* `"molar density liquid key`" ``[string]`` **DOMAIN-molar_density_liquid** 
* `"saturation liquid key`" ``[string]`` **DOMAIN-saturation_liquid** 
* `"cell volume key`" ``[string]`` **DOMAIN-cell_volume**

EVALUATORS:
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
