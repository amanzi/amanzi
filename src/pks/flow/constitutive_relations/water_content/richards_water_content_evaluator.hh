/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
//! RichardsWaterContentEvaluator: water content without vapor

/*
  ATS is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*!
Evaluator type: `"richards water content`"

Evaluates water content in cell E.

.. math::
  \Theta = \phi n_{{liq}} s_{{liq}} |E|

* `"my key`" ``[string]`` **DOMAIN_water_content** Set by code, not user. [mol]
* `"porosity key`" ``[string]`` **DOMAIN_porosity** Names the porosity variable. [-]
* `"saturation liquid key`" ``[string]`` **DOMAIN_saturation_liquid** Names the saturation variable. [-]
* `"molar density liquid key`" ``[string]`` **DOMAIN_molar_density_liquid** Names the density variable. [mol m^-3]
* `"cell volume key`" ``[string]`` **DOMAIN_cell_volume** Names the cell volume variable. [m^3]

Note that in the defaults, DOMAIN is determined from the name of the evaluated data, which is set by the name of the list.

Example:

.. code-block:: xml

  <ParameterList name="water_content">
    <Parameter name="evaluator type" type="string" value="richards water content"/>
  </ParameterList>

*/

#ifndef AMANZI_FLOW_RICHARDS_WATER_CONTENT_EVALUATOR_HH_
#define AMANZI_FLOW_RICHARDS_WATER_CONTENT_EVALUATOR_HH_

#include "factory.hh"
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
