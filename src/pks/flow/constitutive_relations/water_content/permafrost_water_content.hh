/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
//! RichardsWaterContentEvaluator: water content without vapor

/*
  ATS is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*!
Evaluator type: `"permafrost water content`"

Evaluates water content in cell E.

.. math::
  \Theta = \phi (n_{{ice}} s_{{ice}} + n_{{liq}} s_{{liq}} + n_{{gas}} s_{{gas}} \omega) |E|

* `"my key`" ``[string]`` **DOMAIN_water_content** Set by code, not user. [mol]
* `"porosity key`" ``[string]`` **DOMAIN_porosity** Names the porosity variable. [-]
* `"saturation ice key`" ``[string]`` **DOMAIN_saturation_ice** Names the ice saturation variable. [-]
* `"saturation liquid key`" ``[string]`` **DOMAIN_saturation_liquid** Names the liquid saturation variable. [-]
* `"saturation gas key`" ``[string]`` **DOMAIN_saturation_gas** Names the gas saturation variable. [-]
* `"molar density ice key`" ``[string]`` **DOMAIN_molar_density_ice** Names the ice density variable. [mol m^-3]
* `"molar density liquid key`" ``[string]`` **DOMAIN_molar_density_liquid** Names the liquid density variable. [mol m^-3]
* `"molar density gas key`" ``[string]`` **DOMAIN_molar_density_gas** Names the gas density variable. [mol m^-3]
* `"mol fraction vapor in gas key`" ``[string]`` **DOMAIN_mol_frac_gas** Names the molar fraction of water vapor in the gas phase variable. [-]
* `"cell volume key`" ``[string]`` **DOMAIN_cell_volume** Names the cell volume variable. [m^3]

Note that in the defaults, DOMAIN is determined from the name of the evaluated data, which is set by the name of the list.

Example:

.. code-block:: xml

  <ParameterList name="water_content">
    <Parameter name="evaluator type" type="string" value="permafrost water content"/>
  </ParameterList>

*/


#ifndef AMANZI_PERMAFROST_WATER_CONTENT_HH_
#define AMANZI_PERMAFROST_WATER_CONTENT_HH_

#include "Teuchos_ParameterList.hpp"

#include "factory.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {

class PermafrostWaterContent : public SecondaryVariableFieldEvaluator {

 public:
  explicit
  PermafrostWaterContent(Teuchos::ParameterList& wc_plist);
  PermafrostWaterContent(const PermafrostWaterContent& other);

  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

 protected:
  Key phi_key_;
  Key sl_key_;
  Key mdl_key_;
  Key si_key_;
  Key mdi_key_;
  Key sg_key_;
  Key mdg_key_;
  Key mfg_key_;
  Key cv_key_;

private:
  static Utils::RegisteredFactory<FieldEvaluator,PermafrostWaterContent> reg_;

};

} // namespace
} // namespace
} // namespace

#endif
