/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)
*/
//! Energy content evaluator for a standard Richards equation, including energy in the gas phase.
/*!

Calculates energy, in [KJ], via the equation:

.. math::
  E = V * ( \phi (u_l s_l n_l + u_g s_g n_g)  + (1-\phi_0) u_r \rho_r )

Specified with evaluator type: `"liquid+gas energy`"

Note this equation assumes that porosity is compressible, but is based on the
uncompressed rock grain density (not bulk density).  This means that porosity
is the base, uncompressible value when used with the energy in the grain, but
the larger, compressible value when used with the energy in the water.

.. _field-evaluator-type-liquid-gas-energy-spec:
.. admonition:: field-evaluator-type-liquid-gas-energy-spec

   DEPENDENCIES:

   - `"porosity`" The porosity, including any compressibility. [-]
   - `"base porosity`" The uncompressed porosity (note this may be the same as
     porosity for incompressible cases) [-]
   - `"molar density liquid`" [mol m^-3]
   - `"saturation liquid`" [-]
   - `"internal energy liquid`" [KJ mol^-1]
   - `"molar density gas`" [mol m^-3]
   - `"saturation gas`" [-]
   - `"internal energy gas`" [KJ mol^-1]
   - `"density rock`" Units may be either [kg m^-3] or [mol m^-3]
   - `"internal energy rock`" Units may be either [KJ kg^-1] or [KJ mol^-1],
     but must be consistent with the above density.
   - `"cell volume`" [m^3]

*/

#pragma once

#include "Factory.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace Energy {
namespace Relations {

class LiquidGasEnergyModel;

class LiquidGasEnergyEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  explicit
  LiquidGasEnergyEvaluator(Teuchos::ParameterList& plist);
  LiquidGasEnergyEvaluator(const LiquidGasEnergyEvaluator& other);

  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

  Teuchos::RCP<LiquidGasEnergyModel> get_model() { return model_; }

 protected:
  void InitializeFromPlist_();

  Key phi_key_;
  Key phi0_key_;
  Key sl_key_;
  Key nl_key_;
  Key ul_key_;
  Key sg_key_;
  Key ng_key_;
  Key ug_key_;
  Key rho_r_key_;
  Key ur_key_;
  Key cv_key_;

  Teuchos::RCP<LiquidGasEnergyModel> model_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,LiquidGasEnergyEvaluator> reg_;

};

} //namespace
} //namespace
} //namespace

