/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*!

Evaluator for thermal conductivity allows us to select a proper model.
The variety of available models allows to run the energy PK by itself or in
coupling with flow PK.
The single-phase model accepts the following parameters:

.. admonition:: tcm_one_phase-spec

  * `"TCM0`" ``[list]`` defines a reginal model and its parameters.

  * `"thermal conductivity type`" ``[string]`` is the name of conductivity model for
    a mixture of phases. Available two-phase models are `"two-phase Peters-Lidard`",
    and `"two-phase wet/dry`". Available one-phase model is `"one-phase polynomial`".

  * `"regions`" ``[Array(string)]`` specifies regions where model is applied.

  * `"solid phase`" [list] define a thermal conductivity model for solid phase.

    * `"eos type`" ``[stirng]`` defines EOS model. Available options are `"constant`"
      and `"salt`".

    * `"reference conductivity`" ``[double]`` defines constant or reference conductivity.

    * `"reference temperature`" ``[double]`` defines temperature at which reference
      conductivity of liquid is calculated. Default value is 298.15 [K].

    * `"Sutherland constant`" ``[double]`` defines parameter in temperature dependence.

  * `"gas phase`" [list] define a thermal conductivity model for gas phase.

    * `"eos type`" ``[stirng]`` defines EOS model. Available options are `"constant`"
      and `"ideal gas`".

    * `"reference conductivity`" ``[double]`` defines constant or reference conductivity.

    * `"Sutherland constant`" ``[double]`` defines parameter in temperature dependence

  * `"liquid phase`" [list] define a thermal conductivity model for liquid phase.

    * `"eos type`" ``[stirng]`` defines EOS model. Available options are `"constant`"
      and `"liquid water`".

    * `"reference conductivity`" ``[double]`` defines constant or reference conductivity.
      Default value is 0.6065 [W/m/K].

    * `"reference temperature`" ``[double]`` defines temperature at which reference
      conductivity of liquid is calculated. Default value is 298.15 [K].

    * `"polynomial expansion`" ``[Array(double)]`` collect coefficients in the quadratic
      representation of the thermal conductivity of liquid with respect to the dimensionless
      parameter T/Tref.

.. code-block:: xml

  <ParameterList name="_ENERGY">  <!-- parent list -->
  <ParameterList name="thermal conductivity evaluator">
    <ParameterList name="TCM0">
      <Parameter name="regions" type="Array(string)" value="{EntireDomain}"/>
      <Parameter name="thermal conductivity type" type="string" value="one-phase polynomial"/>
      <ParameterList name="solid phase">
        <Parameter name="eos type" type="string" value="constant"/>
        <Parameter name="reference conductivity" type="double" value="0.2"/>
        <Parameter name="reference temperature" type="double" value="298.15"/>
      </ParameterList>
      <ParameterList name="liquid phase">
        <Parameter name="eos type" type="string" value="polynomial"/>
        <Parameter name="reference conductivity" type="double" value="0.6"/>
        <Parameter name="polynomial expansion" type="Array(double)" value="{-1.48445, 4.12292, -1.63866}"/>
      </ParameterList>
      <ParameterList name="gas phase">
      </ParameterList>
    </ParameterList>
  </ParameterList>
  </ParameterList>

*/

#ifndef AMANZI_ENERGY_TCM_EVALUATOR_ONEPHASE_HH_
#define AMANZI_ENERGY_TCM_EVALUATOR_ONEPHASE_HH_

#include "MeshPartition.hh"

#include "EvaluatorSecondaryMonotype.hh"
#include "H2O_ThermalConductivity.hh"

namespace Amanzi {
namespace Energy {

// Equation of State model
class TCMEvaluator_OnePhase
  : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
 public:
  TCMEvaluator_OnePhase(Teuchos::RCP<const AmanziMesh::Mesh> mesh, Teuchos::ParameterList& plist);
  TCMEvaluator_OnePhase(const TCMEvaluator_OnePhase& other);

  // required inteface functions
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override;

  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override;

  virtual void EnsureCompatibility_Units_(State& S) override;

 protected:
  std::vector<Teuchos::RCP<AmanziEOS::EOS_ThermalConductivity>> tc_liq_, tc_solid_;
  Teuchos::RCP<Functions::MeshPartition> partition_;

  // Keys for fields dependencies
  Key pressure_key_, temperature_key_, porosity_key_;
};

} // namespace Energy
} // namespace Amanzi

#endif
