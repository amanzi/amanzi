/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*!

Evaluator for thermal conductivity allows us to select a proper model.
The variety of available models allows to run the energy PK by itself or in
coupling with flow PK.
The structure of the thermal conductivity list resembles that of a field
evaluator list in state.
The two-phase model accepts the following parameters.

* `"thermal conductivity parameters`" [list] defines a model and its parameters.

* `"thermal conductivity type`" [string] is the name of a conductivity model in the
  list of registered models. Available two-phase models are `"two-phase Peters-Lidard`",
  and `"two-phase wet/dry`". Available one-phase model is `"one-phase polynomial`".

* `"thermal conductivity of rock`" [double] defines constant conductivity of rock.

* `"thermal conductivity of gas`" [double] defines constant conductivity of gas.

* `"thermal conductivity of liquid`" [double] defines constant conductivity of fluid.
  Default value is 0.6065 [W/m/K].

* `"unsaturated alpha`" [double] is used to define the Kersten number to interpolate
  between saturated and dry conductivities.

* `"epsilon`" [double] is needed for the case of zero saturation. Default is `"1.0e-10`".

.. code-block:: xml

  <ParameterList name="_ENERGY">  <!-- parent list -->
  <ParameterList name="thermal conductivity evaluator">
    <ParameterList name="thermal conductivity parameters">
      <Parameter name="thermal conductivity type" type="string" value="two-phase Peters-Lidard"/>
      <Parameter name="thermal conductivity of rock" type="double" value="0.2"/>
      <Parameter name="thermal conductivity of gas" type="double" value="0.02"/>
      <Parameter name="thermal conductivity of liquid" type="double" value="0.6065"/>

      <Parameter name="unsaturated alpha" type="double" value="1.0"/>
      <Parameter name="epsilon" type="double" value="1.e-10"/>
    </ParameterList>
  </ParameterList>
  </ParameterList>

The single-phase model accepts some of the parameters defined above (see the example)
and a few additional parameters.

* `"reference temperature`" [double] defines temperature at which reference conductivity
  of liquid is calculated. Default value is 298.15 [K].

* `"polynomial expansion`" [Array(double)] collect coefficients in the quadratic representation of the
  thermal conductivity of liquid with respect to the dimensionless parameter T/Tref.

.. code-block:: xml

  <ParameterList name="_ENERGY">  <!-- parent list -->
  <ParameterList name="thermal conductivity evaluator">
    <ParameterList name="thermal conductivity parameters">
      <Parameter name="thermal conductivity type" type="string" value="one-phase polynomial"/>
      <Parameter name="thermal conductivity of rock" type="double" value="0.2"/>
      <Parameter name="reference temperature" type="double" value="298.15"/>
      <Parameter name="polynomial expansion" type="Array(double)" value="{-1.48445, 4.12292, -1.63866}"/>
    </ParameterList>
  </ParameterList>
  </ParameterList>

*/

#ifndef PK_ENERGY_TCM_FACTORY_TWOPHASE_HH_
#define PK_ENERGY_TCM_FACTORY_TWOPHASE_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "TCM_TwoPhase.hh"
#include "Factory.hh"

namespace Amanzi {
namespace Energy {

class TCMFactory_TwoPhase : public Utils::Factory<TCM_TwoPhase> {
 public:
  Teuchos::RCP<TCM_TwoPhase> CreateTCM(Teuchos::ParameterList& plist);
};

} // namespace Energy
} // namespace Amanzi

#endif
