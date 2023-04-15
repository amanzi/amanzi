/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*!

User defines porosity models in sublist *porosity models*. 
It contains as many sublists, e.g. _SOIL1 and _SOIL2, as there are different soils. 
The porosity models are associated with non-overlapping regions. Each of the sublists (e.g. _SOIL1) 
includes a few mandatory parameters: *regions names*, *model name*, and parameters for the selected model.

* `"porosity model`" [string] specifies a model for the soil.
  The available models are `"compressible`" and `"constant`". 

  * The model `"compressible`" requires `"underformed soil porosity"`" [double],
    `"reference pressure`" [double], and `"pore compressibility`" [string] [Pa^-1].
    Default value for `"reference pressure`" is 101325.0 [Pa].

  * The model `"constant`" requires `"value`" [double].

.. code-block:: xml

  <ParameterList name="flow">  <!-- parent list -->
  <ParameterList name="porosity models">
    <ParameterList name="_SOIL1">
      <Parameter name="regions" type="Array(string)" value="{_TOP HALF}"/>
      <Parameter name="porosity model" type="string" value="constant"/>
      <Parameter name="value" type="double" value="0.2"/>
    </ParameterList>

    <ParameterList name="_SOIL2">
      <Parameter name="regions" type="Array(string)" value="{_BOTTOM HALF}"/>
      <Parameter name="porosity model" type="string" value="compressible"/>
      <Parameter name="underformed soil porosity" type="double" value="0.2"/>
      <Parameter name="reference pressure" type="double" value="101325.0"/>
      <Parameter name="pore compressibility" type="double" value="1e-8"/>
    </ParameterList>
  </ParameterList>
  </ParameterList>

In this example, we define two different porosity models in two soils.

*/

#ifndef AMANZI_FLOW_POROSITY_HH_
#define AMANZI_FLOW_POROSITY_HH_

#include <string>

namespace Amanzi {
namespace Flow {

class Porosity {
 public:
  virtual ~Porosity(){};
  virtual double PorosityValue(double p) = 0;
  virtual double dPorositydPressure(double p) = 0; // derivative wrt to pressure
};

} // namespace Flow
} // namespace Amanzi

#endif
