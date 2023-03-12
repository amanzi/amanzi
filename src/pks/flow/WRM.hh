/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Neil Carlson
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*!

User defines water retention models in sublist *water retention models*. 
It contains as many sublists, e.g. *SOIL_1*, *SOIL_2*, etc, as there are different soils. 
This list is required for the Richards problem only.
 
The water retention models are associated with non-overlapping regions. Each of the sublists (e.g. *Soil 1*) 
includes a few mandatory parameters: region name, model name, and parameters for the selected model.

* `"water retention model`" [string] specifies a model for the soil.
  The available models are `"van Genuchten`", `"Brooks Corey`", and `"fake`". 
  The later is used only to set up a simple analytic solution for convergence study. 

  * The model `"van Genuchten`" requires `"van Genuchten alpha`" [double],
    `"van Genuchten m`" [double], `"van Genuchten l`" [double], `"residual saturation`" [double],
    and `"relative permeability model`" [string].

  * The model `"Brooks-Corey`" requires `"Brooks Corey lambda`" [double], `"Brooks Corey alpha`" [double],
    `"Brooks Corey l`" [double], `"residual saturation`" [double],
    and `"relative permeability model`" [string].

* `"relative permeability model`" [string] The available options are `"Mualem`" (default) 
  and `"Burdine`".

* `"regularization interval`" [double] removes the kink in the water retention curve at the
  saturation point using a cubic spline. The parameter specifies the regularization region [Pa].
  Default value is 0.

Amanzi performs rudimentary checks of validity of the provided parameters. 
The relative permeability curves can be calculated and saved in an ASCI file 
if the list *output* is provided. This list has two mandatory parameters:

* `"file`" [string] is the user defined file name. It should be different for 
  each soil. 

* `"number of points`" [int] is the number of data points. 
  Each file will contain a table with three columns: saturation, relative permeability, and
  capillary pressure. The data points are equidistributed between the residual saturation
  and 1.

.. code-block:: xml

  <ParameterList name="flow">  <!-- parent list -->
  <ParameterList name="water retention models">
    <ParameterList name="_SOIL_1">
      <Parameter name="regions" type="Array(string)" value="{_TOP HALF}"/>
      <Parameter name="water retention model" type="string" value="van Genuchten"/>
      <Parameter name="van Genuchten alpha" type="double" value="0.000194"/>
      <Parameter name="van Genuchten m" type="double" value="0.28571"/>
      <Parameter name="van Genuchten l" type="double" value="0.5"/>
      <Parameter name="residual saturation" type="double" value="0.103"/>
      <Parameter name="regularization interval" type="double" value="100.0"/>
      <Parameter name="relative permeability model" type="string" value="Mualem"/>
      <ParameterList name="output">
        <Parameter name="file" type="string" value="soil1.txt"/>
        <Parameter name="number of points" type="int" value="1000"/>
      </ParameterList>
    </ParameterList>

    <ParameterList name="_SOIL_2">
      <Parameter name="regions" type="Array(string)" value="{_BOTTOM HALF}"/>
      <Parameter name="water retention model" type="string" value="Brooks Corey"/>
      <Parameter name="Brooks Corey lambda" type="double" value="0.0014"/>
      <Parameter name="Brooks Corey alpha" type="double" value="0.000194"/>
      <Parameter name="Brooks Corey l" type="double" value="0.51"/>
      <Parameter name="residual saturation" type="double" value="0.103"/>
      <Parameter name="regularization interval" type="double" value="0.0"/>
      <Parameter name="relative permeability model" type="string" value="Burdine"/>
    </ParameterList>
  </ParameterList>
  </ParameterList>

In this example, we define two different water retention models in two soils.

*/

#ifndef AMANZI_WATER_RETENTION_MODEL_HH_
#define AMANZI_WATER_RETENTION_MODEL_HH_

#include <string>

namespace Amanzi {
namespace Flow {

class WRM {
 public:
  virtual ~WRM(){};

  virtual double k_relative(double pc) const = 0;
  virtual double saturation(double pc) const = 0;
  virtual double
  dSdPc(double pc) const = 0; // derivative of saturation w.r.t. to capillary pressure
  virtual double capillaryPressure(double s) const = 0;
  virtual double residualSaturation() const = 0;
  virtual double dKdPc(double pc) const = 0;
};

typedef double (WRM::*KRelFn)(double pc) const;

} // namespace Flow
} // namespace Amanzi

#endif
