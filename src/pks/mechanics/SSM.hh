/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*!

User defines small strain models in sublist *SSM models*.
It contains as many sublists, e.g. _SOIL1 and _SOIL2, as there are different soils.
The small strain models are associated with non-overlapping regions.
Each of the sublists (e.g. _SOIL1) includes a few mandatory parameters: *regions names*,
*model name*, and parameters for the selected model.

* `"model`" [string] specifies a model for the soil.
  The available models are `"Hardin Drnevich`" and `"unknown`".

  * The model `"Hardin Drnevich`" requires mandatory parameters:
    * `"reference shear strain`" [double] [-]
    * `"maximum shear stress`" [double] [Pa]

.. code-block:: xml

  <ParameterList name="mechanics">  <!-- parent list -->
  <ParameterList name="small strain models">
    <ParameterList name="_SOIL1">
      <Parameter name="regions" type="Array(string)" value="{_CLAY}"/>
      <Parameter name="model" type="string" value="Hardin Drnevich"/>
      <Parameter name="reference shear strain" type="double" value="8.0e-4"/>
      <Parameter name="maximum shear stress" type="double" value="50.0e+6"/>
    </ParameterList>
  </ParameterList>
  </ParameterList>

*/

#ifndef AMANZI_MECHANICS_SMALL_STRAIN_HH_
#define AMANZI_MECHANICS_SMALL_STRAIN_HH_

#include <string>

namespace Amanzi {
namespace Mechanics {

class SSM {
 public:
  virtual ~SSM() {};
  // gamma - shear strain, e - volumetric strain
  virtual double ShearStress(double gamma) = 0;
  virtual double BulkModulus(double e) = 0;
};

} // namespace Mechanics
} // namespace Amanzi

#endif
