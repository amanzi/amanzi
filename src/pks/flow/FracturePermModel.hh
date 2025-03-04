/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*!

A list of permeability models in a fracture network contains similar sublists
that must cover all network.
Each sublist has two paremeters.

* `"region`" [string] defines region where model applies.

* `"model`" [string] specifies the model name. Currently only one parameter
  is available, `"cubic law`".

.. code-block:: xml

  <ParameterList name="flow">  <!-- parent list -->
  <ParameterList name="fracture permeability models">
    <ParameterList name="_ONE FRACTURE LEAVE">
      <Parameter name="region" type="string" value="fracture"/>
      <Parameter name="model" type="string" value="cubic law"/>
    </ParameterList>
  </ParameterList>

*/

#ifndef AMANZI_FLOW_FRACTURE_PERM_MODEL_HH_
#define AMANZI_FLOW_FRACTURE_PERM_MODEL_HH_

#include <string>

namespace Amanzi {
namespace Flow {

class FracturePermModel {
 public:
  virtual ~FracturePermModel(){};
  virtual double Permeability(double aperture) = 0;
};

} // namespace Flow
} // namespace Amanzi

#endif
