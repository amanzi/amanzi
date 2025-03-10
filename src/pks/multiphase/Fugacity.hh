/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*!

 Fugacity model. Units of the returned value is Pa.

.. code-block:: xml

* `"eos type`" [string] specifies model for computing fugacity
* `"component`" [string] is reserved
* `"Henry constant`" [double] used by model `"Henry law`".

  <ParameterList name="fugacity">
    <Parameter name="eos type" type="string" value="ideal gas"/>
    <Parameter name="component" type="string" value="H2"/>
    <Parameter name="Henry constant" type="double" value="1.0"/>
  </ParameterList>

*/

#ifndef AMANZI_FUGACITY_HH_
#define AMANZI_FUGACITY_HH_

namespace Amanzi {
namespace Multiphase {

class Fugacity {
 public:
  Fugacity(){};
  virtual ~Fugacity() {};

  virtual double Value(double T) = 0;
};

} // namespace Multiphase
} // namespace Amanzi

#endif
