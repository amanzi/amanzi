/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*!

Constant aperture
* `"value`" [double] constant aperture

.. code-block:: xml

  <ParameterList name="fracture aperture models"> <!-- parent list -->
    <ParameterList name="FAM 0">
       <Parameter name="fracture aperture model" type="string" value="constant"/>
       <Parameter name="regions" type="Array(string)" value="{RegionBottom}"/>
       <Parameter name="value" type="double" value="1e-5"/>
    </ParameterList>

*/

#ifndef AMANZI_FLOW_APERTURE_MODEL_CONSTANT_HH_
#define AMANZI_FLOW_APERTURE_MODEL_CONSTANT_HH_

#include <cmath>

#include "Teuchos_ParameterList.hpp"

#include "ApertureModel.hh"

namespace Amanzi {
namespace Flow {

class ApertureModel_Constant : public ApertureModel {
 public:
  explicit ApertureModel_Constant(Teuchos::ParameterList& plist)
  {
    a0_ = plist.get<double>("value");
  }
  ~ApertureModel_Constant(){};

  // required methods from the base class
  double Aperture(double p) { return a0_; }

  double dAperturedPressure(double p) { return 0.0; }

 private:
  double a0_;
};

} // namespace Flow
} // namespace Amanzi

#endif
