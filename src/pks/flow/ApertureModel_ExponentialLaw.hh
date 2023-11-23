/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*!

The expoential aperture-stress relationship is defibed by initial aperture
:math:`a_0`, overbuden pressure/normal stress :math:`p_{ov}`, fluid pressure 
:math:`p_f`, and fracture compressibility (fitting parameter) :math:`alpha`:


.. math::
  a = a_0 \exp{\alpha (p_{ov} - p_f)}

* `"undeformed aperture`" [double] aperture at zero effective normal stress 

* `"overburden pressure`" [double] overburden pressure/total normal stress

* `"compressibility`" [double] fracture compressibility, :math:`[Pa^{-1}]`

.. code-block:: xml

  <ParameterList name="fracture aperture models"> <!-- parent list -->
    <ParameterList name="FAM 0">
       <Parameter name="fracture aperture model" type="string" value="exponential law"/>
       <Parameter name="regions" type="Array(string)" value="{RegionBottom}"/>
       <Parameter name="undeformed aperture" type="double" value="1e-5"/>
       <Parameter name="overburden pressure" type="double" value="1e+10"/>
       <Parameter name="compressibility" type="double" value="1.0e-11"/>
    </ParameterList>

*/

#ifndef AMANZI_FLOW_APERTURE_MODEL_EXPONENTIAL_LAW_HH_
#define AMANZI_FLOW_APERTURE_MODEL_EXPONENTIAL_LAW_HH_

#include <cmath>

#include "Teuchos_ParameterList.hpp"

#include "ApertureModel.hh"

namespace Amanzi {
namespace Flow {

class ApertureModel_ExponentialLaw : public ApertureModel {
 public:
  explicit ApertureModel_ExponentialLaw(Teuchos::ParameterList& plist)
  {
    a0_ = plist.get<double>("undeformed aperture");
    pov_ = plist.get<double>("overburden pressure");
    alpha_ = plist.get<double>("compressibility");
  }
  ~ApertureModel_ExponentialLaw(){};

  // required methods from the base class
  double Aperture(double p) { return a0_ * std::exp(-alpha_ * (pov_ - p)); }

  double dAperturedPressure(double p) { return a0_ * alpha_ * std::exp(-alpha_ * (pov_ - p)); }

 private:
  double a0_, pov_, alpha_;
};

} // namespace Flow
} // namespace Amanzi

#endif
