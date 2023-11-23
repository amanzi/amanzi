/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*!

The Barton-Bandis aperture-stress relationship is defibed by initial aperture
:math:`a_0`, overbuden pressure/normal stress :math:`p_{ov}`, fluid pressure 
:math:`p_f`, and two fitting parameters A and B:


.. math::
  a = a_0 - \frac{A (p_{ov} - p_f)}{1 + B (p_{ov} - p_f)}

* `"undeformed aperture`" [double] aperture at zero effective normal stress 

* `"overburden pressure`" [double] overburden pressure/total normal stress

* `"BartonBandis A`" [double] fitting parameter, :math:`[m / Pa]`
* `"BartonBandis B`" [double] fitting parameter, :math:`[Pa^{-1}]`

.. code-block:: xml

  <ParameterList name="fracture aperture models"> <!-- parent list -->
    <ParameterList name="FAM 0">
       <Parameter name="fracture aperture model" type="string" value="Barton Bandis"/>
       <Parameter name="regions" type="Array(string)" value="{RegionBottom}"/>
       <Parameter name="undeformed aperture" type="double" value="1e-5"/>
       <Parameter name="overburden pressure" type="double" value="1e+10"/>
       <Parameter name="BartonBandis A" type="double" value="2.22e-5"/>
       <Parameter name="BartonBandis B" type="double" value="3.47e-2"/>
    </ParameterList>

*/

#ifndef AMANZI_FLOW_APERTURE_MODEL_BARTON_BANDIS_HH_
#define AMANZI_FLOW_APERTURE_MODEL_BARTON_BANDIS_HH_

#include <cmath>

#include "Teuchos_ParameterList.hpp"

#include "ApertureModel.hh"

namespace Amanzi {
namespace Flow {

class ApertureModel_BartonBandis : public ApertureModel {
 public:
  explicit ApertureModel_BartonBandis(Teuchos::ParameterList& plist)
  {
    a0_ = plist.get<double>("undeformed aperture");
    pov_ = plist.get<double>("overburden pressure");
    A_ = plist.get<double>("BartonBandis A");
    B_ = plist.get<double>("BartonBandis B");
  }
  ~ApertureModel_BartonBandis(){};

  // required methods from the base class
  double Aperture(double p) { return a0_ - A_ * (pov_ - p) / (1.0 + B_ * (pov_ - p)); }

  double dAperturedPressure(double p) { return A_ / std::pow(1.0 + B_ * (pov_ - p), 2); }

 private:
  double a0_, pov_, A_, B_;
};

} // namespace Flow
} // namespace Amanzi

#endif
