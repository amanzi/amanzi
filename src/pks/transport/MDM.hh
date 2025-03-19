/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*!

The material properties include dispersivity model and diffusion parameters
for aqueous and gaseous phases.
The dispersivity is defined as a soil property.
The diffusivity is defined independently for each solute.


.. _mdm-typed-spec:
.. admonition:: mdm-typed-spec

   * `"MDM type`" ``[string]`` Defines dispersivity model, choose exactly one of the following:

     - `"scalar`"
     - `"Bear`",
     - `"Burnett-Frind`"
     - `"Lichtner-Kelkar-Robinson`"
     
   * `"region`" ``[Array(string)]`` Defines geometric regions for material SOIL.

   * `"_MDM_type_ parameters`" ``[_MDM_type_-spec]``] See below for the
     required parameter spec for each type.

Three examples are below:

.. code-block:: xml

  <ParameterList name="_TRANSPORT">  <!-- parent list -->
  <ParameterList name="material properties">
    <ParameterList name="_WHITE SOIL">
      <Parameter name="regions" type="Array(string)" value="{_TOP_REGION, _BOTTOM_REGION}"/>
      <Parameter name="model" type="string" value="Bear"/>
      <ParameterList name="Bear parameters">
        <Parameter name="alpha_l" type="double" value="1e-2"/>
        <Parameter name="alpha_t" type="double" value="1e-5"/>
      <ParameterList>
    </ParameterList>

    <ParameterList name="_GREY SOIL">
      <Parameter name="regions" type="Array(string)" value="{_MIDDLE_REGION}"/>
      <Parameter name="model" type="string" value="Burnett-Frind"/>
      <ParameterList name="Burnett-Frind parameters">
        <Parameter name="alpha_l" type="double" value="1e-2"/>
        <Parameter name="alpha_th" type="double" value="1e-3"/>
        <Parameter name="alpha_tv" type="double" value="2e-3"/>
      <ParameterList>
    </ParameterList>
  </ParameterList>
  </ParameterList>

*/

#ifndef AMANZI_MECHANICAL_DISPERSION_MODEL_HH_
#define AMANZI_MECHANICAL_DISPERSION_MODEL_HH_

// Amanzi
#include "Point.hh"
#include "Tensor.hh"

namespace Amanzi {
namespace Transport {

class MDM {
 public:
  virtual ~MDM(){};

  // returns dispersion tensor.
  virtual WhetStone::Tensor mech_dispersion(double t,
                                            const AmanziGeometry::Point& xc,
                                            const AmanziGeometry::Point& u,
                                            int axi_symmetry,
                                            double wc,
                                            double phi) const = 0;

  // -- compatibility version (t, x) are not used
  WhetStone::Tensor
  mech_dispersion(const AmanziGeometry::Point& u, int axi_symmetry, double wc, double phi)
  {
    return mech_dispersion(0.0, u, u, axi_symmetry, wc, phi);
  }

  // The model is valid if at least one parameter is not zero.
  virtual bool is_valid() const = 0;

  // This allows us to set space dimension which could be used for estimating
  // model applicability.
  virtual void set_dim(int dim) { dim_ = dim; }

 protected:
  int dim_;
};

} // namespace Transport
} // namespace Amanzi

#endif
