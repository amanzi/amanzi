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

* _SOIL [list] Defines material properties.
  
  * `"region`" [Array(string)] Defines geometric regions for material SOIL.
  * `"model`" [string] Defines dispersivity model, choose exactly one of the following: `"scalar`", `"Bear`",
    `"Burnett-Frind`", or `"Lichtner-Kelkar-Robinson`".
  * `"parameters for MODEL`" [list] where `"MODEL`" is the model name.
    For model `"scalar`", *only* one of the following options must be specified:

      * `"alpha`" [double] defines dispersivity in all directions, [m].
      * `"dispersion coefficient`" [double] defines dispersion coefficient [m^2/s].

    For model `"Bear`", the following options must be specified:

      * `"alpha_l`" [double] defines dispersion in the direction of Darcy velocity, [m].
      * `"alpha_t`" [double] defines dispersion in the orthogonal direction, [m].
    
    For model `"Burnett-Frind`", the following options must be specified:

      * `"alphaL`" [double] defines the longitudinal dispersion in the direction of Darcy velocity, [m].
      * `"alpha_th`" [double] Defines the transverse dispersion in the horizonla direction orthogonal directions, [m].
      * `"alpha_tv`" [double] Defines dispersion in the orthogonal directions, [m].
        When `"alpha_th`" equals to `"alpha_tv`", we obtain dispersion in the direction of the Darcy velocity.
        This and the above parameters must be defined for `"Burnett-Frind`" and `"Lichtner-Kelkar-Robinson`" models.

    For model `"Lichtner-Kelker-Robinson`", the following options must be specified:

      * `"alpha_lh`" [double] defines the longitudinal dispersion in the horizontal direction, [m].
      * `"alpha_lv`" [double] Defines the longitudinal dispersion in the vertical direction, [m].
        When `"alpha_lh`" equals to `"alpha_lv`", we obtain dispersion in the direction of the Darcy velocity.
        This and the above parameters must be defined for `"Burnett-Frind`" and `"Lichtner-Kelker-Robinson`" models.
      * `"alpha_th`" [double] Defines the transverse dispersion in the horizontal direction orthogonal directions, [m].
      * `"alpha_tv" [double] Defines dispersion in the orthogonal directions.
        When `"alpha_th`" equals to `"alpha_tv`", we obtain dispersion in the direction of the Darcy velocity.
        This and the above parameters must be defined for `"Burnett-Frind`" and `"Lichtner-Kelker-Robinson`" models.

  * `"aqueous tortuosity`" [double] Defines tortuosity for calculating diffusivity of liquid solutes, [-].
  * `"gaseous tortuosity`" [double] Defines tortuosity for calculating diffusivity of gas solutes, [-].
 
Three examples are below:

.. code-block:: xml

  <ParameterList name="_TRANSPORT">  <!-- parent list -->
  <ParameterList name="material properties">
    <ParameterList name="_WHITE SOIL">
      <Parameter name="regions" type="Array(string)" value="{_TOP_REGION, _BOTTOM_REGION}"/>
      <Parameter name="model" type="string" value="Bear"/>
      <ParameterList name="parameters for Bear">
        <Parameter name="alpha_l" type="double" value="1e-2"/>
        <Parameter name="alpha_t" type="double" value="1e-5"/>
      <ParameterList>
      <Parameter name="aqueous tortuosity" type="double" value="1.0"/>
      <Parameter name="gaseous tortuosity" type="double" value="1.0"/>       
    </ParameterList>  
     
    <ParameterList name="_GREY SOIL">
      <Parameter name="regions" type="Array(string)" value="{_MIDDLE_REGION}"/>
      <Parameter name="model" type="string" value="Burnett-Frind"/>
      <ParameterList name="parameters for Burnett-Frind">
        <Parameter name="alpha_l" type="double" value="1e-2"/>
        <Parameter name="alpha_th" type="double" value="1e-3"/>
        <Parameter name="alpha_tv" type="double" value="2e-3"/>
      <ParameterList>
      <Parameter name="aqueous tortuosity" type="double" value="0.5"/>
      <Parameter name="gaseous tortuosity" type="double" value="1.0"/>       
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
  virtual WhetStone::Tensor mech_dispersion(const AmanziGeometry::Point& u,
                                            int axi_symmetry,
                                            double wc,
                                            double phi) const = 0;

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
