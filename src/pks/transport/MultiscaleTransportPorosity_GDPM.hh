/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*!

A generalized dual porosity model, fracture + multiple matrix nodes.
Current naming convention is that the fields used in the single-porosity
model correspond now to the fracture continuum.
Example: tcc = total component concentration in the fracture continuum;
tcc_matrix = total component concentration in the matrix continuum.

* `"number of matrix nodes`" [int] defines number of matrix layers.
* `"matrix depth`" [double] is the characteristic length for matrix continuum.
* `"tortousity`" [double] defines tortuosity to correct diffusivity of a liquid solute.
* `"matrix volume fraction`" [double] defines relative volume of matrix continuum.

.. code-block:: xml

  <ParameterList name="_TRANSPORT">  <!-- parent list -->
  <ParameterList name="multiscale models">
    <ParameterList name="_WHITE SOIL">
      <Parameter name="multiscale model" type="string" value="dual porosity"/>
      <Parameter name="regions" type="Array(string)" value="{_TOP_REGION, _BOTTOM_REGION}"/>
      <ParameterList name="dual porosity parameters">
        <Paramater name="Warren Root parameter" type="double" value="4.0e-5"/>
        <Paramater name="matrix tortuosity" type="double" value="0.95"/>
        <Paramater name="matrix volume fraction" type="double" value="0.9999"/>
      </ParameterList>  
    </ParameterList>  

    <ParameterList name="_GREY SOIL">
      <Parameter name="multiscale model" type="string" value="generalized dual porosity"/>
      <Parameter name="regions" type="Array(string)" value="{_MIDDLE_REGION}"/>
      <ParameterList name="generalized dual porosity parameters">
        <Paramater name="number of matrix nodes" type="int" value="2"/>
        <Paramater name="matrix depth" type="double" value="0.01"/>
        <Paramater name="matrix tortuosity" type="double" value="1.0"/>
      </ParameterList>  
    </ParameterList>  
  </ParameterList>  
  </ParameterList>  

*/

#ifndef MULTISCALE_TRANSPORT_POROSITY_GDPM_HH_
#define MULTISCALE_TRANSPORT_POROSITY_GDPM_HH_

#include <string>
#include <vector>

#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"
#include "Mini_Diffusion1D.hh"

#include "MultiscaleTransportPorosity.hh"

namespace Amanzi {
namespace Transport {

class MultiscaleTransportPorosity_GDPM : public MultiscaleTransportPorosity {
 public:
  MultiscaleTransportPorosity_GDPM(Teuchos::ParameterList& plist);
  ~MultiscaleTransportPorosity_GDPM(){};

  // Compute solute flux: icomp - component id, phi - matrix porosity,
  // tcc_m_aux - vector of concentration values in secondary nodes,
  // wfm[0|1] - fracture water content at initial and final time moments,
  // wcm[0|1] - water content at initial and final time moments
  virtual double ComputeSoluteFlux(double flux_liquid,
                                   double& tcc_f,
                                   WhetStone::DenseVector& tcc_m,
                                   int icomp,
                                   double dt,
                                   double wcf0,
                                   double wcf1,
                                   double wcm0,
                                   double wcm1,
                                   double phi) override;

  // Number of matrix nodes
  virtual int NumberMatrixNodes() override { return matrix_nodes_; }

 private:
  static Utils::RegisteredFactory<MultiscaleTransportPorosity, MultiscaleTransportPorosity_GDPM>
    factory_;
  std::vector<Operators::Mini_Diffusion1D> op_diff_;

  int matrix_nodes_;
  double depth_, tau_;
  std::vector<double> mol_diff_;
};

} // namespace Transport
} // namespace Amanzi

#endif
