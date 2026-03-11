/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*!

.. admonition:: flow_source-spec

  * `"regions`" ``[Array(string)]`` is the list of regions where the source is defined.

  * `"spatial distribution method`" ``[string]`` is the method for distributing
    source Q over the specified regions. The available options are `"volume`" and `"none`".
    For option `"none`", the source term function Q is measured in [J/m^3/s].
    For the other options, it is measured in [J/s].

  * `"submodel`" ``[string]`` refines definition of the source. Available options are

    * `"conductive heat supply`" defines time-dependent supply of heat from an infinite 
      energy reservoir. It requires `"thermal conductivity"` [W/m/K] and 
      `"thermal diffusivity`" [m^2/s] for this reservoir.

.. code-block:: xml

  <ParameterList name="energy">  <!-- parent list -->
  <ParameterList name="source terms">
    <ParameterList name="SRC 0">
      <Parameter name="regions" type="Array(string)" value="{{EntireDomain}}"/>
      <Parameter name="spatial distribution method" type="string" value="node"/>
      <Parameter name="submodel" type="string" value="conductive heat supply"/>
      <Parameter name="thermal conductivity" type="double" value="2.1"/>
      <Parameter name="thermal diffusivity" type="double" value="8.4e-7"/>
      <ParameterList name="temperature at infinity">
        <ParameterList name="function-constant">
          <Parameter name="value" type="double" value="500.0"/>
        </ParameterList>
      </ParameterList>
    </ParameterList>
  </ParameterList>

*/

#ifndef AMANZI_ENERGY_SOURCE_FUNCTION_HH_
#define AMANZI_ENERGY_SOURCE_FUNCTION_HH_

#include <string>
#include <vector>

#include "Teuchos_RCP.hpp"

#include "CompositeVector.hh"
#include "PK_DomainFunction.hh"
#include "State.hh"

namespace Amanzi {
namespace Energy {

class EnergySourceFunction : public PK_DomainFunction {
 public:
  EnergySourceFunction() {};
  EnergySourceFunction(const Teuchos::ParameterList& plist) {
    if (plist.isParameter("submodel")) submodel_ = plist.get<std::string>("submodel");
    name_ = submodel_;
  };

  void ComputeSubmodel(double t_old, double t_new) {};

 private:
  std::string submodel_;
};

} // namespace Energy
} // namespace Amanzi

#endif
