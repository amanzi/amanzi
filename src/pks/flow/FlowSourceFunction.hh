/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*!

The sources and sinks are typically associated with wells.
Negative source means a producing well.
Positive source means an injecting well.
The structure of list *source terms* mimics that of list *boundary conditions*.
Again, constant functions can be replaced by any of the available functions.

.. admonition:: flow_source-spec

  * `"regions`" ``[Array(string)]`` is the list of regions where the source is defined.

  * `"spatial distribution method`" ``[string]`` is the method for distributing
    source Q over the specified regions. The available options are `"volume`",
    `"none`", `"permeability`" and `"simple well`".
    For option `"none`", the source term function Q is measured in [kg/m^3/s].
    For the other options, it is measured in [kg/s].
    When the source function is defined over a few regions, Q is distributed over their union.
    Option `"volume fraction`" can be used when the region geometric
    model support volume fractions. Option `"simple well`" implements the Peaceman model.
    The well flux is defined as :math:`q_w = W (p - p_w)` [kg/s], where :math:`W` is
    the well index and :math:`p_w` is the well pressure. The pressure in a well is assumed
    to be hydrostatic.

  * `"use volume fractions`" instructs the code to use all available volume fractions.
    Note that the region geometric model supports volume fractions only for a few regions.

  * `"submodel`" ``[string]`` refines definition of the source. Available options are

    * `"rate`" defines the source in a natural way as the rate of change, `Q`.
      It requires a rate function.

    * `"integrated source`" defines the indefinite integral, `I`, of the rate of change,
      i.e. the source term is calculated as `Q = dI/dt`.
      It requires a function of `I`.

    * `"bhp"` stands for the bottom hole pressure. It requires `"depth"`, `"well radius`"
      and `"bhp`" function.

    Default submodel is `"rate"`.

    For most distributions methods, two submodels are supported: `"rate`" and
    `"integrated source`". For distribution method `"simple well`", the following two
    submodels are supported: `"rate`" and `"bhp`".

.. code-block:: xml

  <ParameterList name="flow">  <!-- parent list -->
  <ParameterList name="source terms">
    <ParameterList name="_SRC 0">
      <Parameter name="regions" type="Array(string)" value="{{_WELL_EAST}}"/>
      <Parameter name="spatial distribution method" type="string" value="volume"/>
      <Parameter name="submodel" type="string" value="rate"/>
      <ParameterList name="well">
        <ParameterList name="function-constant">
          <Parameter name="value" type="double" value="-0.1"/>
        </ParameterList>
      </ParameterList>
    </ParameterList>

    <ParameterList name="_SRC 1">
      <Parameter name="regions" type="Array(string)" value="{{_WELL_WEST}}"/>
      <Parameter name="spatial distribution method" type="string" value="permeability"/>
      <ParameterList name="well">
        <ParameterList name="function-constant">
          <Parameter name="value" type="double" value="-0.2"/>
        </ParameterList>
      </ParameterList>
    </ParameterList>

    <ParameterList name="_SRC 2">
      <Parameter name="regions" type="Array(string)" value="{{_WELL_NORTH}}"/>
      <Parameter name="spatial distribution method" type="string" value="simple well"/>
      <ParameterList name="well">
        <Parameter name="submodel" type="string" value="bhp"/>
        <Parameter name="depth" type="double" value="-2.5"/>
        <Parameter name="well radius" type="double" value="0.1"/>
        <ParameterList name="bhp">
          <ParameterList name="function-constant">
            <Parameter name="value" type="double" value="10.0"/>
          </ParameterList>
        </ParameterList>
      </ParameterList>
    </ParameterList>

    <ParameterList name="_SRC 3">
      <Parameter name="regions" type="Array(string)" value="{{_WELL_SOUTH}}"/>
        <Parameter name="spatial distribution method" type="string" value="simple well"/>
        <ParameterList name="well">
          <Parameter name="submodel" type="string" value="rate"/>
          <ParameterList name="rate">
            <ParameterList name="function-constant">
              <Parameter name="value" type="double" value="100.0"/>
            </ParameterList>
          </ParameterList>
        </ParameterList>
      </ParameterList>
    </ParameterList>
  </ParameterList>
  </ParameterList>

*/

#ifndef AMANZI_FLOW_SOURCE_FUNCTION_HH_
#define AMANZI_FLOW_SOURCE_FUNCTION_HH_

#include <string>
#include <vector>

#include "Teuchos_RCP.hpp"

#include "CompositeVector.hh"
#include "PK_DomainFunction.hh"
#include "State.hh"

namespace Amanzi {
namespace Flow {

class FlowSourceFunction : public PK_DomainFunction {
 public:
  FlowSourceFunction(){};
  FlowSourceFunction(const Teuchos::ParameterList& plist){};

  void ComputeSubmodel(const Key& key, const State& S)
  {
    if (getName() != "volume" && S.HasRecord(key, Tags::DEFAULT)) {
      auto aperture = *S.Get<CompositeVector>(key, Tags::DEFAULT).ViewComponent("cell", true);
      for (auto it = begin(); it != end(); ++it) { it->second[0] *= aperture[0][it->first]; }
    }
  }
};

} // namespace Flow
} // namespace Amanzi

#endif
