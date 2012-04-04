/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

   energy-bc-list(NAME) is:
  <ParameterList name=NAME>
    energy-condition-list
    ...
    energy-condition-list
  </ParameterList>

where energy-condition-list is one of:
    energy-temperature-list
    energy-enthalpy-flux-list

The parameter list name string NAME is arbitrary, and meaningful only to the
parent parameter list.  Each sublist defines one type of condition for a
energy problem.  Each specific type of sublist can appear at most once.  This
parameter list is given to a boundary condition "factory" which has methods
for instantiating the appropriate boundary condition objects.

1. energy-temperature-list is:
  <ParameterList name="temperature">
    energy-temperature-spec(NAME_1)
    ...
    energy-temperature-spec(NAME_N)
  </ParameterList>

2. energy-enthalpy-flux-list is:
  <ParameterList name="enthalpy flux">
    energy-enthalpy-flux-spec(NAME_1)
    ...
    energy-enthalpy-flux-spec(NAME_N)
  </ParameterList>

Each spec sublist defines one part of the of total boundary condition of that type.
The name strings NAME_1, ..., NAME_N are arbitrary but must be unique
within the parent parameter list.  They may be used in error messages
or diagnostic logging.  Example: "north cribs"
  Each of the following kinds of spec parameter lists have a "regions"
parameter that is an array of 1 or more region names that specify the
portion of the mesh boundary where the condition is to be applied.

energy-temperature-spec(NAME) is:
  <ParameterList name=NAME>
    <Parameter name="regions" type="Array string" value=string-array />
    function-factory-list("boundary temperature")
  </ParameterList>

The function-factory-list should define a function whose argument
is the vector (t, x, y, z).

energy-enthalpy-flux-spec(NAME) is:
  <ParameterList name=NAME>
    <Parameter name="Regions" type="Array string" value=string-array />
    function-factory-list("outward enthalpy flux")
  </ParameterList>

The function-factory-list should define a function whose argument
is the vector (t, x, y, z).

------------------------------------------------------------------------- */

#ifndef _ENERGY_BC_FACTORY_HH_
#define _ENERGY_BC_FACTORY_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Mesh.hh"
#include "boundary-function.hh"

namespace Amanzi {
namespace Energy {

class EnergyBCFactory {

public:
  EnergyBCFactory(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                  Teuchos::ParameterList& plist) :
    mesh_(mesh), plist_(plist) {}

  Teuchos::RCP<BoundaryFunction> CreateTemperature() const;
  Teuchos::RCP<BoundaryFunction> CreateEnthalpyFlux() const;

private:
  void ProcessTemperatureList_(const Teuchos::ParameterList&,
          const Teuchos::RCP<BoundaryFunction>&) const;
  void ProcessTemperatureSpec_(const Teuchos::ParameterList&,
          const Teuchos::RCP<BoundaryFunction>&) const;
  void ProcessEnthalpyFluxList_(const Teuchos::ParameterList&,
          const Teuchos::RCP<BoundaryFunction>&) const;
  void ProcessEnthalpyFluxSpec_(const Teuchos::ParameterList&,
          const Teuchos::RCP<BoundaryFunction>&) const;

  const Teuchos::RCP<const AmanziMesh::Mesh>& mesh_;
  Teuchos::ParameterList plist_;
};

}  // namespace Energy
}  // namespace Amanzi

#endif
