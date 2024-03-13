/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/*!

The `"isotherms`" section is the list of sorption processes that relate the sorbed concentration at
the solid surface to the aqueous concentration in contact with the solid at constant temperature.
It is a function of the free ion primary species concentrations, not total conentrations.
A sorption isotherm may represent equilibrium or kinetic processes depending on the data used to
fit the isotherm.
Each sublist has two parameters:

* `"model`" [string] specifies the model name. The available options are `"linear", `"langmuir`", and `"freundlich`".

* `"parameters`" [Array(double)] is the list of model parameters. The distribution coefficient
  *K* is in the first position.

.. code-block:: xml

  <ParameterList name="isotherms">
    <ParameterList name="A">
      <Parameter name="model" type="string" value="linear"/>
      <Parameter name="parameters" type="Array(double)" value="{10.0}"/>
    </ParameterList>
    <ParameterList name="B">
      <Parameter name="model" type="string" value="langmuir"/>
      <Parameter name="parameters" type="Array(double)" value="{30.0, 0.1}"/>
    </ParameterList>
    <ParameterList name="C">
      <Parameter name="model" type="string" value="freundlich"/>
      <Parameter name="parameters" type="Array(double)" value="{1.5, 1.25}"/>
    </ParameterList>
  </ParameterList>

A few examples are given below.
Each line has three fields: primary species name, adsorption isotherm model, and parameters.
The number of  parameters and their meaning depends on the model; although the first one
is always the distribution coefficient.

.. code-block:: txt

   Pu_238   linear    461168.4
   U_234    linear    329406.0
   Th_230   linear   1482327.0
   Ra_226   linear     41175.75
   Pb_210   linear   3294060.0
   Tc_99    linear       988.218

NOTE: The parameters provided here are *global*.
The state field isotherm_kd overwrites any global data given here.

*/

#ifndef AMANZI_CHEMISTRY_SORP_ISOTHERM_FACTORY_HH_
#define AMANZI_CHEMISTRY_SORP_ISOTHERM_FACTORY_HH_

#include <memory>
#include <string>

#include "Teuchos_ParameterList.hpp"

#include "Species.hh"

namespace Amanzi {
namespace AmanziChemistry {

class SorptionIsotherm;

class SorptionIsothermFactory {
 public:
  SorptionIsothermFactory(){};
  ~SorptionIsothermFactory(){};

  std::shared_ptr<SorptionIsotherm> Create(const Teuchos::ParameterList& plist);

  int VerifySpeciesName(const std::string& species_name, const std::vector<Species>& species) const;

  static const std::string linear;
  static const std::string langmuir;
  static const std::string freundlich;
};

} // namespace AmanziChemistry
} // namespace Amanzi
#endif
