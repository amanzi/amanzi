/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*!

Internal energy model is function of temperature only. Units are J/{mol/kg}.
Internal energy list has a few parameters that allows us to run this PK
in a variety of regimes, e.g. with or without gas phase.

* `"energy key`" [string] specifies name for the internal energy field.
  The default value is `"energy`".

* `"evaluator type`" [string] changes the evaluator for internal energy.
  Available options are `"generic`" and `"constant liquid density`" (default).

* `"vapor diffusion`" [bool] specifies presence of a gas phase.
  The default value is `"true`".

.. code-block:: xml

  <ParameterList name="_ENERGY">  <!-- parent list -->
  <ParameterList name="energy evaluator">
    <Parameter name="energy key" type="string" value="energy"/>
    <Parameter name="evaluator type" type="string" value="constant liquid density"/>
    <Parameter name="vapor diffusion" type="bool" value="true"/>
    <ParameterList name="verbose object">
      <Parameter name="verbosity level" type="string" value="high"/>
    </ParameterList>
  </ParameterList>
  </ParameterList>

*/

#ifndef AMANZI_ENERGY_IEM_HH_
#define AMANZI_ENERGY_IEM_HH_

// #include "Teuchos_ParameterList.hpp"

namespace Amanzi {
namespace Energy {

class IEM {
 public:
  IEM() : ierr_(0){};
  virtual ~IEM() {}

  // IEM(Teuchos::ParameterList& plist);
  virtual double InternalEnergy(double T, double p) = 0;
  virtual double DInternalEnergyDT(double T, double p) = 0;
  virtual double DInternalEnergyDp(double T, double p) = 0;

  // error messages
  int error_code() { return ierr_; }
  std::string error_msg() { return error_msg_; }

 protected:
  int ierr_;
  std::string error_msg_;
};

} // namespace Energy
} // namespace Amanzi

#endif
