/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

//! Simple timestep control based upon previous iteration count.
/*!

This is a simple timestep control mechanism
which sets the next timestep based upon the previous timestep and how many
nonlinear iterations the previous timestep took to converge.

The timestep for step :math:`k+1`, :math:`\Delta t_{k+1}`, is given by:

- if :math:`N_k > N^{max}` then :math:`\Delta t_{k+1} = f_{reduction} * \Delta t_{k}`
- if :math:`N_k < N^{min}` then :math:`\Delta t_{k+1} = f_{increase} * \Delta t_{k}`
- otherwise :math:`\Delta t_{k+1} = \Delta t_{k}`

where :math:`\Delta t_{k}` is the previous timestep and :math:`N_k` is the number of
nonlinear iterations required to solve step :math:`k`:.

.. _timestep_controller_standard-spec:
.. admonition:: timestep_controller_standard-spec

   * `"max iterations`" ``[int]`` :math:`N^{max}`, decrease the timestep if the previous step took more than this.
   * `"min iterations`" ``[int]`` :math:`N^{min}`, increase the timestep if the previous step took less than this.
   * `"timestep reduction factor`" ``[double]`` :math:`f_{reduction}`, reduce the previous timestep by this multiple.
   * `"timestep increase factor`" ``[double]`` :math:`f_{increase}`, increase the previous timestep by this multiple.

   INCLUDES:

   - ``[timestep_controller_recoverable-spec]``

.. code-block:: xml

  <ParameterList name="BDF1"> <!-- parent list -->
    <Parameter name="timestep controller type" type="string" value="standard"/>
    <ParameterList name="timestep controller standard parameters">
      <Parameter name="min iterations" type="int" value="10"/>
      <Parameter name="max iterations" type="int" value="15"/>
      <Parameter name="timestep increase factor" type="double" value="1.2"/>
      <Parameter name="timestep reduction factor" type="double" value="0.5"/>
      <Parameter name="max timestep [s]" type="double" value="1e+9"/>
      <Parameter name="min timestep [s]" type="double" value="0.0"/>
      <Parameter name="initial timestep [s]" type="double" value="86400"/>
    </ParameterList>
  </ParameterList>

In this example, the timestep is increased by factor 1.2 when the nonlinear
solver converges in 10 or less iterations.
The timestep is not changed when the number of nonlinear iterations is
between 11 and 15.
The timestep will be cut twice if the number of nonlinear iterations exceeds 15.

*/

#ifndef AMANZI_STANDARD_TIMESTEP_CONTROLLER_HH_
#define AMANZI_STANDARD_TIMESTEP_CONTROLLER_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "TimestepControllerRecoverable.hh"

namespace Amanzi {

class TimestepControllerStandard : public TimestepControllerRecoverable {
 public:
  TimestepControllerStandard(const std::string& name,
                             Teuchos::ParameterList& plist,
                             const Teuchos::RCP<State>& S);

 protected:
  // single method for timestep control
  double getTimestep_(double dt, int iterations, bool valid) override;

 protected:
  int max_its_;
  int min_its_;
  double reduction_factor_;
  double increase_factor_;
};

} // namespace Amanzi

#endif
