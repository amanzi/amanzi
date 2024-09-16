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

.. _timestep-controller-standard-spec:
.. admonition:: timestep-controller-standard-spec

   * `"max iterations`" ``[int]`` :math:`N^{max}`, decrease the timestep if the previous step took more than this.
   * `"min iterations`" ``[int]`` :math:`N^{min}`, increase the timestep if the previous step took less than this.
   * `"timestep reduction factor`" ``[double]`` :math:`f_{reduction}`, reduce the previous timestep by this multiple.
   * `"timestep increase factor`" ``[double]`` :math:`f_{increase}`, increase the previous timestep by this multiple.

   INCLUDES:

   - ``[timestep-controller-recoverable-spec]``

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

#include "errors.hh"
#include "dbc.hh"
#include "TimestepControllerRecoverable.hh"

namespace Amanzi {


template <class Vector>
class TimestepControllerStandard : public TimestepControllerRecoverable<Vector> {
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


template <class Vector>
TimestepControllerStandard<Vector>::TimestepControllerStandard(const std::string& name,
                             Teuchos::ParameterList& plist,
                             const Teuchos::RCP<State>& S)
  : TimestepControllerRecoverable<Vector>(name, plist, S)
{
  max_its_ = plist.get<int>("max iterations");
  min_its_ = plist.get<int>("min iterations");
  AMANZI_ASSERT(max_its_ > min_its_);
  AMANZI_ASSERT(min_its_ >= 0);

  reduction_factor_ = plist.get<double>("timestep reduction factor");
  AMANZI_ASSERT(reduction_factor_ >= 0.0);
  AMANZI_ASSERT(reduction_factor_ <= 1.0);

  increase_factor_ = plist.get<double>("timestep increase factor");
  AMANZI_ASSERT(increase_factor_ >= 1.0);
}


template <class Vector>
double
TimestepControllerStandard<Vector>::getTimestep_(double dt, int iterations, bool valid)
{
  double dt_next(dt);
  // iterations < 0 implies failed timestep
  if (iterations < 0 || iterations > max_its_ || !valid) {
    dt_next = dt * reduction_factor_;
  } else if (iterations < min_its_) {
    dt_next = dt * increase_factor_;
  }
  return dt_next;
}


} // namespace Amanzi

#endif
