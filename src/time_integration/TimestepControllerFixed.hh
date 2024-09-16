/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*!

A fixed timestep controller simply sets a constant timestep size.

`"timestep controller type`" = `"fixed`"

.. _timestep-controller-fixed-spec:
.. admonition:: timestep-controller-fixed-spec

   * `"initial timestep [s]`" ``[double]`` The fixed timestep size.

*/

#pragma once

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "errors.hh"
#include "State.hh"
#include "TimestepController.hh"

namespace Amanzi {

template<typename Vector>
class TimestepControllerFixed : public TimestepController<Vector> {
 public:
  TimestepControllerFixed(Teuchos::ParameterList& plist);

  // single method for timestep control
  double getTimestep(double dt, int iterations, bool valid) override;

 protected:
  double dt_;
};


template<typename Vector>
TimestepControllerFixed<Vector>::TimestepControllerFixed(Teuchos::ParameterList& plist)
  : TimestepController<Vector>()
{
  if (plist.isParameter("initial timestep [s]")) {
    Errors::Message msg("Deprecated parameter \"initial timestep [s]\" provided to TimestepControllerFixed, use \"timestep [s]\" instead.");
    Exceptions::amanzi_throw(msg);
  }
  if (plist.isParameter("timestep [s]")) {
    dt_ = plist.get<double>("timestep [s]");
  } else {
    dt_ = plist.get<double>("timestep");
  }
}


// single method for timestep control
template<typename Vector>
double
TimestepControllerFixed<Vector>::getTimestep(double dt, int iterations, bool valid)
{
  if (dt > 0 && (iterations < 0 || !valid)) {
    // note that dt < 0 --> initial step, where iterations < 0 is valid.
    Errors::TimestepCrash msg("Timestep failed: fixed timestep size failed.");
    Exceptions::amanzi_throw(msg);
  }
  return dt_;
}

} // namespace Amanzi


