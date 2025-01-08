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

.. _timestep_controller_fixed-spec:
.. admonition:: timestep_controller_fixed-spec

   * `"initial timestep [s]`" ``[double]`` The fixed timestep size.

*/

#pragma once

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "errors.hh"
#include "State.hh"
#include "TimestepController.hh"

namespace Amanzi {

class TimestepControllerFixed : public TimestepController {
 public:
  TimestepControllerFixed(Teuchos::ParameterList& plist);

  // single method for timestep control
  double getTimestep(double dt, int iterations, bool valid) override;

 protected:
  double dt_;
};

} // namespace Amanzi

