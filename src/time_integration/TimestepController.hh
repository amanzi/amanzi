/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*!

A timestep controller is used to abstract across PKs algorithms for computing a
*physical* timestep size that is valid for that PK.  Several controller types
are available which implement common choices.

The timestep controller *fixed* provides a single, uniform dt.

The timestep controller *standard* is a simple timestep control mechanism
which sets the next timestep based upon the previous timestep and how many
nonlinear iterations the previous timestep took to converge.
The next timestep is given by the following rule:

- if :math:`N_k > N^{max}` then :math:`\Delta t_{k+1} = f_{reduction} \Delta t_{k}`
- if :math:`N_k < N^{min}` then :math:`\Delta t_{k+1} = f_{increase} \Delta t_{k}`
- otherwise :math:`\Delta t_{k+1} = \Delta t_{k}`

where :math:`\Delta t_{k}` is the previous timestep and :math:`N_k` is the number of nonlinear
iterations required to solve step :math:`k`.

The timestep controller *smart* is based on *standard*, but also tries to be a bit
smarter to avoid repeated increase/decrease loops where the step size decreases,
converges in few iterations, increases, but then fails again.  It also tries to grow
the timestep geometrically to more quickly recover from tricky nonlinearities.

The timestep controller *from file* loads a timestep history from a file, then
advances the step size with those values.  This is mostly used for testing
purposes, where we need to force the same timestep history as previous runs to
do regression testing.  Otherwise roundoff errors can eventually alter
number of iterations enough to alter the timestep history, resulting in
solutions which are enough different to cause doubt over their correctness.

Note, the timestep controller is also responsible for throwing an error if the
timestep becomes invalid.

*/

#pragma once

#include "Teuchos_RCP.hpp"

namespace Amanzi {

class State;

class TimestepController {
 public:
  // virtual destructor
  virtual ~TimestepController(){};

  // single method for timestep control
  //
  // Note, iterations < 0 implies that the previous step failed to converge,
  // while valid implies that both the step converged and the physics deemed
  // the step valid.
  virtual double getTimestep(double dt_prev,
                             int iterations,
                             bool valid) = 0;

  // default method for initial timestep size
  virtual double getInitialTimestep() {
    return getTimestep(-1.0, 0, true);
  }
};

} // namespace Amanzi
