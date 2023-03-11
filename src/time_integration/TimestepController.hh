/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*!

The time step is controlled by parameter *time step controller type*
and the related list of options.
Nonlinear solver is controlled by parameter *solver type*  and related list of options.
Amanzi supports a few nonlinear solvers described in details in a separate section.

The time step controller *standard* is a simple timestep control mechanism
which sets the next timestep based upon the previous timestep and how many
nonlinear iterations the previous timestep took to converge.
The next time step is given by the following rule:

* if :math:`N_k > N^{max}` then :math:`\Delta t_{k+1} = f_{reduction} \Delta t_{k}`
* if :math:`N_k < N^{min}` then :math:`\Delta t_{k+1} = f_{increase} \Delta t_{k}`
* otherwise :math:`\Delta t_{k+1} = \Delta t_{k}`

where :math:`\Delta t_{k}` is the previous timestep and :math:`N_k` is the number of nonlinear 
iterations required to solve step :math:`k`.

The time step controller *smart* is based on *standard*, but also tries to be a bit 
smarter to avoid repeated increase/decrease loops where the step size decreases, 
converges in few iterations, increases, but then fails again.  It also tries to grow 
the time step geometrically to more quickly recover from tricky nonlinearities.

The time step controller *from file* loads a timestep history from a file, then
advances the step size with those values.  This is mostly used for testing
purposes, where we need to force the same timestep history as previous runs to
do regression testing.  Otherwise roundoff errors can eventually alter
number of iterations enough to alter the timestep history, resulting in
solutions which are enough different to cause doubt over their correctness.

*/

#ifndef AMANZI_TIMESTEP_CONTROLLER_HH_
#define AMANZI_TIMESTEP_CONTROLLER_HH_

#include "Teuchos_RCP.hpp"

namespace Amanzi {

class TimestepController {
 public:
  TimestepController(Teuchos::ParameterList& plist)
    : dt_init_(plist.get<double>("initial time step [s]", 1.0))
  {}

  // virtual destructor
  virtual ~TimestepController(){};

  // single method for timestep control
  virtual double get_timestep(double dt, int iterations) = 0;

  // default method for initial timestep size
  virtual double get_initial_timestep() { return dt_init_; }

 protected:
  double dt_init_;
};

} // namespace Amanzi

#endif
