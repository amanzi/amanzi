/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon
*/

//!  Slightly smarter timestep controller based upon a history of previous timesteps.

/*!

This is based on `Timestep Controller Standard`_, but also tries to be a bit
smarter to avoid repeated increase/decrease loops where the step size
decreases, converges in few iterations, increases, but then fails again.  It
also tries to grow the step geometrically to more quickly recover from tricky
nonlinearities.

.. _timestep-controller-smarter-spec:
.. admonition:: timestep-controller-smarter-spec

    * `"max iterations`" ``[int]`` :math:`N^{max}`, decrease the timestep if the previous step took more than this.
    * `"min iterations`" ``[int]`` :math:`N^{min}`, increase the timestep if the previous step took less than this.
    * `"time step reduction factor`" ``[double]`` :math:`f_{reduction}`, reduce the previous timestep by this multiple.
    * `"time step increase factor`" ``[double]`` :math:`f_{increase}`, increase the previous timestep by this multiple.  Note that this can be modified geometrically in the case of repeated successful steps.
    * `"max time step increase factor`" ``[double]`` **10.** The max :math:`f_{increase}` will ever get.
    * `"growth wait after fail`" ``[int]`` Wait at least this many timesteps before attempting to grow the timestep after a failed timestep.
    * `"count before increasing increase factor`" ``[int]`` Require this many successive increasions before multiplying :math:`f_{increase}` by itself.


*/


#ifndef AMANZI_SMARTER_TIMESTEP_CONTROLLER_HH_
#define AMANZI_SMARTER_TIMESTEP_CONTROLLER_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "TimestepController.hh"

namespace Amanzi {

class State;

class TimestepControllerSmarter : public TimestepController {

 public:
  TimestepControllerSmarter(Teuchos::ParameterList& plist,
                            const Teuchos::RCP<State>& S = Teuchos::null);

  // single method for timestep control
  double get_timestep(double dt, int iterations);

 protected:
  Teuchos::ParameterList plist_;
  std::string name_;

  int max_its_;
  int min_its_;

  double reduction_factor_;

  double increase_factor0_;
  double max_increase_factor_;
  int count_increased_before_increase_;

  double max_dt_;
  double min_dt_;

  int growth_wait_after_fail0_;

  // State variables stored in RCP to keep in state for checkpointing
  Teuchos::RCP<State> S_;
  Teuchos::RCP<double> increase_factor_;
  Teuchos::RCP<int> successive_increases_;
  Teuchos::RCP<int> last_fail_;
  Teuchos::RCP<int> growth_wait_after_fail_;
};

} // namespace

#endif
